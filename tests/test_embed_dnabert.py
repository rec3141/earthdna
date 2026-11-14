#!/usr/bin/env python3
"""Regression tests for pipeline/embed_dnabert.py (code-05 F1).

Exercises pipeline/embed_dnabert.py via subprocess with pure-NumPy stubs for torch and
transformers so no GPU, PyTorch environment, or network access is required. Verifies
desired error/isolation behavior on malformed inputs and exact FSB1 output equivalence
on valid plain and gzipped FASTA controls.
"""
import gzip
import os
import struct
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "pipeline" / "embed_dnabert.py"

STUB_TORCH = textwrap.dedent("""\
    import numpy as np

    float16 = "float16"
    float32 = np.float32
    long = np.int64

    class _NoopCtx:
        def __enter__(self):
            return self
        def __exit__(self, *args):
            return False

    def no_grad():
        return _NoopCtx()

    def autocast(*args, **kwargs):
        return _NoopCtx()

    def _u(x):
        return x.arr if isinstance(x, FakeTensor) else x

    class FakeTensor:
        def __init__(self, arr):
            self.arr = np.asarray(arr)
            self.device = "cpu"
            self.shape = self.arr.shape

        def cuda(self):
            return self

        def cpu(self):
            return self

        def eval(self):
            return self

        def float(self):
            return FakeTensor(self.arr.astype(np.float32))

        def bool(self):
            return FakeTensor(self.arr.astype(bool))

        def numpy(self):
            return self.arr

        def unsqueeze(self, dim):
            return FakeTensor(np.expand_dims(self.arr, axis=dim))

        def sum(self, dim=None):
            return FakeTensor(self.arr.sum(axis=dim))

        def clamp_min(self, val):
            return FakeTensor(np.maximum(self.arr, val))

        def masked_fill(self, mask, val):
            return FakeTensor(np.where(_u(mask), val, self.arr))

        def __getitem__(self, idx):
            return FakeTensor(self.arr[_u(idx)])

        def __mul__(self, other):
            return FakeTensor(self.arr * _u(other))

        def __truediv__(self, other):
            return FakeTensor(self.arr / _u(other))

        def __eq__(self, other):
            return FakeTensor(self.arr == _u(other))

    def tensor(data, dtype=None):
        return FakeTensor(np.array(data, dtype=dtype))

    def zeros(shape, dtype=np.float32):
        return FakeTensor(np.zeros(shape, dtype=dtype))

    def zeros_like(x):
        return FakeTensor(np.zeros_like(_u(x)))

    def where(cond, x, y):
        return FakeTensor(np.where(_u(cond), _u(x), _u(y)))
    """)

STUB_TRANSFORMERS = textwrap.dedent("""\
    import os
    import types
    import torch

    PAD_ID = 0
    CHAR_MAP = {"A": 1, "C": 2, "G": 3, "T": 4, "N": 5}

    class AutoTokenizer:
        @classmethod
        def from_pretrained(cls, model_name, trust_remote_code=False):
            return cls()

        def __call__(self, seqs, return_tensors="pt", padding=True, truncation=True, max_length=2048):
            encoded = [[CHAR_MAP.get(ch, 5) for ch in s[:max_length]] for s in seqs]
            max_len = max((len(x) for x in encoded), default=0)
            ids, masks = [], []
            for row in encoded:
                pad = max_len - len(row)
                ids.append(row + [PAD_ID] * pad)
                masks.append([1] * len(row) + [0] * pad)
            return {
                "input_ids": torch.tensor(ids, dtype=torch.long),
                "attention_mask": torch.tensor(masks, dtype=torch.long),
            }

    class AutoModel:
        def __init__(self):
            self.config = types.SimpleNamespace(hidden_size=4)

        @classmethod
        def from_pretrained(cls, model_name, trust_remote_code=False):
            return cls()

        def cuda(self):
            return self

        def eval(self):
            return self

        def __call__(self, input_ids, attention_mask=None):
            B, L = input_ids.shape
            h = torch.zeros((B, L, self.config.hidden_size), dtype=torch.float32)
            arr = h.numpy()
            ids_np = input_ids.numpy()
            for code in (1, 2, 3, 4):
                arr[:, :, code - 1] = (ids_np == code).astype("float32")
            if os.environ.get("INJECT_NAN_ON_PAD") == "1" and attention_mask is not None:
                h = h.masked_fill(attention_mask.unsqueeze(-1) == 0, float("nan"))
            return (h,)
    """)


def parse_fsb(path: Path):
    raw = path.read_bytes()
    magic = raw[:4]
    d, n = struct.unpack("<iq", raw[4:16])
    mat_end = 16 + 4 * d * n
    assert len(raw) >= mat_end, f"Truncated FSB1 matrix: {len(raw)} < {mat_end}"
    X = np.frombuffer(raw[16:mat_end], dtype=np.float32).reshape(n, d).copy()
    tail = raw[mat_end:].decode("utf-8")
    ids = tail.split("\n") if n > 0 else ([] if tail == "" else tail.split("\n"))
    return magic, d, n, X, ids


class TestEmbedDnabert(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="dnabert-test-")
        self.root = Path(self.tmp.name)
        stub_dir = self.root / "stubs"
        stub_dir.mkdir()
        (stub_dir / "torch.py").write_text(STUB_TORCH)
        (stub_dir / "transformers.py").write_text(STUB_TRANSFORMERS)
        self.stub_dir = stub_dir

    def tearDown(self):
        self.tmp.cleanup()

    def _run(self, fasta: Path, out_fsb: Path, extra_args=(), extra_env=None):
        env = dict(os.environ)
        env["CUDA_VISIBLE_DEVICES"] = ""
        env["PYTHONPATH"] = f"{self.stub_dir}:{env.get('PYTHONPATH', '')}".rstrip(":")
        if extra_env:
            env.update(extra_env)
        return subprocess.run(
            [sys.executable, "-B", str(SCRIPT), str(fasta), str(out_fsb), *extra_args],
            env=env,
            capture_output=True,
            text=True,
        )

    def test_valid_plain_acgt_control_unchanged(self):
        """Unchanged control: valid multi-line plain ACGT FASTA produces identical FSB1 framing and vectors."""
        fasta_txt = (
            ">seq_long description ignored\n"
            "AAAA\n"
            "CCCC\n"
            ">seq_short\n"
            "GG\n"
            ">seq_t\n"
            "TTTT\n"
        )
        plain_fa = self.root / "control.fa"
        plain_fa.write_text(fasta_txt)
        out_plain = self.root / "control.fsb"
        res = self._run(plain_fa, out_plain, extra_args=("--batch", "2"))
        self.assertEqual(res.returncode, 0, res.stderr)

        magic, d, n, X, ids = parse_fsb(out_plain)
        self.assertEqual(magic, b"FSB1")
        self.assertEqual((n, d), (3, 4))
        self.assertEqual(ids, ["seq_long", "seq_short", "seq_t"])
        expected = np.array([
            [0.5, 0.5, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ], dtype=np.float32)
        np.testing.assert_allclose(X, expected, rtol=0, atol=1e-6)

    def test_gzip_input_and_rna_u_normalization(self):
        """Gzipped FASTA input and RNA U bases match equivalent plain T-normalized FASTA output."""
        plain_fa = self.root / "equiv.fa"
        gz_fa = self.root / "equiv.fa.gz"
        plain_fa.write_text(">seq_long\nAAAACCCC\n>seq_rna\nTTTT\n")
        with gzip.open(gz_fa, "wt") as gz:
            gz.write(">seq_long\nAAAACCCC\n>seq_rna\nUUUU\n")

        out_plain = self.root / "equiv_plain.fsb"
        out_gz = self.root / "equiv_gz.fsb"
        self.assertEqual(self._run(plain_fa, out_plain).returncode, 0)
        r_gz = self._run(gz_fa, out_gz)
        self.assertEqual(r_gz.returncode, 0, r_gz.stderr)
        self.assertEqual(out_plain.read_bytes(), out_gz.read_bytes())

    def test_rejects_empty_internal_and_trailing_records_without_writing_corrupt_fsb(self):
        """Empty FASTA records abort cleanly before writing an under-sized FSB1 matrix."""
        for label, content in (
            ("internal_empty", ">seq1\nACGT\n>empty_mid\n>seq2\nGGGG\n"),
            ("whitespace_empty", ">empty_ws\n   \n>seq2\nACGT\n"),
            ("trailing_empty", ">seq1\nACGT\n>empty_end\n"),
        ):
            with self.subTest(case=label):
                fa = self.root / f"{label}.fa"
                out_fsb = self.root / f"{label}.fsb"
                fa.write_text(content)
                res = self._run(fa, out_fsb)
                self.assertNotEqual(res.returncode, 0)
                self.assertIn("empty sequence", res.stderr)
                self.assertFalse(out_fsb.exists())

    def test_padding_nan_does_not_corrupt_shorter_sequences_in_batch(self):
        """Masked padding positions containing NaN do not propagate NaN into pooled embeddings."""
        fa = self.root / "pad_nan.fa"
        fa.write_text(">short\nAA\n>long\nCCCCCCCC\n")
        out_fsb = self.root / "pad_nan.fsb"
        res = self._run(fa, out_fsb, extra_args=("--batch", "2"), extra_env={"INJECT_NAN_ON_PAD": "1"})
        self.assertEqual(res.returncode, 0, res.stderr)
        _, d, n, X, ids = parse_fsb(out_fsb)
        self.assertEqual(ids, ["short", "long"])
        self.assertTrue(np.isfinite(X).all(), f"Non-finite pooled output: {X}")
        np.testing.assert_allclose(X[0], [1.0, 0.0, 0.0, 0.0], rtol=0, atol=1e-6)
        np.testing.assert_allclose(X[1], [0.0, 1.0, 0.0, 0.0], rtol=0, atol=1e-6)


if __name__ == "__main__":
    unittest.main()
