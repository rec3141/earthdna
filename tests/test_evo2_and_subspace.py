#!/usr/bin/env python3
"""Regression tests for pipeline/models_prep.py (code-05 F3) and pipeline/embed_evo2.py (code-05 F4).

Runs production scripts via subprocess with pure-NumPy stubs for torch/evo2 so no GPU or
external model weights are required. Includes unchanged-case controls alongside desired
error/isolation checks.
"""
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
MODELS_PREP = REPO / "pipeline" / "models_prep.py"
EMBED_EVO2 = REPO / "pipeline" / "embed_evo2.py"

STUB_TORCH = textwrap.dedent("""\
    import os
    import numpy as np

    builtins_int = __builtins__["int"] if isinstance(__builtins__, dict) else __builtins__.int
    int = np.int32

    class OutOfMemoryError(RuntimeError):
        pass

    class _Ctx:
        def __enter__(self): return self
        def __exit__(self, *a): return False
        def __call__(self, fn):
            def w(*a, **kw):
                with self: return fn(*a, **kw)
            return w

    def inference_mode(): return _Ctx()

    class cuda:
        @staticmethod
        def device_count(): return 1
        @staticmethod
        def empty_cache(): pass

    def _u(x):
        return x.arr if isinstance(x, FakeTensor) else x

    class FakeTensor:
        def __init__(self, arr):
            self.arr = np.asarray(arr)
            self.device = "cpu"
            self.shape = self.arr.shape

        def to(self, _d): return self
        def cpu(self): return self
        def detach(self): return FakeTensor(self.arr.copy())
        def clone(self): return FakeTensor(self.arr.copy())
        def float(self): return FakeTensor(self.arr.astype(np.float32))
        def numpy(self): return self.arr
        def sum(self, dim=None): return FakeTensor(self.arr.sum(axis=dim))
        def norm(self, dim=-1, keepdim=False): return FakeTensor(np.linalg.norm(self.arr, axis=dim, keepdims=keepdim))
        def clamp_min(self, val): return FakeTensor(np.maximum(self.arr, val))
        def __getitem__(self, idx):
            if isinstance(idx, tuple): idx = tuple(_u(i) for i in idx)
            else: idx = _u(idx)
            return FakeTensor(self.arr[idx])
        def __setitem__(self, idx, val):
            if isinstance(idx, tuple): idx = tuple(_u(i) for i in idx)
            else: idx = _u(idx)
            self.arr[idx] = _u(val)
        def __mul__(self, o): return FakeTensor(self.arr * _u(o))
        def __truediv__(self, o): return FakeTensor(self.arr / _u(o))
        def __floordiv__(self, o): return FakeTensor(self.arr // _u(o))
        def __lt__(self, o): return FakeTensor(self.arr < _u(o))
        def __ge__(self, o): return FakeTensor(self.arr >= _u(o))
        def __and__(self, o): return FakeTensor(self.arr & _u(o))

    def tensor(data, dtype=None, device=None):
        dt = np.int32 if dtype is int else dtype
        return FakeTensor(np.array(data, dtype=dt))

    def full(shape, fill_value, dtype=None, device=None):
        dt = np.int32 if dtype is int else dtype
        return FakeTensor(np.full(shape, fill_value, dtype=dt))

    def arange(end, device=None): return FakeTensor(np.arange(end))
    def zeros_like(x): return FakeTensor(np.zeros_like(_u(x)))
    def where(c, x, y): return FakeTensor(np.where(_u(c), _u(x), _u(y)))
    """)

STUB_EVO2 = textwrap.dedent("""\
    import os
    import numpy as np
    import torch

    class _Handle:
        def __init__(self, hooks, fn): self._h, self._fn = hooks, fn
        def remove(self): self._h.remove(self._fn)

    class _Layer:
        def __init__(self): self.hooks = []
        def register_forward_hook(self, fn):
            self.hooks.append(fn)
            return _Handle(self.hooks, fn)

    class _Tokenizer:
        def tokenize(self, s): return [ord(c) for c in s]

    class _InnerModel:
        def __init__(self): self.layers = {}
        def get_submodule(self, name):
            return self.layers.setdefault(name, _Layer())
        def forward(self, x):
            arr = x.numpy()
            B, L = arr.shape
            mode = os.environ.get("EVO2_STUB_MODE", "normal")
            h = np.zeros((B, L, 4), dtype=np.float32)
            for idx, ch in enumerate("ACGT"):
                h[:, :, idx] = (arr == ord(ch)).astype(np.float32)
            if mode == "zeros":
                h.fill(0.0)
            elif mode == "pad_nan":
                h[arr == ord("N"), :] = np.nan
            t = torch.FakeTensor(h)
            for lyr in self.layers.values():
                for fn in list(lyr.hooks):
                    fn(lyr, (x,), t)
            return t, None

    class Evo2:
        def __init__(self, model_name, local_path=None, use_kernels=False):
            self.tokenizer = _Tokenizer()
            self.model = _InnerModel()
    """)


def write_fsb(path: Path, X: np.ndarray, ids: list[str]):
    X32 = np.asarray(X, dtype=np.float32)
    with path.open("wb") as fh:
        fh.write(b"FSB1" + struct.pack("<iq", X32.shape[1], len(ids)))
        fh.write(X32.tobytes())
        fh.write("\n".join(ids).encode("utf-8"))


def read_fsb(path: Path):
    raw = path.read_bytes()
    magic = raw[:4]
    d, n = struct.unpack("<iq", raw[4:16])
    mat_end = 16 + 4 * d * n
    assert len(raw) >= mat_end, f"Truncated FSB1: {len(raw)} < {mat_end}"
    X = np.frombuffer(raw[16:mat_end], dtype=np.float32).reshape(n, d).copy()
    tail = raw[mat_end:].decode("utf-8")
    ids = tail.split("\n") if n > 0 else ([] if tail == "" else tail.split("\n"))
    return magic, d, n, X, ids


class TestModelsPrepAndEmbedEvo2(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="earthdna-f3f4-")
        self.root = Path(self.tmp.name)
        self.stubs = self.root / "stubs"
        self.stubs.mkdir()
        (self.stubs / "torch.py").write_text(STUB_TORCH)
        (self.stubs / "evo2.py").write_text(STUB_EVO2)

    def tearDown(self):
        self.tmp.cleanup()

    def _env(self, extra=None):
        env = dict(os.environ)
        env["CUDA_VISIBLE_DEVICES"] = ""
        env["PYTHONPATH"] = f"{self.stubs}:{env.get('PYTHONPATH', '')}".rstrip(":")
        if extra:
            env.update(extra)
        return env

    def test_models_prep_valid_subspace_control_unchanged(self):
        """Unchanged control: valid positive-sd, non-zero-norm subspace produces unit-norm finite rows."""
        tax = self.root / "set.tax.tsv"
        tax.write_text("id\tdomain\tgroup\tphylum\tclass\torder\nseq1\tBacteria\tP1\tP1\tC1\tO1\nseq2\tArchaea\tP2\tP2\tC2\tO2\n")
        in_fsb = self.root / "in.fsb"
        write_fsb(in_fsb, np.array([[2.0, 0.0], [0.0, 4.0]], dtype=np.float32), ["seq1", "seq2"])
        mdl = self.root / "sub.npz"
        np.savez(mdl, mu=np.array([0.0, 0.0]), sd=np.array([1.0, 2.0]), b=np.array([0.0, 0.0]), W=np.eye(2))
        out_pref = self.root / "out_ok"
        res = subprocess.run(
            [sys.executable, "-B", str(MODELS_PREP), str(tax), str(out_pref), "--subspace", f"sub={in_fsb}={mdl}"],
            env=self._env(),
            capture_output=True,
            text=True,
        )
        self.assertEqual(res.returncode, 0, res.stderr)
        magic, d, n, X, ids = read_fsb(Path(f"{out_pref}.sub.fsb"))
        self.assertEqual(magic, b"FSB1")
        self.assertEqual(ids, ["seq1", "seq2"])
        np.testing.assert_allclose(X, [[1.0, 0.0], [0.0, 1.0]], rtol=0, atol=1e-6)

    def test_models_prep_rejects_zero_sd_and_zero_row_norm(self):
        """F3: models_prep.py fails fast instead of writing silent NaN rows/columns into .fsb."""
        tax = self.root / "set.tax.tsv"
        tax.write_text("id\tdomain\tgroup\tphylum\tclass\torder\nseq1\tBacteria\tP1\tP1\tC1\tO1\nseq2\tArchaea\tP2\tP2\tC2\tO2\n")
        in_fsb = self.root / "in.fsb"
        write_fsb(in_fsb, np.array([[0.0, 0.0], [1.0, 2.0]], dtype=np.float32), ["seq1", "seq2"])

        zero_norm_mdl = self.root / "zero_norm.npz"
        np.savez(zero_norm_mdl, mu=np.zeros(2), sd=np.ones(2), b=np.zeros(2), W=np.eye(2))
        out_zn = self.root / "out_zn"
        r_zn = subprocess.run(
            [sys.executable, "-B", str(MODELS_PREP), str(tax), str(out_zn), "--subspace", f"s={in_fsb}={zero_norm_mdl}"],
            env=self._env(),
            capture_output=True,
            text=True,
        )
        self.assertNotEqual(r_zn.returncode, 0)
        self.assertFalse(Path(f"{out_zn}.s.fsb").exists())

        zero_sd_mdl = self.root / "zero_sd.npz"
        np.savez(zero_sd_mdl, mu=np.zeros(2), sd=np.array([1.0, 0.0]), b=np.zeros(2), W=np.eye(2))
        out_zs = self.root / "out_zs"
        r_zs = subprocess.run(
            [sys.executable, "-B", str(MODELS_PREP), str(tax), str(out_zs), "--subspace", f"s={in_fsb}={zero_sd_mdl}"],
            env=self._env(),
            capture_output=True,
            text=True,
        )
        self.assertNotEqual(r_zs.returncode, 0)
        self.assertFalse(Path(f"{out_zs}.s.fsb").exists())

    def test_embed_evo2_valid_control_unchanged(self):
        """Unchanged control: finite non-empty sequences embed into chunks and assemble into valid FSB1."""
        fa = self.root / "ok.fa"
        fa.write_text(">s1\nAAAA\n>s2\nCCCCGGGG\n")
        out_pref = self.root / "evo_ok"
        r1 = subprocess.run(
            [sys.executable, "-B", str(EMBED_EVO2), str(fa), str(out_pref), "--layers", "L1", "--pools", "mean", "--chunk", "2", "--check", "2"],
            env=self._env(),
            capture_output=True,
            text=True,
        )
        self.assertEqual(r1.returncode, 0, r1.stderr)
        self.assertNotIn("WARN", r1.stdout)
        r2 = subprocess.run(
            [sys.executable, "-B", str(EMBED_EVO2), str(fa), str(out_pref), "--layers", "L1", "--pools", "mean", "--chunk", "2", "--assemble"],
            env=self._env(),
            capture_output=True,
            text=True,
        )
        self.assertEqual(r2.returncode, 0, r2.stderr)
        magic, d, n, X, ids = read_fsb(Path(f"{out_pref}.L1.mean.fsb"))
        self.assertEqual((magic, d, n, ids), (b"FSB1", 4, 2, ["s1", "s2"]))
        np.testing.assert_allclose(X, [[1.0, 0.0, 0.0, 0.0], [0.0, 0.5, 0.5, 0.0]], rtol=0, atol=1e-6)

    def test_embed_evo2_check_warns_on_zero_vectors_and_isolates_pad_nan(self):
        """F4: --check warns on zero/NaN cosines and right-padded NaNs do not poison shorter batched sequences."""
        fa = self.root / "two.fa"
        fa.write_text(">short\nAA\n>long\nCCCCGGGG\n")

        out_z = self.root / "evo_zero"
        rz = subprocess.run(
            [sys.executable, "-B", str(EMBED_EVO2), str(fa), str(out_z), "--layers", "L1", "--pools", "mean", "--check", "2", "--end", "0"],
            env=self._env({"EVO2_STUB_MODE": "zeros"}),
            capture_output=True,
            text=True,
        )
        self.assertEqual(rz.returncode, 0, rz.stderr)
        self.assertIn("WARN batched vs unbatched cosine", rz.stdout)

        out_p = self.root / "evo_padnan"
        rp = subprocess.run(
            [sys.executable, "-B", str(EMBED_EVO2), str(fa), str(out_p), "--layers", "L1", "--pools", "mean", "--check", "2"],
            env=self._env({"EVO2_STUB_MODE": "pad_nan"}),
            capture_output=True,
            text=True,
        )
        self.assertEqual(rp.returncode, 0, rp.stderr)
        self.assertNotIn("WARN", rp.stdout)
        ra = subprocess.run(
            [sys.executable, "-B", str(EMBED_EVO2), str(fa), str(out_p), "--layers", "L1", "--pools", "mean", "--assemble"],
            env=self._env(),
            capture_output=True,
            text=True,
        )
        self.assertEqual(ra.returncode, 0, ra.stderr)
        _, _, _, X, _ = read_fsb(Path(f"{out_p}.L1.mean.fsb"))
        self.assertTrue(np.isfinite(X).all())
        np.testing.assert_allclose(X[0], [1.0, 0.0, 0.0, 0.0], rtol=0, atol=1e-6)

    def test_embed_evo2_empty_assemble_header_and_stale_chunk_rejection(self):
        """F4: empty FASTA --assemble writes valid 0-row FSB1 header; stale mismatched chunk row count aborts."""
        empty_fa = self.root / "empty.fa"
        empty_fa.write_text("")
        out_empty = self.root / "evo_empty"
        re_res = subprocess.run(
            [sys.executable, "-B", str(EMBED_EVO2), str(empty_fa), str(out_empty), "--layers", "L1", "--pools", "mean", "--assemble"],
            env=self._env(),
            capture_output=True,
            text=True,
        )
        self.assertEqual(re_res.returncode, 0, re_res.stderr)
        magic, d, n, X, ids = read_fsb(Path(f"{out_empty}.L1.mean.fsb"))
        self.assertEqual((magic, d, n, ids), (b"FSB1", 0, 0, []))
        self.assertEqual(X.shape, (0, 0))

        fa3 = self.root / "three.fa"
        fa3.write_text(">s1\nAA\n>s2\nCC\n>s3\nGG\n")
        out_stale = self.root / "evo_stale"
        parts = Path(f"{out_stale}.parts")
        parts.mkdir()
        np.savez(parts / "000000.npz", **{"L1.mean": np.zeros((1, 4), dtype=np.float32)})
        rs = subprocess.run(
            [sys.executable, "-B", str(EMBED_EVO2), str(fa3), str(out_stale), "--layers", "L1", "--pools", "mean", "--chunk", "4", "--assemble"],
            env=self._env(),
            capture_output=True,
            text=True,
        )
        self.assertNotEqual(rs.returncode, 0)
        self.assertFalse(Path(f"{out_stale}.L1.mean.fsb").exists())


if __name__ == "__main__":
    unittest.main()
