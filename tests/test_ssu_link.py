#!/usr/bin/env python3
"""Tests for ssu_link_set.py and ssu_link_score.py (code-05 F5)."""
import gzip
import os
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PIPELINE = ROOT / "pipeline"
PYTHON = sys.executable

STUB_TORCH = textwrap.dedent("""\
    import numpy as np

    int = np.int32

    def no_grad():
        def deco(fn):
            return fn
        return deco

    def _u(x):
        return x.arr if isinstance(x, FakeTensor) else x

    class FakeTensor:
        def __init__(self, arr):
            self.arr = np.asarray(arr)

        @property
        def device(self):
            return "cpu"

        def unsqueeze(self, dim):
            return FakeTensor(np.expand_dims(self.arr, axis=dim))

        def squeeze(self, dim=None):
            return FakeTensor(np.squeeze(self.arr, axis=dim))

        def to(self, *args, **kwargs):
            return self

        def float(self):
            return FakeTensor(self.arr.astype(np.float32))

        def long(self):
            return FakeTensor(self.arr.astype(np.int64))

        def gather(self, dim, index):
            return FakeTensor(np.take_along_axis(self.arr, _u(index), axis=dim))

        def sum(self, dim=None):
            return FakeTensor(self.arr.sum(axis=dim))

        def __float__(self):
            return float(self.arr.item())

        def __getitem__(self, idx):
            if isinstance(idx, tuple):
                idx = tuple(_u(i) for i in idx)
            else:
                idx = _u(idx)
            return FakeTensor(self.arr[idx])

    def tensor(data, dtype=None, device=None):
        return FakeTensor(np.array(data, dtype=dtype))

    def log_softmax(x, dim=-1):
        arr = _u(x).astype(np.float64)
        m = np.max(arr, axis=dim, keepdims=True)
        shifted = arr - m
        lse = np.log(np.sum(np.exp(shifted), axis=dim, keepdims=True))
        return FakeTensor((shifted - lse).astype(np.float32))
""")

STUB_EVO2 = textwrap.dedent("""\
    import os
    import numpy as np
    import torch

    class _InnerModel:
        def parameters(self):
            yield torch.FakeTensor(np.zeros(1, dtype=np.float32))

        def forward(self, x):
            arr = x.arr
            b, l = arr.shape
            logits = np.zeros((b, max(l, 1), 256), dtype=np.float32)
            mode = os.environ.get("MOCK_EVO2_MODE", "acgt_motif")
            for bi in range(b):
                seq = "".join(chr(int(c)) for c in arr[bi])
                if mode == "acgt_motif" and "ATGCATGC" in seq:
                    for ch in "ACGT":
                        logits[bi, :, ord(ch)] = 5.0
                elif mode == "iupac_rc" and "YYYYKMVBHD" in seq:
                    for ch in "ACGT":
                        logits[bi, :, ord(ch)] = 6.0
            return torch.FakeTensor(logits), None

    class _Tokenizer:
        def tokenize(self, s):
            return [ord(ch) for ch in s]

    class Evo2:
        def __init__(self, model_name="evo2_40b", local_path=None):
            self.tokenizer = _Tokenizer()
            self.model = _InnerModel()
""")


def install_stubs(td: Path) -> Path:
    stubs = td / "stubs"
    (stubs / "evo2").mkdir(parents=True, exist_ok=True)
    (stubs / "torch.py").write_text(STUB_TORCH)
    (stubs / "evo2" / "__init__.py").write_text(STUB_EVO2)
    return stubs


def run_script(name: str, args: list[str], stubs: Path | None = None, extra_env: dict[str, str] | None = None):
    env = os.environ.copy()
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    if stubs is not None:
        env["PYTHONPATH"] = f"{stubs}:{env.get('PYTHONPATH', '')}".rstrip(":")
    if extra_env:
        env.update(extra_env)
    return subprocess.run([PYTHON, "-B", str(PIPELINE / name), *args], text=True, capture_output=True, env=env)


def init_gtdb_fixture(td: Path, ssu_records: str, genome_records: str) -> Path:
    for rel in ("214", "220"):
        (td / f"release{rel}").mkdir(parents=True, exist_ok=True)
        for dom in ("bac120", "ar53"):
            gzip.open(td / f"release{rel}" / f"{dom}_taxonomy_r{rel}.tsv.gz", "wt").close()
    r226 = td / "release226"
    for sub in ("metadata", "ssu", "skani/g"):
        (r226 / sub).mkdir(parents=True, exist_ok=True)
    with gzip.open(r226 / "metadata" / "bac120_metadata_r226.tsv.gz", "wt") as f:
        f.write(
            "accession\tgtdb_representative\tgtdb_taxonomy\n"
            "GB_GCA_000000001.1\tt\td__Bacteria;p__P1;c__C1;o__O1;f__F1;g__G1;s__S1\n"
            "GB_GCA_000000002.1\tt\td__Bacteria;p__P1;c__C1;o__O1;f__F1;g__G2;s__S2\n"
            "GB_GCA_000000003.1\tt\td__Bacteria;p__P1;c__C1;o__O1;f__F2;g__G3;s__S3\n"
        )
    with gzip.open(r226 / "ssu" / "ssu_all_r226.fna.gz", "wt") as f:
        f.write(ssu_records)
    (r226 / "skani" / "genome_paths.tsv").write_text("GCA_000000001.1_genomic.fna.gz g\n")
    with gzip.open(r226 / "skani" / "g" / "GCA_000000001.1_genomic.fna.gz", "wt") as f:
        f.write(genome_records)
    return r226


class TestSsuLink(unittest.TestCase):
    def test_single_operon_clean_control_unchanged(self):
        """Unchanged control: single-operon clean ACGT contig selection and strand scoring stay identical."""
        with tempfile.TemporaryDirectory() as raw:
            td = Path(raw)
            stubs = install_stubs(td)
            ssu_seq = "ACGTACGTACGTACGTACGT"
            clean_ctg = ("A" * 20) + "GCATGCAT" + ("C" * 32)  # rc contains ATGCATGC
            r226 = init_gtdb_fixture(
                td,
                ssu_records=(
                    f">GB_GCA_000000001.1~ctg1\n{ssu_seq}\n"
                    f">GB_GCA_000000002.1~ctg1\n{ssu_seq}\n"
                    f">GB_GCA_000000003.1~ctg1\n{ssu_seq}\n"
                ),
                genome_records=f">ctg1\n{ssu_seq}\n>ctg2\n{clean_ctg}\n",
            )
            out_set = str(td / "ctrl_set")
            r_set = run_script(
                "ssu_link_set.py",
                [str(r226), out_set, "--n", "1", "--contig-len", "50", "--min-contig", "50", "--min-ssu", "10", "--n-family", "1", "--n-order", "1"],
            )
            self.assertEqual(r_set.returncode, 0, r_set.stderr)
            extracted = (td / "ctrl_set.contigs.fa").read_text().splitlines()
            self.assertEqual(extracted[0], ">GCA_000000001.1")
            self.assertEqual(len(extracted[1]), 50)
            self.assertIn("GCATGCAT", extracted[1])

            out_score = str(td / "ctrl_score")
            r_sc = run_script(
                "ssu_link_score.py",
                [out_set, out_score, "--local_path", "dummy.pt"],
                stubs=stubs,
                extra_env={"MOCK_EVO2_MODE": "acgt_motif"},
            )
            self.assertEqual(r_sc.returncode, 0, r_sc.stderr)
            cols = (td / "ctrl_score.tsv").read_text().splitlines()[1].split("\t")
            uncond, cond_fwd, cond_rev, gain = map(float, cols[4:8])
            self.assertAlmostEqual(cond_fwd, uncond, places=3)
            self.assertGreater(cond_rev, uncond + 5.0)
            self.assertAlmostEqual(gain, cond_rev - uncond, places=3)

    def test_ssu_link_set_excludes_secondary_partial_ssu_and_gap_contigs(self):
        """F5: contigs carrying secondary/partial 16S hits or internal N-runs are excluded in favor of clean non-SSU contigs."""
        with tempfile.TemporaryDirectory() as raw:
            td = Path(raw)
            ssu_primary = "ACGTACGTACGTACGTACGT"      # 20 bp on ctg1
            ssu_secondary = "ACGTACGTACGTACGTAC"      # 18 bp on ctg2 (>= min_ssu 10)
            ssu_partial = "ACGTACGT"                  # 8 bp on ctg3 (< min_ssu 10)
            ctg2 = ("T" * 25) + ssu_secondary + ("T" * 27)   # 70 bp (longest, carries secondary 16S)
            ctg3 = ("G" * 28) + ssu_partial + ("G" * 29)     # 65 bp (carries partial 16S)
            ctg4 = ("C" * 25) + ("N" * 8) + ("C" * 27)       # 60 bp (carries assembly gap in central window)
            ctg5 = "A" * 52                                  # 52 bp (clean non-16S contig)
            r226 = init_gtdb_fixture(
                td,
                ssu_records=(
                    f">GB_GCA_000000001.1~ctg1\n{ssu_primary}\n"
                    f">GB_GCA_000000001.1~ctg2\n{ssu_secondary}\n"
                    f">GB_GCA_000000001.1~ctg3\n{ssu_partial}\n"
                    f">GB_GCA_000000002.1~ctg1\n{ssu_primary}\n"
                    f">GB_GCA_000000003.1~ctg1\n{ssu_primary}\n"
                ),
                genome_records=(
                    f">ctg1\n{ssu_primary}\n"
                    f">ctg2\n{ctg2}\n"
                    f">ctg3\n{ctg3}\n"
                    f">ctg4\n{ctg4}\n"
                    f">ctg5\n{ctg5}\n"
                ),
            )
            out_set = str(td / "multi_set")
            r_set = run_script(
                "ssu_link_set.py",
                [str(r226), out_set, "--n", "1", "--contig-len", "50", "--min-contig", "50", "--min-ssu", "10", "--n-family", "1", "--n-order", "1"],
            )
            self.assertEqual(r_set.returncode, 0, r_set.stderr)
            seq = (td / "multi_set.contigs.fa").read_text().splitlines()[1]
            self.assertEqual(seq, "A" * 50)

    def test_ssu_link_score_iupac_reverse_complement(self):
        """F5: ssu_link_score.py complements full IUPAC ambiguity codes on the reverse strand."""
        with tempfile.TemporaryDirectory() as raw:
            td = Path(raw)
            stubs = install_stubs(td)
            # Complement of "YYYYKMVBHD" is "RRRRMKBVDH"; reversing yields "HDVBKMRRRR", so rc(fwd_contig) contains "YYYYKMVBHD".
            fwd_contig = "AAAA" + "HDVBKMRRRR" + "TTTT"
            (td / "iupac.contigs.fa").write_text(f">g1\n{fwd_contig}\n")
            (td / "iupac.ssu.fa").write_text(">g1\nACGTACGTAC\n")
            (td / "iupac.pairs.tsv").write_text("test\tcandidate\trelation\ttaxonomy\ng1\tg1\tself\td__Bacteria\n")

            out_score = str(td / "iupac_score")
            r_sc = run_script(
                "ssu_link_score.py",
                [str(td / "iupac"), out_score, "--local_path", "dummy.pt"],
                stubs=stubs,
                extra_env={"MOCK_EVO2_MODE": "iupac_rc"},
            )
            self.assertEqual(r_sc.returncode, 0, r_sc.stderr)
            cols = (td / "iupac_score.tsv").read_text().splitlines()[1].split("\t")
            uncond, cond_fwd, cond_rev, gain = map(float, cols[4:8])
            self.assertAlmostEqual(cond_fwd, uncond, places=3)
            self.assertGreater(cond_rev, uncond + 5.0)
            self.assertAlmostEqual(gain, cond_rev - uncond, places=3)


if __name__ == "__main__":
    unittest.main(verbosity=2)
