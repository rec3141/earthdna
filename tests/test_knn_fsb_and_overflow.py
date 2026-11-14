#!/usr/bin/env python3
"""Tests for knn_eval.py, gtdb_reselect.py, gtdb_prep.py, gtdb_read_refs.py, and kmers.py (code-05 F2, F7, F8)."""
import gzip
import os
import struct
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PIPELINE = ROOT / "pipeline"
PYTHON = sys.executable


def run_script(name: str, args: list[str], extra_env: dict[str, str] | None = None):
    env = os.environ.copy()
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    if extra_env:
        env.update(extra_env)
    return subprocess.run([PYTHON, "-B", str(PIPELINE / name), *args], text=True, capture_output=True, env=env)


def write_fsb(path: Path, X: np.ndarray, ids: list[str], magic: bytes = b"FSB1") -> None:
    arr = np.asarray(X, dtype=np.float32)
    with open(path, "wb") as f:
        f.write(magic + struct.pack("<iq", arr.shape[1], len(ids)))
        f.write(arr.tobytes())
        f.write("\n".join(ids).encode())


def read_fsb_raw(path: Path):
    raw = path.read_bytes()
    magic = raw[:4]
    d, n = struct.unpack("<iq", raw[4:16])
    body = raw[16:]
    return magic, d, n, body


class TestKnnFsbAndOverflow(unittest.TestCase):
    def test_valid_controls_unchanged(self):
        """Unchanged controls: finite knn_eval (:sqrt/:z), populated gtdb_reselect, k=5 prep, and k=4 palindrome conventions."""
        with tempfile.TemporaryDirectory() as raw:
            td = Path(raw)
            tax = td / "ctrl.tax.tsv"
            tax.write_text(
                "id\tdomain\tgroup\tphylum\tclass\torder\tfamily\tgenus\n"
                "s0\tB\tP1\tP1\tC1\tO1\tF1\tG1\n"
                "s1\tB\tP1\tP1\tC1\tO1\tF1\tG1\n"
                "s2\tB\tP2\tP2\tC2\tO2\tF2\tG2\n"
                "s3\tB\tP2\tP2\tC2\tO2\tF2\tG2\n"
            )
            X = np.array([[1.0, 0.0], [1.1, 0.0], [0.0, 4.0], [0.0, 4.2]], dtype=np.float32)
            fsb = td / "ctrl.fsb"
            write_fsb(fsb, X, ["s0", "s1", "s2", "s3"])
            rk = run_script("knn_eval.py", [str(tax), str(fsb), f"{fsb}:sqrt", f"{fsb}:z", "--k", "1", "--pca", "0"])
            self.assertEqual(rk.returncode, 0, rk.stderr)
            lines = rk.stdout.strip().splitlines()[1:]
            self.assertEqual(len(lines), 3)
            for ln in lines:
                self.assertEqual(ln.split("\t")[1:], ["1.000"] * 6)

            # Palindrome counting conventions remain distinct (--canonical -> 1.0, --rc -> 2.0 on even k=4 AATT).
            fa_pal = td / "pal.fa"
            fa_pal.write_text(">pal\nAATT\n")
            c4, r4 = td / "c4.fsb", td / "r4.fsb"
            self.assertEqual(run_script("kmers.py", [str(fa_pal), str(c4), "--k", "4", "--canonical"]).returncode, 0)
            self.assertEqual(run_script("kmers.py", [str(fa_pal), str(r4), "--k", "4", "--rc"]).returncode, 0)
            _, dc, nc, bc = read_fsb_raw(c4)
            _, dr, nr, br = read_fsb_raw(r4)
            Xc = np.frombuffer(bc[: 4 * dc * nc], np.float32)
            Xr = np.frombuffer(br[: 4 * dr * nr], np.float32)
            self.assertEqual(float(Xc.max()), 1.0)
            self.assertEqual(float(Xr.max()), 2.0)

            # Populated shard in gtdb_reselect.py selects min-per-genome rows cleanly.
            shard = td / "pop.fsb"
            write_fsb(shard, np.ones((2, 4), np.float32), ["GCA_000000001.1_ctg1_0", "GCA_000000001.1_ctg1_100"])
            shards_txt = td / "shards.txt"
            shards_txt.write_text(f"{shard}\n")
            genomes = td / "genomes.tsv"
            genomes.write_text("accession\tdomain\tgroup\tphylum\tclass\torder\nGCA_000000001.1\tB\tP\tP\tC\tO\n")
            rr = run_script("gtdb_reselect.py", [str(shards_txt), str(genomes), str(td / "re"), "--ref", "0.0", "--min-per-genome", "2"])
            self.assertEqual(rr.returncode, 0, rr.stderr)
            mg, drs, nrs, _ = read_fsb_raw(td / "re.ref.fsb")
            self.assertEqual((mg, drs, nrs), (b"FSB1", 4, 2))

    def test_knn_eval_rejects_nan_neg_sqrt_bad_magic_and_unknown_transform(self):
        """F2, F7: knn_eval.py rejects NaN matrices, negative :sqrt values, bad FSB1 headers, and unknown :transforms."""
        with tempfile.TemporaryDirectory() as raw:
            td = Path(raw)
            tax = td / "u.tax.tsv"
            tax.write_text(
                "id\tdomain\tgroup\tphylum\tclass\torder\tfamily\tgenus\n"
                + "".join(f"s{i}\tB\tP\tP\tC\tO\tF\tG{i}\n" for i in range(10))
            )
            ids = [f"s{i}" for i in range(10)]

            X_nan = np.ones((10, 4), np.float32)
            X_nan[0, 0] = np.nan
            nan_fsb = td / "nan.fsb"
            write_fsb(nan_fsb, X_nan, ids)
            self.assertNotEqual(run_script("knn_eval.py", [str(tax), str(nan_fsb), "--k", "2", "--pca", "0"]).returncode, 0)

            X_neg = np.ones((10, 4), np.float32)
            X_neg[0, 0] = -1.0
            neg_fsb = td / "neg.fsb"
            write_fsb(neg_fsb, X_neg, ids)
            self.assertNotEqual(run_script("knn_eval.py", [str(tax), f"{neg_fsb}:sqrt", "--k", "2", "--pca", "0"]).returncode, 0)

            bad_fsb = td / "bad.fsb"
            write_fsb(bad_fsb, np.ones((10, 4), np.float32), ids, magic=b"NOPE")
            self.assertNotEqual(run_script("knn_eval.py", [str(tax), str(bad_fsb), "--k", "2", "--pca", "0"]).returncode, 0)

            ok_fsb = td / "ok.fsb"
            write_fsb(ok_fsb, np.ones((10, 4), np.float32), ids)
            self.assertNotEqual(run_script("knn_eval.py", [str(tax), f"{ok_fsb}:zscore", "--k", "2", "--pca", "0"]).returncode, 0)

    def test_knn_eval_float64_distance_avoids_norm_cancellation(self):
        """F2: float64 squared distances preserve true 1-NN ordering on +/-1e4 clusters where float32 cancels."""
        with tempfile.TemporaryDirectory() as raw:
            td = Path(raw)
            tax = td / "cancel.tax.tsv"
            genera = ["GA", "GB", "GA", "GC", "GD", "GE"]
            tax.write_text(
                "id\tdomain\tgroup\tphylum\tclass\torder\tfamily\tgenus\n"
                + "".join(f"s{i}\tB\tP\tP\tC\tO\tF\t{g}\n" for i, g in enumerate(genera))
            )
            X = np.array(
                [
                    [10000.0, 0.0],
                    [10000.0, 0.20],
                    [10000.0, 0.05],
                    [-10000.0, 0.0],
                    [-10000.0, -0.20],
                    [-10000.0, -0.05],
                ],
                dtype=np.float32,
            )
            fsb = td / "cancel.fsb"
            write_fsb(fsb, X, [f"s{i}" for i in range(6)])
            res = run_script("knn_eval.py", [str(tax), str(fsb), "--k", "1", "--pca", "0"])
            self.assertEqual(res.returncode, 0, res.stderr)
            genus_acc = float(res.stdout.strip().splitlines()[-1].split("\t")[-1])
            self.assertAlmostEqual(genus_acc, 2.0 / 6.0, places=3)

    def test_gtdb_reselect_handles_empty_shard_without_keyerror_or_fd_leak(self):
        """F7: gtdb_reselect.py processes 0-row FSB1 shards cleanly with zero ResourceWarnings."""
        with tempfile.TemporaryDirectory() as raw:
            td = Path(raw)
            short_fa = td / "short.fa"
            short_fa.write_text(">s0\nACGT\n")
            empty_fsb = td / "empty.fsb"
            self.assertEqual(run_script("kmers.py", [str(short_fa), str(empty_fsb), "--k", "2", "--min-len", "99"]).returncode, 0)

            pop_fsb = td / "pop.fsb"
            write_fsb(pop_fsb, np.ones((1, 16), np.float32), ["GCA_000000001.1_ctg1_0"])
            shards = td / "shards.txt"
            shards.write_text(f"{empty_fsb}\n{pop_fsb}\n")
            genomes = td / "genomes.tsv"
            genomes.write_text("accession\tdomain\tgroup\tphylum\tclass\torder\nGCA_000000001.1\tB\tP\tP\tC\tO\n")

            rr = run_script(
                "gtdb_reselect.py",
                [str(shards), str(genomes), str(td / "out"), "--ref", "0.0", "--min-per-genome", "1"],
                extra_env={"PYTHONWARNINGS": "error::ResourceWarning"},
            )
            self.assertEqual(rr.returncode, 0, rr.stderr)
            mg, d, n, body = read_fsb_raw(td / "out.ref.fsb")
            self.assertEqual((mg, d, n), (b"FSB1", 16, 1))
            self.assertEqual(body[4 * d * n :].decode().splitlines(), ["GCA_000000001.1_ctg1_0"])

    def test_gtdb_prep_k9_no_uint16_overflow_and_cli_guards(self):
        """F8: gtdb_prep.py --k 9 matches kmers.py --k 9 --canonical without uint16 wrap; CLI guards reject invalid args."""
        with tempfile.TemporaryDirectory() as raw:
            td = Path(raw)
            gtdb = td / "gtdb9"
            (gtdb / "taxonomy").mkdir(parents=True)
            (gtdb / "skani" / "g").mkdir(parents=True)
            (gtdb / "taxonomy" / "gtdb_taxonomy.tsv").write_text("GB_GCA_000000001.1\td__Bacteria;p__P;c__C;o__O;f__F;g__G;s__S\n")
            (gtdb / "skani" / "genome_paths.tsv").write_text("GCA_000000001.1_genomic.fna.gz g\n")
            seq = "TTTTTTTTTTTTACGT"  # 16 bp window containing TTTTTTTTT (code 262143, canonical AAAAAAAAA = col 0)
            with gzip.open(gtdb / "skani" / "g" / "GCA_000000001.1_genomic.fna.gz", "wt") as f:
                f.write(f">ctg1\n{seq}\n")

            prep_out = td / "prep9"
            rp = run_script("gtdb_prep.py", [str(gtdb), str(prep_out), "--k", "9", "--frag", "16", "--ref", "1.0", "--workers", "1"])
            self.assertEqual(rp.returncode, 0, rp.stderr)

            ref_fa = td / "single.fa"
            ref_fa.write_text(f">w0\n{seq}\n")
            km_fsb = td / "km9.fsb"
            rk = run_script("kmers.py", [str(ref_fa), str(km_fsb), "--k", "9", "--canonical"])
            self.assertEqual(rk.returncode, 0, rk.stderr)

            _, dp, np_, bp = read_fsb_raw(Path(f"{prep_out}.all.000.fsb"))
            _, dk, nk, bk = read_fsb_raw(km_fsb)
            Xp = np.frombuffer(bp[: 4 * dp * np_], np.float32).reshape(np_, dp)
            Xk = np.frombuffer(bk[: 4 * dk * nk], np.float32).reshape(nk, dk)
            self.assertEqual(float(Xp[0, 0]), 4.0)
            np.testing.assert_array_equal(Xp, Xk)

            # Mutually exclusive --canonical and --rc in kmers.py:
            rboth = run_script("kmers.py", [str(ref_fa), str(td / "both.fsb"), "--k", "5", "--canonical", "--rc"])
            self.assertNotEqual(rboth.returncode, 0)

            # Clean SystemExit (no ValueError traceback) when --len < --k in gtdb_read_refs.py:
            ref_gz = td / "ref_short.fa.gz"
            with gzip.open(ref_gz, "wt") as f:
                f.write(">p1\nACGTACGT\n")
            rshort = run_script("gtdb_read_refs.py", [str(ref_gz), str(td / "rrefs"), "--len", "3", "--step", "3", "--k", "5", "--workers", "1"])
            self.assertNotEqual(rshort.returncode, 0)
            self.assertNotIn("ValueError", rshort.stderr)


if __name__ == "__main__":
    unittest.main(verbosity=2)
