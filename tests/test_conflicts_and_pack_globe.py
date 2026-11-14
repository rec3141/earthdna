#!/usr/bin/env python3
"""Tests for conflicts.py and pack_globe.py (code-05 F9)."""
import base64
import json
import os
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


def decode_flags(globe_dir: Path, n: int):
    raw = np.frombuffer(base64.b64decode((globe_dir / "flags.b64.txt").read_text()), dtype=np.uint8)
    perm = np.random.default_rng(3).permutation(n)
    inv = np.argsort(perm)
    return raw[:n][inv], raw[n:][inv]


class TestConflictsAndPackGlobe(unittest.TestCase):
    def test_named_group_and_domain_conflict_controls_unchanged(self):
        """Unchanged controls: flag=0, flag=1 phylum conflict, and single-phylum flag=2 domain conflict stay identical."""
        with tempfile.TemporaryDirectory() as raw:
            td = Path(raw)
            tax = td / "ctrl.tax.tsv"
            lines = ["id\tdomain\tgroup\tphylum\tclass\torder\tfamily\tgenus\tspecies\tlineage_rank\tlineage\torganism\n"]
            # q0: Bacteria;P0 surrounded by 15 Bacteria;P1 -> flag=1, Bacteria, P1
            lines.append("q0\tBacteria\tP0\tP0\tC\tO\tF\tG\tS\tp\tP0\tOrg\n")
            # q1: Eukaryota;Chloroplastida surrounded by 15 Bacteria;P1 -> flag=2, Bacteria, P1
            lines.append("q1\tEukaryota\tChloroplastida\tChloroplastida\tC\tO\tF\tG\tS\tp\tChloroplastida\tOrg\n")
            for i in range(2, 17):
                lines.append(f"q{i}\tBacteria\tP1\tP1\tC\tO\tF\tG\tS\tp\tP1\tOrg\n")
            tax.write_text("".join(lines))

            N = 17
            graph = np.tile(np.arange(2, 17, dtype=np.int32), (N, 1))
            knn = td / "ctrl.knn"
            graph.tofile(knn)
            ctsv = td / "ctrl.conflicts.tsv"
            rc = run_script("conflicts.py", [str(tax), str(knn), "15", str(ctsv)])
            self.assertEqual(rc.returncode, 0, rc.stderr)

            rows = {l.split("\t")[0]: l.rstrip("\n").split("\t") for l in ctsv.read_text().splitlines()[1:]}
            self.assertEqual(rows["q0"], ["q0", "1", "Bacteria", "P1"])
            self.assertEqual(rows["q1"], ["q1", "2", "Bacteria", "P1"])
            self.assertEqual(rows["q2"], ["q2", "0", "Bacteria", "P1"])

            layout = td / "map.tsv"
            layout.write_text("id\tlon\tlat\n" + "".join(f"q{i}\t0.0\t0.0\n" for i in range(N)))
            gdir = td / "globe"
            rp = run_script("pack_globe.py", [str(tax), str(gdir), f"--layout=similarity={layout}", f"--flags={ctsv}"])
            self.assertEqual(rp.returncode, 0, rp.stderr)

            meta = json.loads((gdir / "globe.json").read_text())
            flags_arr, maj_arr = decode_flags(gdir, N)
            p1_idx = meta["groups"].index("Bacteria · P1")
            self.assertEqual(int(flags_arr[0]), 1)
            self.assertEqual(int(maj_arr[0]), p1_idx)
            self.assertEqual(int(flags_arr[1]), 2)
            self.assertEqual(int(maj_arr[1]), p1_idx)
            self.assertEqual(int(flags_arr[2]), 0)
            self.assertEqual(int(maj_arr[2]), 255)

    def test_domain_conflict_preserves_majority_domain_and_ungrouped_globe_index(self):
        """F9: multi-phylum domain conflict keeps majority domain (maj=255); ungrouped domain conflict maps to 'Domain · ?'."""
        with tempfile.TemporaryDirectory() as raw:
            td = Path(raw)
            tax = td / "c.tax.tsv"
            lines = ["id\tdomain\tgroup\tphylum\tclass\torder\tfamily\tgenus\tspecies\tlineage_rank\tlineage\torganism\n"]
            for i in range(4):
                lines.append(f"q{i}\tEukaryota\tChloroplastida\tChloroplastida\tC\tO\tF\tG\tS\tp\tChloroplastida\tOrg\n")
            for i in range(4, 16):
                ph = f"P{(i - 4) // 2}"
                lines.append(f"q{i}\tBacteria\t{ph}\t{ph}\tC\tO\tF\tG\tS\tp\t{ph}\tOrg\n")
            lines.append("q16\tEukaryota\tChloroplastida\tChloroplastida\tC\tO\tF\tG\tS\tp\tChloroplastida\tOrg\n")
            for i in range(17, 32):
                lines.append(f"q{i}\tBacteria\t\t\tC\tO\tF\tG\tS\tp\t\tOrg\n")
            tax.write_text("".join(lines))

            N = 32
            graph = np.zeros((N, 15), dtype=np.int32)
            graph[:] = np.arange(1, 16, dtype=np.int32)
            graph[16:] = np.arange(17, 32, dtype=np.int32)
            knn = td / "c.knn"
            graph.tofile(knn)
            ctsv = td / "c.conflicts.tsv"
            rc = run_script("conflicts.py", [str(tax), str(knn), "15", str(ctsv)])
            self.assertEqual(rc.returncode, 0, rc.stderr)

            rows = {l.split("\t")[0]: l.rstrip("\n").split("\t") for l in ctsv.read_text().splitlines()[1:]}
            # q0: 12/15 Bacteria across 6 phyla + 3/15 Eukaryota;Chloroplastida -> Bacteria with no single majority phylum.
            self.assertEqual(rows["q0"], ["q0", "2", "Bacteria", ""])
            # q16: 15/15 ungrouped Bacteria -> Bacteria with '?' ungrouped majority.
            self.assertEqual(rows["q16"], ["q16", "2", "Bacteria", "?"])
            # Ungrouped Bacteria surrounded by ungrouped Bacteria remain flag=0.
            self.assertEqual(rows["q17"][1], "0")

            layout = td / "map.tsv"
            layout.write_text("id\tlon\tlat\n" + "".join(f"q{i}\t0.0\t0.0\n" for i in range(N)))
            gdir = td / "globe"
            rp = run_script("pack_globe.py", [str(tax), str(gdir), f"--layout=similarity={layout}", f"--flags={ctsv}"], extra_env={"PYTHONWARNINGS": "error::ResourceWarning"})
            self.assertEqual(rp.returncode, 0, rp.stderr)

            meta = json.loads((gdir / "globe.json").read_text())
            unk_bac_idx = meta["groups"].index("Bacteria · ?")
            flags_arr, maj_arr = decode_flags(gdir, N)
            self.assertEqual(int(flags_arr[0]), 2)
            self.assertEqual(int(maj_arr[0]), 255)
            self.assertEqual(int(flags_arr[16]), 2)
            self.assertEqual(int(maj_arr[16]), unk_bac_idx)


if __name__ == "__main__":
    unittest.main(verbosity=2)
