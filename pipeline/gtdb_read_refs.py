#!/usr/bin/env python3
"""GTDB reference fragments -> read-length pieces that carry their parent's map position, for placing reads.

The counterpart of region_refs.py for genome globes: a 3-5 kb nanopore read's k-mer profile is not comparable with
a 10 kb fragment's, so reads are projected with a model built from the reference fragments cut to read length.
Each reference fragment (gtdb_prep.py's OUT.ref.fasta.gz) gives pieces of --len bp at offsets 0, --step, 2*step ...
that fit inside it; a piece keeps the parent's position in every layout (fastsne --fixed-y), so a read lands where
the fragment it resembles sits on the globe.

  OUT.fsb               canonical k-mer counts per piece (gtdb_prep.py's columns), ids <parent>@<offset>
  OUT.<layout>.tsv      fixed-y map per --layout NAME=parent_layout.tsv: piece id, parent lon, parent lat

then, per layout:
  fastsne -i OUT.fsb -o /dev/null --transform sqrt --pca 128 --sphere 35 --fixed-y OUT.sculpt.tsv --save-model M
  kmers.py reads.fa reads.fsb --k 5 --canonical;  fastsne --project M -i reads.fsb -o reads.placed.tsv

usage: gtdb_read_refs.py REF.fasta.gz OUT [--len 4000] [--step 3000] [--k 5] [--workers 16] [--maps-only]
                         --layout sculpt=gtdb226.sculpt.tsv [--layout similarity=gtdb226.sphere.tsv]
"""
import argparse, gzip, os, struct
from multiprocessing import Pool
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("ref"); ap.add_argument("out")
ap.add_argument("--len", type=int, default=4000); ap.add_argument("--step", type=int, default=3000)
ap.add_argument("--k", type=int, default=5); ap.add_argument("--workers", type=int, default=os.cpu_count())
ap.add_argument("--layout", action="append", default=[])
ap.add_argument("--batch", type=int, default=20000, help="fragments per worker task")
ap.add_argument("--maps-only", action="store_true", help="reuse OUT.fsb (its ids) and only write the layout maps")
a = ap.parse_args()
k = a.k

# canonical columns exactly as gtdb_prep.py builds them
lut = np.full(256, 255, np.uint8)
for i, c in enumerate(b"ACGT"): lut[c] = i; lut[c + 32] = i
pw = (4 ** np.arange(k - 1, -1, -1)).astype(np.int64)
idx = np.arange(4 ** k)
digits = (idx[:, None] // pw[None, :]) % 4
rc = ((3 - digits[:, ::-1]) * pw[None, :]).sum(1)
canon = np.unique(np.minimum(idx, rc))
col = np.zeros(4 ** k, np.int64); col[canon] = np.arange(len(canon)); col[rc[canon]] = np.arange(len(canon))
D = len(canon)

def counts(s):
    c = lut[np.frombuffer(s, np.uint8)]
    bad = c > 3
    c = np.where(bad, 0, c).astype(np.int64)
    w = np.lib.stride_tricks.sliding_window_view(c, k) @ pw
    ok = ~np.lib.stride_tricks.sliding_window_view(bad, k).any(1)
    return np.bincount(col[w[ok]], minlength=D).astype(np.float32)

def work(batch):
    ids, rows = [], []
    for pid, seq in batch:
        for off in range(0, len(seq) - a.len + 1, a.step):
            ids.append(f"{pid}@{off}"); rows.append(counts(seq[off:off + a.len]))
    return ids, (np.stack(rows) if rows else np.zeros((0, D), np.float32))

def batches():
    b, pid, seq = [], None, []
    with gzip.open(a.ref, "rb") as f:
        for line in f:
            if line[:1] == b">":
                if pid is not None: b.append((pid, b"".join(seq)))
                if len(b) >= a.batch: yield b; b = []
                pid, seq = line[1:].split()[0].decode(), []
            else:
                seq.append(line.strip())
    if pid is not None: b.append((pid, b"".join(seq)))
    if b: yield b

all_ids = []
if a.maps_only:
    with open(a.out + ".fsb", "rb") as f:
        f.read(4); d, n = struct.unpack("<iq", f.read(12)); f.seek(16 + 4 * d * n)
        all_ids = f.read().decode().split("\n")
else:
  with open(a.out + ".fsb", "wb") as o, Pool(a.workers) as pool:
      o.write(b"FSB1" + struct.pack("<iq", D, 0))
      for ids, X in pool.imap(work, batches()):
          o.write(X.tobytes()); all_ids += ids
          if len(all_ids) % 1000000 < len(ids): print(f"{len(all_ids):,} pieces", flush=True)
      o.write("\n".join(all_ids).encode())
      o.seek(8); o.write(struct.pack("<q", len(all_ids)))
  print(f"{len(all_ids):,} pieces of {a.len} bp ({D} canonical {k}-mers) -> {a.out}.fsb")

for spec in a.layout:
    name, _, path = spec.partition("=")
    pos = {}
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p[0] != "id": pos[p[0]] = (p[1], p[2])
    miss = 0
    with open(f"{a.out}.{name}.tsv", "w") as o:
        o.write("id\tlon\tlat\n")
        for pid in all_ids:
            q = pos.get(pid.rsplit("@", 1)[0])
            if q is None: miss += 1; continue
            o.write(f"{pid}\t{q[0]}\t{q[1]}\n")
    print(f"{name}: {len(all_ids) - miss:,} pieces placed from {path}" + (f", {miss} parents missing" if miss else ""))
