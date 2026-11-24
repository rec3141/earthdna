#!/usr/bin/env python3
"""GC content of each row of k-mer count tables (.fsb from kmers.py / gtdb_prep.py), for pack_globe.py --color.

A k-mer and its reverse complement hold the same number of G+C, so this works on plain (4^k columns) and canonical
tables alike. Writes "id<TAB>gc" (a fraction) with an "id gc" header; --only keeps the ids in a tax table's first
column (e.g. the globe's point set), which is what makes 35M GTDB fragments a 3.9M-row file.

usage: fsb_gc.py OUT.tsv TABLE.fsb [TABLE2.fsb ...] [--only GLOBE.tax.tsv]      ('@list.txt' names a file of paths)
"""
import argparse, struct
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("out"); ap.add_argument("fsb", nargs="+"); ap.add_argument("--only")
a = ap.parse_args()
paths = [q for p in a.fsb for q in ([l.strip() for l in open(p[1:]) if l.strip()] if p.startswith("@") else [p])]
only = {l.split("\t", 1)[0] for l in list(open(a.only))[1:]} if a.only else None

def gc_weights(D):
    k = next(k for k in range(1, 12) if D in (4 ** k, len({min(w, rcw(w, k)) for w in range(4 ** k)})))
    w = np.arange(4 ** k)
    gc = sum(np.isin((w >> (2 * j)) & 3, (1, 2)) for j in range(k))               # A C G T = 0 1 2 3
    if D != 4 ** k: gc = gc[np.unique(np.minimum(w, [rcw(x, k) for x in w]))]     # canonical columns, in index order
    return gc.astype(np.float32) / k

def rcw(w, k): return sum((3 - ((w >> (2 * j)) & 3)) * 4 ** (k - 1 - j) for j in range(k))

n = 0
with open(a.out, "w") as o:
    o.write("id\tgc\n")
    for p in paths:
        with open(p, "rb") as f:
            vt = {b"FSB1": np.float32, b"FSB2": np.uint16}[f.read(4)]   # float32, or uint16 counts
            vb = np.dtype(vt).itemsize
            D, N = struct.unpack("<iq", f.read(12))
            wt = gc_weights(D)
            f.seek(16 + vb * D * N); ids = f.read().decode().split("\n")
            f.seek(16)
            for s in range(0, N, 200000):
                X = np.frombuffer(f.read(vb * D * min(200000, N - s)), vt).reshape(-1, D).astype(np.float32)
                g = X @ wt / np.maximum(X.sum(1), 1)
                for i, v in zip(ids[s:s + len(g)], g):
                    if only is None or i in only: o.write(f"{i}\t{v:.4f}\n"); n += 1
print(f"{n} rows -> {a.out}")
