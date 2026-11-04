#!/usr/bin/env python3
"""FASTA (plain or .gz) -> k-mer count table as a fastsne .fsb file (float32 rows, ids after the data).
K-mers spanning ambiguity codes are skipped; U counts as T.

usage: kmers.py in.fasta[.gz] out.fsb [--k 5]
"""
import argparse, gzip, struct
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("fasta"); ap.add_argument("out"); ap.add_argument("--k", type=int, default=5)
a = ap.parse_args()
k = a.k
lut = np.full(256, 255, np.uint8)
for i, c in enumerate(b"ACGT"): lut[c] = i; lut[c + 32] = i
lut[ord("U")] = lut[ord("u")] = 3
pw = (4 ** np.arange(k - 1, -1, -1)).astype(np.int64)

def counts(s):
    c = lut[np.frombuffer(s, np.uint8)]
    bad = c > 3
    c = np.where(bad, 0, c).astype(np.int64)
    if len(c) < k: return np.zeros(4 ** k, np.float32)
    w = np.lib.stride_tricks.sliding_window_view(c, k) @ pw
    ok = ~np.lib.stride_tricks.sliding_window_view(bad, k).any(1)
    return np.bincount(w[ok], minlength=4 ** k).astype(np.float32)

op = gzip.open if a.fasta.endswith(".gz") else open
ids, seq = [], None
with op(a.fasta, "rb") as f, open(a.out, "wb") as o:
    o.write(b"FSB1" + struct.pack("<iq", 4 ** k, 0))
    for line in f:
        if line[:1] == b">":
            if seq is not None: o.write(counts(b"".join(seq)).tobytes())
            ids.append(line[1:].split()[0].decode()); seq = []
        else:
            seq.append(line.strip())
    if seq is not None: o.write(counts(b"".join(seq)).tobytes())
    o.write("\n".join(ids).encode())
    o.seek(8); o.write(struct.pack("<q", len(ids)))
print(f"{len(ids)} sequences -> {a.out}")
