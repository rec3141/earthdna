#!/usr/bin/env python3
"""FASTA (plain or .gz) -> k-mer count table as a fastsne .fsb file (float32 rows, ids after the data).
K-mers spanning ambiguity codes are skipped; U counts as T. --rc adds the reverse complement's counts, making
the profile strand-independent (for reads, which come from either strand). --canonical does the same in half the
columns: a k-mer and its reverse complement share one (512 for k = 5, the GTDB globe's features, gtdb_prep.py).

usage: kmers.py in.fasta[.gz] out.fsb [--k 5] [--rc | --canonical] [--min-len N]
"""
import argparse, gzip, struct
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("fasta"); ap.add_argument("out"); ap.add_argument("--k", type=int, default=5)
ap.add_argument("--rc", action="store_true")
ap.add_argument("--canonical", action="store_true")
ap.add_argument("--min-len", type=int, default=0, help="skip shorter sequences")
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

D = 4 ** k
if a.canonical:
    idx = np.arange(4 ** k)
    rcx = ((3 - ((idx[:, None] // pw[None, :]) % 4)[:, ::-1]) * pw[None, :]).sum(1)
    canon = np.unique(np.minimum(idx, rcx))
    col = np.zeros(4 ** k, np.int64); col[canon] = np.arange(len(canon)); col[rcx[canon]] = np.arange(len(canon))
    D, full = len(canon), counts
    counts = lambda s: np.bincount(col, weights=full(s), minlength=D).astype(np.float32)
elif a.rc:
    fwd, comp = counts, bytes.maketrans(b"ACGTUacgtu", b"TGCAAtgcaa")
    counts = lambda s: fwd(s) + fwd(s.translate(comp)[::-1])

op = gzip.open if a.fasta.endswith(".gz") else open
ids, sid, seq = [], None, None
def put(o):
    s = b"".join(seq)
    if len(s) >= a.min_len: o.write(counts(s).tobytes()); ids.append(sid)
with op(a.fasta, "rb") as f, open(a.out, "wb") as o:
    o.write(b"FSB1" + struct.pack("<iq", D, 0))
    for line in f:
        if line[:1] == b">":
            if seq is not None: put(o)
            sid = line[1:].split()[0].decode(); seq = []
        else:
            seq.append(line.strip())
    if seq is not None: put(o)
    o.write("\n".join(ids).encode())
    o.seek(8); o.write(struct.pack("<q", len(ids)))
print(f"{len(ids)} sequences -> {a.out}")
