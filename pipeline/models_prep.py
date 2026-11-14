#!/usr/bin/env python3
"""Inputs fastsne needs for the DNA-model globe (build_models.sh), from the shared set (models_set.py):

  OUT.anchor              spring groups per window (domain, group, class, order), as gtdb_reselect.py writes them
  OUT.<name>.fsb          per --zscore table: channels standardised (language-model activations have a few huge
                          channels that would otherwise own the PCA)
  OUT.<name>.fsb          per --subspace NAME=features.fsb=model.npz: the 10-d phylogeny subspace (train_subspace's
                          W, b, channel mean / sd), rows normalised (the subspace's distance is an angle)

usage: models_prep.py SET.tax.tsv OUT [--zscore NAME=in.fsb ...] [--subspace NAME=in.fsb=model.npz ...]
"""
import argparse, struct
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("tax"); ap.add_argument("out")
ap.add_argument("--zscore", action="append", default=[]); ap.add_argument("--subspace", action="append", default=[])
a = ap.parse_args()
rows = [l.rstrip("\n").split("\t") for l in open(a.tax)][1:]
ids = [r[0] for r in rows]

def read(p):
    with open(p, "rb") as f:
        magic = f.read(4)
        if magic != b"FSB1":
            raise SystemExit(f"{p}: invalid FSB1 magic {magic!r}")
        d, n = struct.unpack("<iq", f.read(12))
        buf = f.read(4 * d * n)
        if len(buf) != 4 * d * n:
            raise SystemExit(f"{p}: truncated FSB1 matrix ({len(buf)} < {4 * d * n})")
        X = np.frombuffer(buf, np.float32).reshape(n, d)
        tail = f.read().decode()
        got = tail.split("\n")[:n] if n > 0 else ([] if tail == "" else tail.split("\n"))
    assert got == ids, f"{p}: rows are not in {a.tax} order"
    if not np.isfinite(X).all():
        raise SystemExit(f"{p}: non-finite values in input features")
    return X.astype(np.float64)

def write(p, X):
    with open(p, "wb") as o:
        o.write(b"FSB1" + struct.pack("<iq", X.shape[1], len(ids))); o.write(X.astype(np.float32).tobytes())
        o.write("\n".join(ids).encode())
    print(f"{X.shape[0]} x {X.shape[1]} -> {p}")

t = [r[1:] for r in rows]   # domain, group, phylum, class, order, ...
codes = []
for j, deep in ((0, 1), (1, 2), (3, 4), (4, 5)):
    ix = {}; codes.append([ix.setdefault(";".join(r[:deep]), len(ix)) if r[j] else -1 for r in t])
with open(a.out + ".anchor", "w") as f:
    for row in zip(*codes): f.write(" ".join(map(str, row)) + "\n")
for spec in a.zscore:
    name, p = spec.split("=", 1); X = read(p)
    write(f"{a.out}.{name}.fsb", (X - X.mean(0)) / (X.std(0) + 1e-9))
for spec in a.subspace:
    name, p, m = spec.split("="); X = read(p); M = np.load(m)
    sd = M["sd"]
    if not (np.isfinite(sd).all() and np.all(sd > 0)):
        raise SystemExit(f"{m}: non-positive or non-finite channel sd")
    Z = ((X - M["mu"]) / sd - M["b"]) @ M["W"].T
    norms = np.linalg.norm(Z, axis=1, keepdims=True)
    if not (np.isfinite(norms).all() and np.all(norms > 0)):
        raise SystemExit(f"{p}: zero or non-finite subspace row norm")
    write(f"{a.out}.{name}.fsb", Z / norms)
