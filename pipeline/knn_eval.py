#!/usr/bin/env python3
"""How often a sequence's nearest neighbours share its label, per rank, for one or more feature tables.

Each table goes through fastsne's own preprocessing (optional sqrt, centre, PCA to --pca components), then exact
Euclidean k-NN. Prints, per table, the mean fraction of the k neighbours that share the query's
domain / phylum / class / order / family / genus (queries with an empty label at that rank are skipped).

usage: knn_eval.py OUT.tax.tsv a.fsb[:sqrt|:z] b.fsb ... [--k 15 --pca 128]
"""
import argparse, struct
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("tax"); ap.add_argument("tables", nargs="+")
ap.add_argument("--k", type=int, default=15); ap.add_argument("--pca", type=int, default=128)
a = ap.parse_args()
with open(a.tax) as f:
    rows = [l.rstrip("\n").split("\t") for l in f][1:]
taxid = {r[0]: r for r in rows}
ranks = {"domain": 1, "phylum": 3, "class": 4, "order": 5, "family": 6, "genus": 7}

def read_fsb(path):
    with open(path, "rb") as f:
        hdr = f.read(16)
        if len(hdr) != 16 or hdr[:4] != b"FSB1":
            raise SystemExit(f"{path}: invalid FSB1 header")
        d, n = struct.unpack("<iq", hdr[4:])
        if d < 0 or n < 0:
            raise SystemExit(f"{path}: negative dimensions ({n}, {d})")
        mbytes = f.read(4 * d * n)
        if len(mbytes) != 4 * d * n:
            raise SystemExit(f"{path}: truncated FSB1 matrix ({len(mbytes)} of {4 * d * n} bytes)")
        X = np.frombuffer(mbytes, np.float32).reshape(n, d)
        tail = f.read().decode()
        ids = tail.splitlines() if n else ([] if tail == "" else tail.splitlines())
        if len(ids) != n:
            raise SystemExit(f"{path}: id count {len(ids)} != header n {n}")
    return X, ids

print("table\t" + "\t".join(ranks))
for spec in a.tables:
    path, _, tr = spec.partition(":")
    if tr not in ("", "sqrt", "z"):
        raise SystemExit(f"{spec}: unknown transform {tr!r} (expected '', 'sqrt', or 'z')")
    X, ids = read_fsb(path)
    sel = [i for i, s in enumerate(ids) if s in taxid]
    if len(sel) <= a.k:
        raise SystemExit(f"{spec}: matched {len(sel)} taxonomy rows, need > k ({a.k})")
    X = X[sel].astype(np.float64)
    if not np.isfinite(X).all():
        raise SystemExit(f"{spec}: non-finite values in input matrix")
    lab = [taxid[ids[i]] for i in sel]
    if tr == "sqrt":
        if (X < 0).any():
            raise SystemExit(f"{spec}: negative values cannot use :sqrt")
        X = np.sqrt(X)
    if tr == "z": X = (X - X.mean(0)) / (X.std(0) + 1e-9)   # standardise channels (LM activations)
    X -= X.mean(0)
    if a.pca and a.pca < X.shape[1]:
        _, _, Vt = np.linalg.svd(X[np.random.default_rng(0).choice(len(X), min(len(X), 50000), False)],
                                 full_matrices=False)
        X = X @ Vt[:a.pca].T
    sq = (X * X).sum(1)
    nn = np.zeros((len(X), a.k), np.int64)
    for s in range(0, len(X), 2048):
        d = sq[s:s + 2048, None] + sq[None, :] - 2 * (X[s:s + 2048] @ X.T)
        d[np.arange(d.shape[0]), np.arange(s, s + d.shape[0])] = np.inf
        nn[s:s + 2048] = np.argpartition(d, a.k, 1)[:, :a.k]
    out = []
    for r, c in ranks.items():
        L = np.array([x[c] for x in lab], dtype=object)
        ok = L != ""
        out.append(f"{(L[nn][ok] == L[ok, None]).mean():.3f}")
    print(spec + "\t" + "\t".join(out), flush=True)
