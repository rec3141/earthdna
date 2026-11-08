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
rows = [l.rstrip("\n").split("\t") for l in open(a.tax)][1:]
taxid = {r[0]: r for r in rows}
ranks = {"domain": 1, "phylum": 3, "class": 4, "order": 5, "family": 6, "genus": 7}

def read_fsb(path):
    with open(path, "rb") as f:
        magic, d, n = f.read(4), *struct.unpack("<iq", f.read(12))
        X = np.frombuffer(f.read(4 * d * n), np.float32).reshape(n, d)
        ids = f.read().decode().split("\n")
    return X, ids

print("table\t" + "\t".join(ranks))
for spec in a.tables:
    path, _, tr = spec.partition(":")
    X, ids = read_fsb(path)
    sel = [i for i, s in enumerate(ids) if s in taxid]
    X = X[sel].astype(np.float64)
    lab = [taxid[ids[i]] for i in sel]
    if tr == "sqrt": X = np.sqrt(X)
    if tr == "z": X = (X - X.mean(0)) / (X.std(0) + 1e-9)   # standardise channels (LM activations)
    X -= X.mean(0)
    if a.pca and a.pca < X.shape[1]:
        _, _, Vt = np.linalg.svd(X[np.random.default_rng(0).choice(len(X), min(len(X), 50000), False)],
                                 full_matrices=False)
        X = X @ Vt[:a.pca].T
    X = X.astype(np.float32)
    sq = (X * X).sum(1)
    nn = np.zeros((len(X), a.k), np.int64)
    for s in range(0, len(X), 2048):
        d = sq[s:s + 2048, None] + sq[None, :] - 2 * X[s:s + 2048] @ X.T
        d[np.arange(d.shape[0]), np.arange(s, s + d.shape[0])] = np.inf
        nn[s:s + 2048] = np.argpartition(d, a.k, 1)[:, :a.k]
    out = []
    for r, c in ranks.items():
        L = np.array([x[c] for x in lab], dtype=object)
        ok = L != ""
        out.append(f"{(L[nn][ok] == L[ok, None]).mean():.3f}")
    print(spec + "\t" + "\t".join(out), flush=True)
