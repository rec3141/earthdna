#!/usr/bin/env python3
"""Sequences whose label disagrees with their sequence neighbours (likely annotation errors).

A sequence is flagged when at most 20% of its 15 nearest neighbours (fastsne --write-knn graph) share its label and
at least 60% share one other label: domain conflicts first, then group (phylum, or eukaryote major clade).
Writes id, flag (0 none, 1 group, 2 domain), neighbours' majority domain and group.

usage: conflicts.py OUT.tax.tsv graph.knn K out.tsv
"""
import sys
import numpy as np

taxf, knnf, K, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
rows = [l.rstrip("\n").split("\t") for l in open(taxf)][1:]
N = len(rows)
I = np.fromfile(knnf, np.int32).reshape(N, K)[:, :15]

def check(labels):
    names = sorted(set(labels) - {""})
    ix = {n: i for i, n in enumerate(names)}
    L = np.array([ix.get(x, -1) for x in labels])
    NL = L[I]
    same = (NL == L[:, None]).mean(1)
    maj = np.array([np.bincount(r[r >= 0]).argmax() if (r >= 0).any() else -1 for r in NL])
    majf = (NL == maj[:, None]).mean(1)
    return (L >= 0) & (same <= 0.2) & (majf >= 0.6), [names[m] if m >= 0 else "" for m in maj]

dom = [r[1] for r in rows]
cd, md = check(dom)
cg, mg = check([r[1] + ";" + r[2] if r[2] else "" for r in rows])
flag = np.where(cd, 2, np.where(cg, 1, 0))
with open(out, "w") as f:
    f.write("id\tflag\tneighbour_domain\tneighbour_group\n")
    for r, fl, a, b in zip(rows, flag, md, mg):
        d, _, g = b.partition(";")
        f.write(f"{r[0]}\t{fl}\t{d or a}\t{g}\n")
print(f"{(flag == 2).sum()} domain and {(flag == 1).sum()} group conflicts of {N}")
