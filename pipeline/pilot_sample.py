#!/usr/bin/env python3
"""Genus-stratified subset of a SILVA build for comparing feature sets (k-mers vs Evo 2 ...).

Takes up to --per sequences from each of --genera random genera that have at least --per sequences, so every
sequence has same-genus neighbours to find. Writes OUT.fasta and OUT.tax.tsv (rows of the build's tax table).

usage: pilot_sample.py nr99.tax.tsv SILVA_..._trunc.fasta.gz OUT [--genera 2000 --per 10 --seed 1]
"""
import argparse, gzip, random
from collections import defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("tax"); ap.add_argument("fasta"); ap.add_argument("out")
ap.add_argument("--genera", type=int, default=2000); ap.add_argument("--per", type=int, default=10)
ap.add_argument("--seed", type=int, default=1)
a = ap.parse_args()
lines = open(a.tax).read().splitlines()
head, rows = lines[0], [l.split("\t") for l in lines[1:]]
by = defaultdict(list)
for r in rows:
    if r[7]: by[(r[1], r[7])].append(r)
rng = random.Random(a.seed)
ok = sorted(g for g, v in by.items() if len(v) >= a.per)
pick = rng.sample(ok, min(a.genera, len(ok)))
keep = {r[0]: r for g in pick for r in rng.sample(by[g], a.per)}
with open(a.out + ".tax.tsv", "w") as t, open(a.out + ".fasta", "w") as o, gzip.open(a.fasta, "rt") as f:
    t.write(head + "\n")
    on = False
    for line in f:
        if line[0] == ">":
            sid = line[1:].split()[0]
            on = sid in keep
            if on: o.write(f">{sid}\n"); t.write("\t".join(keep[sid]) + "\n")
        elif on:
            o.write(line)
print(f"{len(keep)} sequences from {len(pick)} of {len(ok)} genera with >= {a.per} -> {a.out}.fasta")
