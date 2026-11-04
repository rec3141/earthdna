#!/usr/bin/env python3
"""Pack fastsne sphere layouts + taxonomy into the files viewer/globe.html loads (all served as text, since
artifact hosts and some static hosts refuse binary types):

  globe.json        names (groups, genera, lineages, species), counts, layouts
  pos.b64.txt       first layout lon, lat (uint16 each) + group (uint8), base64
  tax.b64.txt       genus, lineage, species codes (uint16; 65535 = none), base64
  layout_NAME.b64.txt  further layouts: lon, lat
  flags.b64.txt     conflict flag (uint8) + neighbours' majority group (uint8; 255 = none), base64
  ranks.b64.txt     phylum, class, order, family codes (uint16 each; 65535 = none), base64, for drilling down
                    rank by rank (domain comes from the group, genus and species from tax.b64.txt)

Rows are shuffled once (fixed seed) so points draw in random order.

usage: pack_globe.py OUT.tax.tsv OUTDIR --layout similarity=map.tsv [--layout sculpt=map2.tsv] [--flags conflicts.tsv]
                     [--title "SILVA 144 SSU Ref NR99"]
"""
import argparse, base64, collections, json, os
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("tax"); ap.add_argument("outdir")
ap.add_argument("--layout", action="append", required=True)
ap.add_argument("--flags"); ap.add_argument("--title", default="")
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)

rows = [l.rstrip("\n").split("\t") for l in open(a.tax)][1:]
N = len(rows)
perm = np.random.default_rng(3).permutation(N)
C = dict(dom=1, group=2, phylum=3, **{"class": 4}, order=5, family=6, genus=7, species=8, lrank=9, lineage=10)

def codes(values, limit=65535):
    c = collections.Counter(v for v in values if v)
    names = [k for k, _ in c.most_common()]
    if len(names) >= limit: raise SystemExit(f"{len(names)} names do not fit a {limit}-code field")
    ix = {k: i for i, k in enumerate(names)}
    return names, [c[k] for k in names], np.array([ix.get(v, limit) for v in values])

groups, gcount, gi = codes([r[C["dom"]] + " · " + (r[C["group"]] or "?") for r in rows], 255)
genera, _, ge = codes([r[C["genus"]] for r in rows])
lineages, _, li = codes([r[C["lrank"]][:1] + r[C["lineage"]] if r[C["lineage"]] else "" for r in rows])
species, _, sp = codes([r[C["species"]] for r in rows])
RANKS = ["phylum", "class", "order", "family"]
rank_names, rank_codes = {}, []
for rk in RANKS:
    rank_names[rk], _, c = codes([r[C[rk]] for r in rows])
    rank_codes.append(c)
b64 = lambda *arrs: base64.b64encode(b"".join(np.ascontiguousarray(x).tobytes() for x in arrs)).decode()

def lonlat(path):
    ids = [r[0] for r in rows]
    pos = {}
    with open(path) as f:
        next(f)
        for line in f:
            i, lo, la = line.split("\t")[:3]
            pos[i] = (float(lo), float(la))
    M = np.array([pos[i] for i in ids])[perm]
    return (np.round((M[:, 0] + 180) / 360 * 65535).astype("<u2"), np.round((M[:, 1] + 90) / 180 * 65535).astype("<u2"))

layouts = []
for k, spec in enumerate(a.layout):
    name, path = spec.split("=", 1)
    lo, la = lonlat(path)
    if k == 0:
        open(os.path.join(a.outdir, "pos.b64.txt"), "w").write(b64(lo, la, gi[perm].astype(np.uint8)))
    else:
        open(os.path.join(a.outdir, f"layout_{name}.b64.txt"), "w").write(b64(lo, la))
    layouts.append(name)
open(os.path.join(a.outdir, "tax.b64.txt"), "w").write(b64(*(x[perm].astype("<u2") for x in (ge, li, sp))))
open(os.path.join(a.outdir, "ranks.b64.txt"), "w").write(b64(*(x[perm].astype("<u2") for x in rank_codes)))
if a.flags:
    fr = {l.split("\t")[0]: l.rstrip("\n").split("\t") for l in list(open(a.flags))[1:]}
    gix = {g: i for i, g in enumerate(groups)}
    flag = np.array([int(fr[r[0]][1]) for r in rows], np.uint8)
    maj = np.array([gix.get(fr[r[0]][2] + " · " + fr[r[0]][3], 255) if int(fr[r[0]][1]) else 255 for r in rows], np.uint8)
    open(os.path.join(a.outdir, "flags.b64.txt"), "w").write(b64(flag[perm], maj[perm]))
json.dump(dict(n=N, title=a.title, layouts=layouts, flags=bool(a.flags), groups=groups, groupCounts=gcount,
               genera=genera, lineages=[x[1:] for x in lineages], lineageRank="".join(x[0] for x in lineages),
               species=species, ranks=dict(levels=["domain"] + RANKS + ["genus", "species"], names=rank_names)),
          open(os.path.join(a.outdir, "globe.json"), "w"), separators=(",", ":"))
print(f"{N} points, {len(groups)} groups, {len(genera)} genera, {len(species)} species -> {a.outdir}")
