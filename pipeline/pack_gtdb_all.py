#!/usr/bin/env python3
"""Every GTDB fragment on one globe: the embedded reference plus all projected fragments, packed straight into the
binary web files without the per-row taxonomy table gtdb_globe_set.py + pack_globe.py go through, which for 35M
fragments needs more memory than a desktop has. Reference fragments keep their embedded positions and conflict flags,
projected ones take their projected positions and flag 0.

Taxonomy and traits belong to genomes (accession = the first two '_' fields of the id), so they are stored once per
genome: records.bin gives each point its genome (uint32), and tax.bin, ranks.bin, traits.bin hold one entry per genome
in pack_globe.py's layout (globe.json `records` = the number of genomes). pos.bin, layout_NAME.bin and flags.bin stay
per point. This halves the download and the viewer's memory against per-point taxonomy.

usage: pack_gtdb_all.py GENOMES.tsv REF.conflicts.tsv TRAITS.tsv OUTDIR \
           --layout similarity=ref.sphere.tsv:all.sphere.tsv [--layout sculpt=ref.sculpt.tsv:all.sculpt.tsv] [--title T]
  every all.*.tsv lists the same ids in the same order
"""
import argparse, collections, json, os
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("genomes"); ap.add_argument("conflicts"); ap.add_argument("traits"); ap.add_argument("outdir")
ap.add_argument("--layout", action="append", required=True); ap.add_argument("--title", default="")
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)

with open(a.genomes) as f:
    h = next(f).rstrip("\n").split("\t")
    G = [l.rstrip("\n").split("\t") for l in f]
gix = {r[0]: k for k, r in enumerate(G)}
col = lambda name: [r[h.index(name)] for r in G]

# pass 1 (first layout's all file): genome of every point, its projected position, and the rows of reference ids
first = a.layout[0].split("=", 1)[1].split(":")
refrow = {}
with open(first[0]) as f:
    next(f)
    for line in f: refrow[line.split("\t", 1)[0]] = -1
with open(first[1]) as f:
    N = sum(1 for _ in f) - 1
gen = np.empty(N, np.int32)
def quant(lo, la): return round((float(lo) + 180) / 360 * 65535), round((float(la) + 90) / 180 * 65535)

def read_layout(spec, k):
    name, (refp, allp) = spec.split("=", 1)[0], spec.split("=", 1)[1].split(":")
    ll = np.empty((N, 2), np.uint16)
    with open(allp) as f:
        next(f)
        for n, line in enumerate(f):
            i, lo, la = line.split("\t", 3)[:3]
            ll[n] = quant(lo, la)
            if k == 0:
                gen[n] = gix["_".join(i.split("_", 2)[:2])]
                if i in refrow: refrow[i] = n
    assert n == N - 1, f"{allp}: {n + 1} rows, expected {N}"
    m = 0
    with open(refp) as f:
        next(f)
        for line in f:
            i, lo, la = line.split("\t", 3)[:3]
            r = refrow[i]
            if r >= 0: ll[r] = quant(lo, la); m += 1
    print(f"{name}: {N} points, {m} at their embedded position", flush=True)
    return name, ll

perm = np.random.default_rng(3).permutation(N)
layouts = []
for k, spec in enumerate(a.layout):
    name, ll = read_layout(spec, k)
    if k == 0:
        lost = sum(1 for r in refrow.values() if r < 0)
        assert lost == 0, f"{lost} reference ids are missing from the all file"
    ll = ll[perm]
    layouts.append(name)
    if k == 0: first_ll = ll
    else: open(os.path.join(a.outdir, f"layout_{name}.bin"), "wb").write(ll[:, 0].tobytes() + ll[:, 1].tobytes())
gen = gen[perm]
pts = np.bincount(gen, minlength=len(G))   # points per genome

def codes(values, limit=65535):
    """per-genome values -> names by point count, counts, and the code of every genome"""
    c = collections.Counter()
    for v, p in zip(values, pts):
        if v and p: c[v] += int(p)
    names = [k for k, _ in c.most_common()]
    if len(names) >= limit: raise SystemExit(f"{len(names)} names do not fit a {limit}-code field")
    ix = {k: i for i, k in enumerate(names)}
    return names, [c[k] for k in names], np.array([ix.get(v, limit) for v in values], np.int64)

groups, gcount, gi = codes([d + " · " + (g or "?") for d, g in zip(col("domain"), col("group"))], 255)
open(os.path.join(a.outdir, "pos.bin"), "wb").write(first_ll[:, 0].tobytes() + first_ll[:, 1].tobytes() + gi[gen].astype(np.uint8).tobytes())
open(os.path.join(a.outdir, "records.bin"), "wb").write(gen.astype("<u4").tobytes())
del first_ll
genera, _, ge = codes(col("genus"))
lineages, _, li = codes([r[:1] + l if l else "" for r, l in zip(col("lineage_rank"), col("lineage"))])
species_bytes = 2 if len(set(v for v in col("species") if v)) < 65535 else 4
species, _, sp = codes(col("species"), 65535 if species_bytes == 2 else 2 ** 32 - 1)
open(os.path.join(a.outdir, "tax.bin"), "wb").write(ge.astype("<u2").tobytes() + li.astype("<u2").tobytes()
                                                   + sp.astype("<u2" if species_bytes == 2 else "<u4").tobytes())
del ge, li, sp
RANKS = ["phylum", "class", "order", "family"]
rank_names = {}
with open(os.path.join(a.outdir, "ranks.bin"), "wb") as o:
    for rk in RANKS:
        rank_names[rk], _, c = codes(col(rk))
        o.write(c.astype("<u2").tobytes())

flag, maj = np.zeros(N, np.uint8), np.full(N, 255, np.uint8)
inv = np.empty(N, np.int64); inv[perm] = np.arange(N)   # row in the all file -> row on the globe
gpos = {g: i for i, g in enumerate(groups)}
with open(a.conflicts) as f:
    next(f)
    for line in f:
        i, fl, nd, ng = line.rstrip("\n").split("\t")[:4]
        r = refrow.get(i, -1)
        if r >= 0 and int(fl):
            flag[inv[r]] = int(fl); maj[inv[r]] = gpos.get(nd + " · " + ng, 255)
open(os.path.join(a.outdir, "flags.bin"), "wb").write(flag.tobytes() + maj.tobytes())

with open(a.traits) as f:
    th = next(f).rstrip("\n").split("\t")
    key, trait_names = th[0], th[1:]
    bits = {r[0]: sum(1 << k for k, v in enumerate(r[1:]) if v == "1") for r in (l.rstrip("\n").split("\t") for l in f)}
tr = np.array([bits.get(v, 0) for v in col(key)], np.uint8)
open(os.path.join(a.outdir, "traits.bin"), "wb").write(tr.tobytes())
print("traits (points): " + ", ".join(f"{n} {int(((tr >> k & 1) * pts).sum())}" for k, n in enumerate(trait_names)))

json.dump(dict(n=N, records=len(G), title=a.title, layouts=layouts, flags=True, traits=trait_names, groups=groups, groupCounts=gcount,
               genera=genera, lineages=[x[1:] for x in lineages], lineageRank="".join(x[0] for x in lineages),
               species=species, speciesBytes=species_bytes, ranks=dict(levels=["domain"] + RANKS + ["genus", "species"], names=rank_names)),
          open(os.path.join(a.outdir, "globe.json"), "w"), separators=(",", ":"))
print(f"{N} points of {len(G)} genomes, {len(groups)} groups, {len(genera)} genera, {len(species)} species -> {a.outdir}")
