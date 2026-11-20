#!/usr/bin/env python3
"""Pack fastsne sphere layouts + taxonomy into the files viewer/globe.html loads (all served as text, since
artifact hosts and some static hosts refuse binary types):

  globe.json        names (groups, genera, lineages, species), counts, layouts
  pos.b64.txt       first layout lon, lat (uint16 each) + group (uint8), base64
  tax.b64.txt       genus, lineage (uint16; 65535 = none), species codes (uint16, or uint32 with 2^32-1 = none
                    when a release has more species than fit; globe.json speciesBytes says which), base64
  layout_NAME.b64.txt  further layouts: lon, lat
  flags.b64.txt     conflict flag (uint8) + neighbours' majority group (uint8; 255 = none), base64
  ranks.b64.txt     phylum, class, order, family codes (uint16 each; 65535 = none), base64, for drilling down
                    rank by rank (domain comes from the group, genus and species from tax.b64.txt)
  traits.b64.txt    with --traits: one byte per point, bit k set when trait k applies (globe.json traits names them)
  color_NAME.b64.txt  with --color: one byte per point for a colour scheme (globe.json colors describes each): a number
                    scaled so 0..254 spans its min..max, or a category code; 255 = no value

Rows are shuffled once (fixed seed) so points draw in random order.

usage: pack_globe.py OUT.tax.tsv OUTDIR --layout similarity=map.tsv [--layout sculpt=map2.tsv] [--flags conflicts.tsv]
                     [--traits traits.tsv] [--color gc=gc.tsv ...] [--title "SILVA 144 SSU Ref NR99"]
                     [--layout-label "NAME=Button text|Status line"] [--about "subtitle"]
                     [--labels "KEY=Display name=labels.tsv"] [--labels-info KEY=info.tsv] [--labels-tip "KEY=Mostly in"]
                     [--default-labels KEY] [--geo-layout NAME]
  The first --layout is the one the viewer opens with. --layout-label names a layout's button and the line shown when
  it is chosen (the viewer has defaults for similarity and sculpt); --about replaces the viewer's subtitle, which
  otherwise describes an SSU globe.
  --labels adds a label set the viewer can switch to in place of the taxonomy (e.g. environment): labels.tsv has a
  header naming a column of OUT.tax.tsv and then one column per level, top level first, and a row per key; an empty
  value means no label at that level. Written as labels_KEY.b64.txt, one uint16 array per level (65535 = none),
  named in globe.json labelSets. --labels-info gives a label set a tooltip: info.tsv has a header naming a column of
  OUT.tax.tsv and then the columns to show (the first is the title; lat and lon are shown as coordinates), written as
  info_KEY.txt, one tab-separated line per point; --labels-tip words go before the label path in that tooltip.
  --default-labels opens the globe on that label set instead of the taxonomy. --geo-layout marks a layout whose
  lon/lat are places on Earth, which the viewer must show as they are (it otherwise turns each layout to suit it).
  --color NAME=values.tsv: a header naming a column of OUT.tax.tsv (id, organism, ...) and one value column, then a
  row per key (fsb_gc.py writes id/gc). All-numeric values make a continuous scheme, anything else a categorical
  one (at most 255 categories, the commonest first). The viewer's schemes say how each is coloured.
  traits.tsv: a header naming a column of OUT.tax.tsv (e.g. organism) and up to 8 traits, then one row per key with
  0/1 per trait (gtdb_traits.py writes one from GTDB's metadata); points whose key is missing get no traits
"""
import argparse, base64, collections, json, os
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("tax"); ap.add_argument("outdir")
ap.add_argument("--layout", action="append", required=True)
ap.add_argument("--flags"); ap.add_argument("--title", default="")
ap.add_argument("--traits"); ap.add_argument("--color", action="append", default=[])
ap.add_argument("--layout-label", action="append", default=[]); ap.add_argument("--about", default="")
ap.add_argument("--labels", action="append", default=[])
ap.add_argument("--labels-info", action="append", default=[]); ap.add_argument("--labels-tip", action="append", default=[])
ap.add_argument("--default-labels", default="")
ap.add_argument("--geo-layout", action="append", default=[])
ap.add_argument("--models", help="models globe: KEY=label;... (layouts named <shape>[-KEY]; the first KEY is empty), e.g. \"=Evo 2 7B;kmer=5-mers\"")
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)

with open(a.tax) as tf:
    lines = [l.rstrip("\n").split("\t") for l in tf]
header, rows = lines[0], lines[1:]
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
sp_values = [r[C["species"]] for r in rows]
species_bytes = 2 if len(set(v for v in sp_values if v)) < 65535 else 4
species, _, sp = codes(sp_values, 65535 if species_bytes == 2 else 2 ** 32 - 1)
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
        with open(os.path.join(a.outdir, "pos.b64.txt"), "w") as f:
            f.write(b64(lo, la, gi[perm].astype(np.uint8)))
    else:
        with open(os.path.join(a.outdir, f"layout_{name}.b64.txt"), "w") as f:
            f.write(b64(lo, la))
    layouts.append(name)
with open(os.path.join(a.outdir, "tax.b64.txt"), "w") as f:
    f.write(b64(ge[perm].astype("<u2"), li[perm].astype("<u2"), sp[perm].astype("<u2" if species_bytes == 2 else "<u4")))
with open(os.path.join(a.outdir, "ranks.b64.txt"), "w") as f:
    f.write(b64(*(x[perm].astype("<u2") for x in rank_codes)))
if a.flags:
    with open(a.flags) as ff:
        fr = {l.split("\t")[0]: l.rstrip("\n").split("\t") for l in list(ff)[1:]}
    gix = {g: i for i, g in enumerate(groups)}
    flag = np.array([int(fr[r[0]][1]) for r in rows], np.uint8)
    def maj_group(rid):
        rec = fr[rid]
        if not int(rec[1]):
            return 255
        grp = rec[3] if len(rec) > 3 else ""
        return gix.get(rec[2] + " · " + grp, 255) if grp else 255
    maj = np.array([maj_group(r[0]) for r in rows], np.uint8)
    with open(os.path.join(a.outdir, "flags.b64.txt"), "w") as f:
        f.write(b64(flag[perm], maj[perm]))
trait_names = []
if a.traits:
    with open(a.traits) as trf:
        tl = [l.rstrip("\n").split("\t") for l in trf]
    key, trait_names = tl[0][0], tl[0][1:]
    assert len(trait_names) <= 8, "at most 8 traits fit in a byte"
    kc = header.index(key)
    bits = {r[0]: sum(1 << k for k, v in enumerate(r[1:]) if v == "1") for r in tl[1:]}
    tr = np.array([bits.get(r[kc], 0) for r in rows], np.uint8)
    with open(os.path.join(a.outdir, "traits.b64.txt"), "w") as f:
        f.write(b64(tr[perm]))
    print("traits: " + ", ".join(f"{n} {int((tr >> k & 1).sum())}" for k, n in enumerate(trait_names)))
colors = []
for spec in a.color:
    name, path = spec.split("=", 1)
    with open(path) as cf:
        cl = [l.rstrip("\n").split("\t") for l in cf]
    kc, val = header.index(cl[0][0]), {r[0]: r[1] for r in cl[1:] if len(r) > 1 and r[1] != ""}
    v = [val.get(r[kc]) for r in rows]
    try:
        x = np.array([np.nan if t is None else float(t) for t in v])
        lo, hi = float(np.nanmin(x)), float(np.nanmax(x))
        code = np.where(np.isnan(x), 255, np.round((np.nan_to_num(x, nan=lo) - lo) / max(hi - lo, 1e-12) * 254)).astype(np.uint8)
        colors.append(dict(name=name, kind="continuous", min=lo, max=hi))
    except ValueError:
        labels = [k for k, _ in collections.Counter(t for t in v if t is not None).most_common(255)]
        ix = {k: i for i, k in enumerate(labels)}
        code = np.array([ix.get(t, 255) for t in v], np.uint8)
        colors.append(dict(name=name, kind="categorical", labels=labels))
    with open(os.path.join(a.outdir, f"color_{name}.b64.txt"), "w") as f:
        f.write(b64(code[perm]))
    print(f"color {name}: {colors[-1]['kind']}, {int((code != 255).sum())} of {N} points have a value")
label_sets = []
for spec in a.labels:
    key, name, path = spec.split("=", 2)
    with open(path) as lf:
        ll = [l.rstrip("\n").split("\t") for l in lf]
    kc, levels = header.index(ll[0][0]), ll[0][1:]
    by = {r[0]: r[1:] for r in ll[1:]}
    names, arrs = {}, []
    for k, lv in enumerate(levels):
        vals = [(by.get(r[kc]) or [""] * len(levels))[k] for r in rows]
        names[lv], _, c = codes(vals)
        arrs.append(c[perm].astype("<u2"))
    with open(os.path.join(a.outdir, f"labels_{key}.b64.txt"), "w") as f:
        f.write(b64(*arrs))
    label_sets.append(dict(key=key, name=name, levels=levels, names=names))
    print(f"labels {key}: " + ", ".join(f"{lv} {len(names[lv])}" for lv in levels) + f"; {sum(1 for r in rows if r[kc] in by)} of {N} points keyed")
for spec in a.labels_info:
    key, path = spec.split("=", 1)
    ls = next(x for x in label_sets if x["key"] == key)
    with open(path) as lf:
        il = [l.rstrip("\n").split("\t") for l in lf]
    kc, cols = header.index(il[0][0]), il[0][1:]
    by = {r[0]: r[1:] for r in il[1:]}
    empty = [""] * len(cols)
    with open(os.path.join(a.outdir, f"info_{key}.txt"), "w") as f:
        f.write("\n".join("\t".join(by.get(rows[i][kc], empty)) for i in perm) + "\n")
    ls["info"] = cols
for spec in a.labels_tip:
    key, text = spec.split("=", 1)
    next(x for x in label_sets if x["key"] == key)["tip"] = text
assert not a.default_labels or any(x["key"] == a.default_labels for x in label_sets), "--default-labels: no such --labels"
layout_info = {}
for spec in a.layout_label:
    name, text = spec.split("=", 1)
    assert name in layouts, f"--layout-label {name}: no such --layout"
    label, _, about = text.partition("|")
    layout_info[name] = dict(label=label, about=about)
extra = dict(models=[dict(key=k, label=v) for k, v in (m.split("=", 1) for m in a.models.split(";"))]) if a.models else {}
with open(os.path.join(a.outdir, "globe.json"), "w") as gf:
    json.dump(dict(**extra, colors=colors, n=N, title=a.title, about=a.about, layouts=layouts, layoutInfo=layout_info, labelSets=label_sets, defaultLabels=a.default_labels, geoLayouts=a.geo_layout, flags=bool(a.flags), traits=trait_names, groups=groups, groupCounts=gcount,
                   genera=genera, lineages=[x[1:] for x in lineages], lineageRank="".join(x[0] for x in lineages),
                   species=species, speciesBytes=species_bytes, ranks=dict(levels=["domain"] + RANKS + ["genus", "species"], names=rank_names)),
              gf, separators=(",", ":"))
print(f"{N} points, {len(groups)} groups, {len(genera)} genera, {len(species)} species -> {a.outdir}")
