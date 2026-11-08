#!/usr/bin/env python3
"""Placed sequences (fastsne --project output, one file per layout) -> an overlay the globe viewer draws on top.

  overlay_NAME.bin   per layout, in --layout order: lon then lat as uint16 (the globe's encoding), all points
  overlays.json      [{name, label, colour, n, layouts, file}], updated in place (one entry per NAME)

Rows must be in the same order in every layout file (fastsne --project keeps input order). Points with nn_dist above
--max-nn are dropped (far from every reference: placement means little).

usage: pack_overlay.py GLOBE_DIR NAME --label "Freshwater · FBA28698 bc16" --colour "#ffd166" \
           --layout sculpt=placed.sculpt.tsv --layout similarity=placed.similarity.tsv [--max-nn 0.5]
"""
import argparse, json, os
import numpy as np, pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("globe"); ap.add_argument("name")
ap.add_argument("--label", required=True); ap.add_argument("--colour", default="#ffd166")
ap.add_argument("--layout", action="append", required=True); ap.add_argument("--max-nn", type=float, default=0)
a = ap.parse_args()

names, cols, ids, keep = [], [], None, None
for spec in a.layout:
    name, path = spec.split("=", 1)
    d = pd.read_csv(path, sep="\t")
    if ids is None:
        ids = d.iloc[:, 0].to_numpy()
        keep = np.ones(len(d), bool) if not a.max_nn else (d.iloc[:, 3].to_numpy() <= a.max_nn)
    elif not np.array_equal(ids, d.iloc[:, 0].to_numpy()):
        raise SystemExit(f"{path}: rows differ from the first layout")
    lon, lat = d.iloc[:, 1].to_numpy()[keep], d.iloc[:, 2].to_numpy()[keep]
    cols += [np.round((lon + 180) / 360 * 65535).astype("<u2"), np.round((lat + 90) / 180 * 65535).astype("<u2")]
    names.append(name)
f = f"overlay_{a.name}.bin"
open(os.path.join(a.globe, f), "wb").write(b"".join(c.tobytes() for c in cols))
jp = os.path.join(a.globe, "overlays.json")
ov = [o for o in (json.load(open(jp)) if os.path.exists(jp) else []) if o["name"] != a.name]
ov.append(dict(name=a.name, label=a.label, colour=a.colour, n=int(keep.sum()), layouts=names, file=f))
json.dump(ov, open(jp, "w"), indent=1)
print(f"{a.name}: {keep.sum():,} of {len(keep):,} points, layouts {names} -> {f}")
