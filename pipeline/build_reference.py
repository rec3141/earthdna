#!/usr/bin/env python3
"""Build the EarthDNA reference database that amplicon pipelines place ASVs with (danaSeq illumina_amplicon
PLACE_GLOBE, or pipeline/place_asvs.py):

  manifest.json      release, radius, regions, file list
  map.tsv            full-length reference positions (id, lon, lat) from build.sh
  panel.fa           every Nth aligned reference, ungapped, for detecting which region an ASV covers
  panel_cols.npy     their alignment columns, concatenated; panel_offs.npy: row offsets
  regions/C0-C1.model  placement models: references cut to alignment columns C0-C1, 5-mer counts, positions
                     from map.tsv (fastsne --fixed-y)

usage: build_reference.py SILVA_..._full_align_trunc.fasta.gz map.tsv OUTDIR --region C0-C1 [--region ...]
                          [--fastsne fastsne] [--pca 64] [--radius 35] [--release 144] [--panel-every 40]
Region columns are what pipeline/place_asvs.py reports for a set of ASVs (e.g. 14964-22534 for 16S V4).
"""
import argparse, gzip, json, os, shutil, subprocess
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("aln"); ap.add_argument("map"); ap.add_argument("out")
ap.add_argument("--region", action="append", required=True)
ap.add_argument("--fastsne", default="fastsne"); ap.add_argument("--pca", type=int, default=64)
ap.add_argument("--radius", type=float, default=35); ap.add_argument("--release", default="")
ap.add_argument("--panel-every", type=int, default=40)
a = ap.parse_args()
os.makedirs(os.path.join(a.out, "regions"), exist_ok=True)
here = os.path.dirname(os.path.abspath(__file__))
regions = [tuple(int(x) for x in r.split("-")) for r in a.region]

def aligned():
    name, seq = None, []
    with gzip.open(a.aln, "rt") as f:
        for line in f:
            if line[0] == ">":
                if name: yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
    if name: yield name, "".join(seq)

def base_cols(s):
    arr = np.frombuffer(s.encode(), np.uint8)
    return np.where(((arr >= 65) & (arr <= 90)) | ((arr >= 97) & (arr <= 122)))[0].astype(np.int32)

# one pass: panel + every region's reference set
cols, offs = [], [0]
outs = {r: open(os.path.join(a.out, f"regions/{r[0]}-{r[1]}.fa"), "w") for r in regions}
with open(os.path.join(a.out, "panel.fa"), "w") as panel:
    for n, (name, s) in enumerate(aligned()):
        c = base_cols(s)
        if len(c) == 0: continue
        if n % a.panel_every == 0:
            panel.write(f">{name}\n{''.join(s[i] for i in c).upper().replace('U', 'T')}\n")
            cols.append(c); offs.append(offs[-1] + len(c))
        for (c0, c1), f in outs.items():
            if c[0] > c0 or c[-1] < c1: continue
            reg = s[c0:c1 + 1].replace("-", "").replace(".", "")
            if len(reg) >= 50: f.write(f">{name}\n{reg}\n")
for f in outs.values(): f.close()
np.save(os.path.join(a.out, "panel_cols.npy"), np.concatenate(cols))
np.save(os.path.join(a.out, "panel_offs.npy"), np.array(offs, np.int64))
shutil.copy(a.map, os.path.join(a.out, "map.tsv"))

for c0, c1 in regions:
    base = os.path.join(a.out, f"regions/{c0}-{c1}")
    subprocess.run(["python3", os.path.join(here, "kmers.py"), base + ".fa", base + ".fsb"], check=True)
    subprocess.run([a.fastsne, "-i", base + ".fsb", "-o", "/dev/null", "--transform", "sqrt", "--pca", str(a.pca),
                    "--sphere", str(a.radius), "--fixed-y", a.map, "--save-model", base + ".model"],
                   check=True, stderr=subprocess.DEVNULL)
    os.remove(base + ".fsb"); os.remove(base + ".fa")
    print(f"region {c0}-{c1}: {os.path.getsize(base + '.model') / 1e6:.0f} MB")

json.dump(dict(release=a.release, radius=a.radius, pca=a.pca, transform="sqrt", k=5,
               regions=[f"{c0}-{c1}" for c0, c1 in regions], panel_every=a.panel_every),
          open(os.path.join(a.out, "manifest.json"), "w"), indent=1)
print("wrote", a.out)
