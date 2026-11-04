#!/usr/bin/env python3
"""Place amplicon sequence variants (any SSU region, any primers) on a fastsne globe of full-length references.

1. Map every ASV (both strands) to a panel of aligned references with vsearch and read off the alignment columns
   it covers, so no primer list is needed and mixed runs (e.g. 16S V4 + 18S V9) sort themselves. At least 80% of the
   ASV must align, which keeps out non-SSU amplicons (ITS, 28S LSU) that share short conserved stretches.
2. Group the column ranges into regions (ranges seen in at least --min-asvs ASVs); each ASV joins the region its
   range overlaps best (Jaccard >= --min-overlap), the rest stay unplaced (e.g. ITS, or nothing similar enough).
3. Per region, cut every reference to those columns (one pass over the alignment), count 5-mers and build a
   placement model whose positions come from the full-length map (fastsne --fixed-y).
4. Project each region's ASVs (forward strand) onto the globe: fastsne --project.

usage: place_asvs.py asvs.fasta SILVA_..._full_align_trunc.fasta.gz map.tsv WORKDIR [--fastsne PATH] [--radius 35]
output: WORKDIR/placed.tsv  id, lon, lat, nn_dist, region, strand
"""
import argparse, collections, gzip, os, subprocess
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("asvs"); ap.add_argument("aln"); ap.add_argument("map"); ap.add_argument("work")
ap.add_argument("--fastsne", default="fastsne"); ap.add_argument("--radius", type=float, default=35)
ap.add_argument("--panel-every", type=int, default=40, help="every Nth aligned reference goes in the panel")
ap.add_argument("--min-asvs", type=int, default=20); ap.add_argument("--min-overlap", type=float, default=0.6)
ap.add_argument("--threads", type=int, default=8)
a = ap.parse_args()
os.makedirs(a.work, exist_ok=True)
W = lambda f: os.path.join(a.work, f)
here = os.path.dirname(os.path.abspath(__file__))

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

# 1. panel and column ranges
if not os.path.exists(W("ranges.tsv")):
    cols = {}
    with open(W("panel.fa"), "w") as out:
        for n, (name, s) in enumerate(aligned()):
            if n % a.panel_every: continue
            c = base_cols(s); cols[name] = c
            out.write(f">{name}\n{''.join(s[i] for i in c).upper().replace('U', 'T')}\n")
    subprocess.run(["vsearch", "--usearch_global", a.asvs, "--db", W("panel.fa"), "--id", "0.6", "--strand", "both",
                    "--maxaccepts", "4", "--maxrejects", "64", "--top_hits_only", "--query_cov", "0.8", "--threads", str(a.threads),
                    "--userout", W("panel_hits.tsv"), "--userfields", "query+target+id+qstrand+tilo+tihi", "--quiet"],
                   check=True)
    best = {}
    for line in open(W("panel_hits.tsv")):
        q, t, i, st, lo, hi = line.split("\t")
        if q not in best or float(i) > best[q][0]: best[q] = (float(i), t, st, int(lo), int(hi))
    with open(W("ranges.tsv"), "w") as f:
        for q, (i, t, st, lo, hi) in best.items():
            c = cols[t]
            f.write(f"{q}\t{st}\t{c[lo - 1]}\t{c[hi - 1]}\t{i}\n")

ranges = {}
for line in open(W("ranges.tsv")):
    q, st, c0, c1, _ = line.rstrip("\n").split("\t")
    ranges[q] = (st, int(c0), int(c1))

# 2. regions: frequent column ranges (in 250-column bins), then each ASV to its best-overlapping region
bins = collections.Counter((c0 // 250, c1 // 250) for _, c0, c1 in ranges.values())
regions = []
for (b0, b1), n in bins.most_common():
    if n < a.min_asvs: break
    members = [(c0, c1) for _, c0, c1 in ranges.values() if (c0 // 250, c1 // 250) == (b0, b1)]
    r = (int(np.median([m[0] for m in members])), int(np.median([m[1] for m in members])))
    if all(min(r[1], o[1]) - max(r[0], o[0]) < 0.9 * (r[1] - r[0]) for o in regions): regions.append(r)
def jacc(x, y): return max(0, min(x[1], y[1]) - max(x[0], y[0])) / (max(x[1], y[1]) - min(x[0], y[0]))
assign = {}
for q, (st, c0, c1) in ranges.items():
    j = [jacc((c0, c1), r) for r in regions]
    if j and max(j) >= a.min_overlap: assign[q] = int(np.argmax(j))
cnt = collections.Counter(assign.values())
for k, r in enumerate(regions): print(f"region {k}: alignment columns {r[0]}-{r[1]}, {cnt[k]} ASVs")
print(f"{len(assign)} of {len(ranges)} mapped ASVs assigned to a region")

# 3. reference sets cut to each region, in one pass over the alignment
todo = [k for k in range(len(regions)) if not os.path.exists(W(f"region{k}.model"))]
if todo:
    outs = {k: open(W(f"region{k}_refs.fa"), "w") for k in todo}
    for name, s in aligned():
        c = base_cols(s)
        if len(c) == 0: continue
        for k in todo:
            c0, c1 = regions[k]
            if c[0] > c0 or c[-1] < c1: continue
            reg = s[c0:c1 + 1].replace("-", "").replace(".", "")
            if len(reg) >= 50: outs[k].write(f">{name}\n{reg}\n")
    for f in outs.values(): f.close()
    for k in todo:
        subprocess.run(["python3", os.path.join(here, "kmers.py"), W(f"region{k}_refs.fa"), W(f"region{k}_refs.fsb")], check=True)
        subprocess.run([a.fastsne, "-i", W(f"region{k}_refs.fsb"), "-o", "/dev/null", "--transform", "sqrt",
                        "--pca", "128", "--sphere", str(a.radius), "--fixed-y", a.map,
                        "--save-model", W(f"region{k}.model")], check=True, stderr=subprocess.DEVNULL)
        os.remove(W(f"region{k}_refs.fsb"))

# 4. project each region's ASVs, forward strand
rc = str.maketrans("ACGTURYKMBDHVN", "TGCAAYRMKVHDBN")
seqs, name = {}, None
for line in open(a.asvs):
    if line[0] == ">": name = line[1:].split()[0]; seqs[name] = []
    else: seqs[name].append(line.strip())
with open(W("placed.tsv"), "w") as out:
    out.write("id\tlon\tlat\tnn_dist\tregion\tstrand\n")
    for k in range(len(regions)):
        qs = [q for q, r in assign.items() if r == k]
        if not qs: continue
        with open(W(f"region{k}_asvs.fa"), "w") as f:
            for q in qs:
                s = "".join(seqs[q]).upper()
                if ranges[q][0] == "-": s = s.translate(rc)[::-1]
                f.write(f">{q}\n{s}\n")
        subprocess.run(["python3", os.path.join(here, "kmers.py"), W(f"region{k}_asvs.fa"), W(f"region{k}_asvs.fsb")], check=True)
        subprocess.run([a.fastsne, "--project", W(f"region{k}.model"), "-i", W(f"region{k}_asvs.fsb"),
                        "-o", W(f"region{k}_placed.tsv")], check=True, stderr=subprocess.DEVNULL)
        for line in list(open(W(f"region{k}_placed.tsv")))[1:]:
            q, lo, la, nd = line.rstrip("\n").split("\t")
            out.write(f"{q}\t{lo}\t{la}\t{nd}\t{regions[k][0]}-{regions[k][1]}\t{ranges[q][0]}\n")
print("wrote", W("placed.tsv"))
