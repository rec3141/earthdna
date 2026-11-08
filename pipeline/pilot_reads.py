#!/usr/bin/env python3
"""Genus-stratified set of labelled nanopore reads from the dana archive, for comparing feature sets.

Labels are each read's sendsketch hit (GTDB-style lineage d__;p__;...;g__;s__ with ANI) from the barcode's
dana.duckdb. Uses one random barcode per flowcell (--all: every barcode), reads with ANI >= --min_ani and length in --len, then up to
--per reads per genus, at most --per_sample of them from one sample. Writes OUT.fasta and OUT.tax.tsv with the build's tax columns (group = phylum)
plus ani and sample.

usage: pilot_reads.py /vistara/dana/out_dana_bc OUT [--min_ani 95 --len 3000,5000 --per 10 --per_sample 3 --seed 1 --all]
needs duckdb (python)
"""
import argparse, glob, os, random
from collections import defaultdict
import duckdb

ap = argparse.ArgumentParser()
ap.add_argument("archive"); ap.add_argument("out")
ap.add_argument("--min_ani", type=float, default=95); ap.add_argument("--len", default="3000,5000")
ap.add_argument("--per", type=int, default=10); ap.add_argument("--seed", type=int, default=1)
ap.add_argument("--flowcells", type=int, default=0, help="limit to this many flowcells (0 = all)")
ap.add_argument("--all", action="store_true", help="every barcode, not one random barcode per flowcell")
ap.add_argument("--per_sample", type=int, default=3, help="at most this many reads of a genus from one sample")
a = ap.parse_args()
lo, hi = map(int, a.len.split(","))
rng = random.Random(a.seed)

by_fc = defaultdict(list)
for p in sorted(glob.glob(f"{a.archive}/*/barcode*/dana.duckdb")):
    by_fc[p.split("/")[-3]].append(p)
fcs = sorted(by_fc)
if a.flowcells: fcs = rng.sample(fcs, min(a.flowcells, len(fcs)))
cand = []
for fc in fcs:
    for db in (by_fc[fc] if a.all else [rng.choice(by_fc[fc])]):
        sample = "/".join(db.split("/")[-3:-1])
        try:
            con = duckdb.connect(db, read_only=True)
            rows = con.execute("""select s.read_id, s.ani, s.lineage, i.fileid from sendsketch s
                join stats st on st.seqid = s.read_id join sequence_index i on i.seqid = s.read_id
                where s.status = 'C' and s.ani >= ? and st.length between ? and ?""", [a.min_ani, lo, hi]).fetchall()
            con.close()
        except Exception as e:
            print(f"skip {sample}: {e}"); continue
        cand += [(r, ani, lin, os.path.dirname(db) + "/fa/" + fid + ".fa", sample) for r, ani, lin, fid in rows]
        print(f"{sample}: {len(rows)} labelled reads", flush=True)

def ranks(lin):
    d = {}
    for f in lin.split(";"):
        k, _, v = f.partition("__")
        d[k] = v
    return d

by_genus = defaultdict(list)
for c in cand:
    r = ranks(c[2])
    if r.get("g"): by_genus[(r.get("d", ""), r["g"])].append(c)
keep = {}
for g in sorted(by_genus):
    v = by_genus[g]
    rng.shuffle(v)
    took, per_s = [], defaultdict(int)
    for c in v:
        if per_s[c[4]] < a.per_sample: took.append(c); per_s[c[4]] += 1
        if len(took) == a.per: break
    if len(took) == a.per:
        for c in took: keep[c[0]] = c
files = defaultdict(set)
for rid, c in keep.items(): files[c[3]].add(rid)
seqs = {}
for fa, want in files.items():
    cur = None
    for line in open(fa):
        if line[0] == ">":
            cur = line[1:].split()[0]
            cur = cur if cur in want else None
            if cur: seqs[cur] = []
        elif cur:
            seqs[cur].append(line.strip())
with open(a.out + ".fasta", "w") as o, open(a.out + ".tax.tsv", "w") as t:
    t.write("id\tdomain\tgroup\tphylum\tclass\torder\tfamily\tgenus\tspecies\tlineage_rank\tlineage\torganism\tani\tsample\n")
    for rid, (_, ani, lin, _, sample) in keep.items():
        if rid not in seqs: continue
        r = ranks(lin)
        o.write(f">{rid}\n{''.join(seqs[rid])}\n")
        t.write("\t".join([rid, r.get("d", ""), r.get("p", ""), r.get("p", ""), r.get("c", ""), r.get("o", ""),
                           r.get("f", ""), r.get("g", ""), r.get("s", ""), "", lin, "", f"{ani:.2f}", sample]) + "\n")
ng = len({tuple(ranks(c[2]).get(k) for k in 'dg') for c in keep.values()})
print(f"{len(seqs)} reads from {ng} genera ({len(cand)} labelled candidates, {len(fcs)} flowcells) -> {a.out}.fasta")
