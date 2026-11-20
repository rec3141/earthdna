#!/usr/bin/env python3
"""One genome's genes and sequence, for the globe viewer's fragment -> genes -> sequence views.

The viewer can open any drawn fragment of a packed genome: the point becomes an open loop of the fragment's genes,
and zooming on from there shows the sequence. Run once per genome to feature (a show needs a handful, not GTDB):

  genome_ACC.json     name, fragment length, contigs (id, length, offset into the sequence), genes
                      [contig, start, end, strand, type, name, product, attributes], and the globe points that are
                      fragments of this genome [point index, contig, start]
  genome_ACC.seq.txt  the contigs' sequence, upper case, concatenated
  genomes.json        the list of packed genomes (accession, name, number of points), updated in place

GLOBE.tax.tsv is the table pack_globe.py was given (the same file, so the point order matches); fragment ids are
<ACC>_<contig>_<start> as gtdb_prep.py writes them. The GFF is NCBI's or Prokka's/Bakta's: gene-level features
(CDS, rRNA, tRNA, ...) with Name/gene, locus_tag, product; any other attribute named with --attr is kept for colour
schemes (e.g. --attr cog --attr ec_number), as is every Dbxref "COG:..." id.

usage: pack_genome.py GLOBE.tax.tsv GLOBEDIR ACC genome.fna[.gz] genome.gff[.gz] [--name "Thermus aquaticus Y51MC23"]
                      [--frag 10000] [--attr KEY ...]
"""
import argparse, gzip, json, os, re, urllib.parse
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("tax"); ap.add_argument("outdir"); ap.add_argument("acc"); ap.add_argument("fasta"); ap.add_argument("gff")
ap.add_argument("--name", default=""); ap.add_argument("--frag", type=int, default=10000)
ap.add_argument("--attr", action="append", default=[])
a = ap.parse_args()
op = lambda p: gzip.open(p, "rt") if p.endswith(".gz") else open(p)

contigs, seqs = [], {}
with op(a.fasta) as f:
    cid = None
    for line in f:
        if line.startswith(">"): cid = line[1:].split()[0]; contigs.append(cid); seqs[cid] = []
        elif cid: seqs[cid].append(line.strip().upper())
seqs = {c: re.sub("[^ACGT]", "N", "".join(s)) for c, s in seqs.items()}
off, table = 0, []
for c in contigs: table.append(dict(id=c, len=len(seqs[c]), off=off)); off += len(seqs[c])
cix = {c: i for i, c in enumerate(contigs)}

KEEP = {"CDS", "rRNA", "tRNA", "tmRNA", "ncRNA", "pseudogene", "repeat_region", "mobile_genetic_element", "CRISPR"}
genes, seen = [], set()
with op(a.gff) as f:
    for line in f:
        if line.startswith("##FASTA"): break
        if line.startswith("#"): continue
        r = line.rstrip("\n").split("\t")
        if len(r) < 9 or r[2] not in KEEP or r[0] not in cix: continue
        at = {k: urllib.parse.unquote(v) for k, v in (x.split("=", 1) for x in r[8].split(";") if "=" in x)}
        key = (r[0], r[3], r[4], r[2])
        if key in seen: continue
        seen.add(key)
        extra = {k: at[k] for k in a.attr if k in at}
        cog = re.search(r"COG:?(COG\d+|[A-Z])\b", at.get("Dbxref", "") + "," + at.get("db_xref", ""))
        if cog and "cog" not in extra: extra["cog"] = cog.group(1)
        if at.get("locus_tag"): extra["locus_tag"] = at["locus_tag"]
        genes.append([cix[r[0]], int(r[3]) - 1, int(r[4]), -1 if r[6] == "-" else 1, r[2],
                      at.get("gene") or at.get("Name") or "", at.get("product", ""), extra])
genes.sort(key=lambda g: (g[0], g[1]))

# the globe's points that belong to this genome: pack_globe.py draws row perm[j] of the tax table as point j
ids = [l.split("\t", 1)[0] for l in list(open(a.tax))[1:]]
perm = np.random.default_rng(3).permutation(len(ids))
points = []
for j, r in enumerate(perm):
    i = ids[r]
    if not i.startswith(a.acc + "_"): continue
    c, s = i[len(a.acc) + 1:].rsplit("_", 1)
    if c in cix: points.append([j, cix[c], int(s)])
if not points: print(f"warning: no point of the globe is a fragment of {a.acc}; the viewer will have nothing to open")

os.makedirs(a.outdir, exist_ok=True)
name = a.name or a.acc
json.dump(dict(acc=a.acc, name=name, frag=a.frag, contigs=table, genes=genes, points=points),
          open(os.path.join(a.outdir, f"genome_{a.acc}.json"), "w"), separators=(",", ":"))
open(os.path.join(a.outdir, f"genome_{a.acc}.seq.txt"), "w").write("".join(seqs[c] for c in contigs))
lp = os.path.join(a.outdir, "genomes.json")
lst = [g for g in (json.load(open(lp)) if os.path.exists(lp) else []) if g["acc"] != a.acc]
lst.append(dict(acc=a.acc, name=name, points=len(points)))
json.dump(lst, open(lp, "w"), separators=(",", ":"))
print(f"{a.acc}: {len(contigs)} contigs, {off} bp, {len(genes)} genes, {len(points)} points on the globe")
