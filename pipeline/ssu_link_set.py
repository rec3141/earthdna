#!/usr/bin/env python3
"""Test set for linking a genome contig to its own 16S with Evo 2 (does the model know which 16S goes with which
genome?). Only genomes Evo 2 cannot have seen: GTDB r226 bacterial species representatives absent from r214 and r220
(Evo 2 trained on r214.1 representatives plus genomes new in r220).

For each of N test genomes (at most one per family, chosen by CRC32 of the accession):
  contig   the central --contig-len bp of its longest contig that does not carry the 16S (at least --min-contig)
  true     its own full-length 16S (GTDB ssu_all, >= --min-ssu bp)
  decoys   --n-family 16S from other genera of the same family, --n-order from other families of the same order
           (any r226 representative with a full-length 16S; lowest CRC32 first)

  OUT.contigs.fa   >acc
  OUT.ssu.fa       >acc (one per genome used: tests and decoys)
  OUT.pairs.tsv    test_acc  candidate_acc  relation (self / family / order)  candidate taxonomy

usage: ssu_link_set.py GTDB_RELEASE_DIR OUT [--n 50] [--contig-len 20000]
  GTDB_RELEASE_DIR: release226 with metadata/, skani/genome_paths.tsv, ssu/ssu_all_r226.fna.gz; ../release214 and
  ../release220 with {bac120,ar53}_taxonomy_r*.tsv.gz
"""
import argparse, gzip, os, zlib, collections

ap = argparse.ArgumentParser()
ap.add_argument("gtdb"); ap.add_argument("out")
ap.add_argument("--n", type=int, default=50); ap.add_argument("--contig-len", type=int, default=20000)
ap.add_argument("--min-contig", type=int, default=30000); ap.add_argument("--min-ssu", type=int, default=1400)
ap.add_argument("--n-family", type=int, default=10); ap.add_argument("--n-order", type=int, default=10)
a = ap.parse_args()
crc = lambda s: zlib.crc32(s.encode())
strip = lambda acc: acc.split("_", 1)[1] if acc[:3] in ("RS_", "GB_") else acc

seen = set()
for rel in ("214", "220"):
    for dom in ("bac120", "ar53"):
        with gzip.open(os.path.join(a.gtdb, "..", f"release{rel}", f"{dom}_taxonomy_r{rel}.tsv.gz"), "rt") as f:
            seen.update(strip(l.split("\t")[0]) for l in f)
# also drop any genome whose version-less accession was seen (GCA/GCF twins and re-versioned assemblies)
seen_base = {s.split(".")[0][4:] for s in seen}

tax = {}
with gzip.open(os.path.join(a.gtdb, "metadata", "bac120_metadata_r226.tsv.gz"), "rt") as f:
    h = f.readline().rstrip("\n").split("\t"); c = {k: i for i, k in enumerate(h)}
    for l in f:
        r = l.rstrip("\n").split("\t")
        if r[c["gtdb_representative"]] == "t": tax[strip(r[0])] = r[c["gtdb_taxonomy"]]

ssu = {}  # genome -> (contig, seq) of its longest reported 16S
with gzip.open(os.path.join(a.gtdb, "ssu", "ssu_all_r226.fna.gz"), "rt") as f:
    name = None; buf = []
    def put():
        if name is None: return
        g, contig = name.split("~", 1); g = strip(g); s = "".join(buf).upper()
        if g in tax and len(s) >= a.min_ssu and (g not in ssu or len(s) > len(ssu[g][1])): ssu[g] = (contig, s)
    for l in f:
        if l[0] == ">": put(); name = l[1:].split()[0]; buf = []
        else: buf.append(l.strip())
    put()
print(f"{len(tax):,} bacterial representatives, {len(ssu):,} with a 16S >= {a.min_ssu} bp")

rank = lambda t, k: ";".join(t.split(";")[:k])          # d;p;c;o;f;g;s: k=4 to order, 5 family, 6 genus
new = sorted((g for g in ssu if g not in seen and g.split(".")[0][4:] not in seen_base), key=crc)
print(f"{len(new):,} of them are new in r226 (not in r214 or r220)")
by_fam, by_ord = collections.defaultdict(list), collections.defaultdict(list)
for g in sorted(ssu, key=crc):
    by_fam[rank(tax[g], 5)].append(g); by_ord[rank(tax[g], 4)].append(g)

paths = {}   # accession -> genome FASTA (genome_paths.tsv: file name, directory under skani/)
for l in open(os.path.join(a.gtdb, "skani", "genome_paths.tsv")):
    p = l.split()
    if len(p) >= 2: paths[p[0].split("_genomic")[0]] = os.path.join(a.gtdb, "skani", p[1], p[0])

def contigs(path):
    with gzip.open(path, "rt") if path.endswith(".gz") else open(path) as f:
        name, buf = None, []
        for l in f:
            if l[0] == ">":
                if name: yield name, "".join(buf)
                name, buf = l[1:].split()[0], []
            else: buf.append(l.strip())
        if name: yield name, "".join(buf)

picked, fams = [], set()
cf, sf, pairs = open(a.out + ".contigs.fa", "w"), open(a.out + ".ssu.fa", "w"), open(a.out + ".pairs.tsv", "w")
pairs.write("test\tcandidate\trelation\ttaxonomy\n")
ssu_written = set()
def write_ssu(g):
    if g not in ssu_written: sf.write(f">{g}\n{ssu[g][1]}\n"); ssu_written.add(g)
for g in new:
    if len(picked) >= a.n: break
    fam, odr, gen = rank(tax[g], 5), rank(tax[g], 4), rank(tax[g], 6)
    if fam in fams: continue
    fdec = [x for x in by_fam[fam] if rank(tax[x], 6) != gen][:a.n_family]
    odec = [x for x in by_ord[odr] if rank(tax[x], 5) != fam][:a.n_order]
    if len(fdec) < a.n_family or len(odec) < a.n_order: continue
    path = paths.get(g)
    if not path or not os.path.exists(path): continue
    best = max(((n, s) for n, s in contigs(path) if n != ssu[g][0] and len(s) >= a.min_contig),
               key=lambda x: len(x[1]), default=None)
    if not best: continue
    s = best[1]; m = (len(s) - a.contig_len) // 2
    cf.write(f">{g}\n{s[m:m + a.contig_len].upper()}\n")
    for x, rel in [(g, "self")] + [(x, "family") for x in fdec] + [(x, "order") for x in odec]:
        write_ssu(x); pairs.write(f"{g}\t{x}\t{rel}\t{tax[x]}\n")
    picked.append(g); fams.add(fam)
print(f"{len(picked)} test genomes (one per family), {len(ssu_written)} 16S sequences -> {a.out}.*")
