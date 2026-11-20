#!/usr/bin/env python3
"""mOTUs 4.1 marker-gene database -> EarthDNA inputs (the files silva_prep.py writes, minus the k-mer table).

  OUT.fasta     the chosen marker-gene sequences, id = mOTUs header (mOTU.MGC_<COG>_<cluster>.MG_<n>)
  OUT.tax.tsv   id, domain, group, phylum, class, order, family, genus, species, lineage_rank, lineage, organism
                (silva_prep.py's columns; rows in OUT.fasta order; organism = the mOTU id, the key for --traits)
  OUT.anchor    group ids for fastsne --anchor-file: domain, group, class, order (-1 = rank missing), nested as in
                silva_prep.py so a name repeated under two parents gets two ids
  OUT.source.tsv id, mOTU, COG, length, source genome (map.tsv ORIG), whether it is the mOTU's representative genome
  OUT.stats.tsv sequences, mOTUs and length quartiles per COG among the chosen sequences

The k-mer table comes from EarthDNA's own pipeline/kmers.py run on OUT.fasta (forward strand: marker genes are
stored in gene orientation), with no --min-len, so its rows line up with OUT.tax.tsv; length filtering happens here.

Taxonomy is the mOTU's GTDB R226 lineage by 80% majority vote over its genomes
(mOTUsv4.1.gtdb.taxonomy.80mv.tsv.gz). For prokaryotes `group` is the phylum, as in silva_prep.py. Placeholder
species ("Unknown Prevotella mOTUv4.0_000001") are left blank; GTDB species such as "UBA11524 sp000437595" are kept,
being real GTDB species clusters.

Sequences: the "mOTUv4.0_unassigned" bin (marker genes in no species-level mOTU) and the database's blocklist are
skipped. --cog picks one of the ten marker-gene families, or "all". --per-motu K keeps K sequences per mOTU and COG
(0 = all): the one from the mOTU's representative genome first (mOTUsv4.1.gtdb.taxonomy.rep.tsv.gz GENOME, matched
through map.tsv ORIG), then the rest in file order. One COG at one sequence per mOTU gives one point per species,
the counterpart of SILVA NR99's one point per SSU cluster.

With --concat (and --cog all) each mOTU becomes one record, its markers of all ten COGs joined by N: every mOTU gets
a point, and the profile pools ~15 kb of marker sequence instead of one ~1.1 kb gene. The stats then count mOTUs
by the set of COGs they carry.

usage: python3 -I motus_prep.py DB_DIR OUT [--cog COG0012|all] [--per-motu 1] [--min-len 300] [--max-len 0]
                                [--concat]
  DB_DIR holds the extracted db_mOTU/ files: mOTUsv4.1.db.fna.gz, mOTUsv4.1.gtdb.taxonomy.80mv.tsv.gz,
  mOTUsv4.1.gtdb.taxonomy.rep.tsv.gz, mOTUsv4.1.map.tsv.gz, mOTUsv4.1.db.blocklist.gz
"""
import argparse, collections, gzip, os, re

ap = argparse.ArgumentParser()
ap.add_argument("db"); ap.add_argument("out")
ap.add_argument("--cog", default="COG0012")
ap.add_argument("--per-motu", type=int, default=1)
ap.add_argument("--min-len", type=int, default=300)
ap.add_argument("--max-len", type=int, default=0, help="0 = no upper limit")
ap.add_argument("--concat", action="store_true", help="one record per mOTU: its chosen sequences of every COG in COG "
                "order, joined by N (k-mers spanning a junction are skipped by kmers.py); id = the mOTU")
a = ap.parse_args()
if a.concat and not a.per_motu: ap.error("--concat needs --per-motu >= 1")
P = lambda name: os.path.join(a.db, "mOTUsv4.1." + name)
want = lambda cog: a.cog == "all" or cog == a.cog

def table(path):
    with gzip.open(path, "rt") as f:
        head = next(f).lstrip("#").rstrip("\n").split("\t")
        for line in f:
            yield dict(zip(head, line.rstrip("\n").split("\t")))

tax = {}
for r in table(P("gtdb.taxonomy.80mv.tsv.gz")):
    d = {}
    for part in r["GTDB"].split(";"):
        if len(part) > 3 and part[1:3] == "__": d[part[0]] = part[3:].strip()
    tax[r["MOTU"]] = d
rep_genome = {r["MOTU"]: r["GENOME"] for r in table(P("gtdb.taxonomy.rep.tsv.gz"))}
with gzip.open(P("db.blocklist.gz"), "rt") as f: blocked = {l.split()[0] for l in f if l.strip()}
origin, from_rep = {}, set()
for r in table(P("map.tsv.gz")):
    if not want(r["COG"]): continue
    g = r["ORIG"].rsplit(".", 1)[0]
    origin[r["MG"]] = g
    if rep_genome.get(r["MOTU"]) == g: from_rep.add(r["MG"])

def row(motu):
    d = tax.get(motu, {})
    dom, ph, cl, od, fa, ge = (d.get(k, "") for k in "dpcofg")
    sp = d.get("s", "")
    if not sp or sp.startswith("Unknown") or motu in sp: sp = ""
    lin_rank, lin = "", ""
    for rk, v in (("family", fa), ("order", od), ("class", cl), ("phylum", ph)):
        if v: lin_rank, lin = rk, v; break
    return [dom, ph, ph, cl, od, fa, ge, sp, lin_rank, lin, motu]

HDR = re.compile(r"^(?P<motu>\S+?)\.MGC_(?P<cog>COG\d+)_\d+\.MG_\d+$")
n = collections.Counter()
chosen = collections.defaultdict(list)          # (motu, cog) -> [(rank, order, sid, seq)], when --per-motu > 0
order = []                                      # (motu, cog) keys in first-seen order
fo = open(a.out + ".fasta", "w")
kept = []                                       # (sid, motu, cog, length)

def write(sid, s, motu, cog):
    kept.append((sid, motu, cog, len(s)))
    fo.write(">" + sid + "\n")
    for i in range(0, len(s), 80): fo.write(s[i:i + 80] + "\n")

def consider(sid, seq):
    n["read"] += 1
    m = HDR.match(sid)
    if not m: n["unparsed header"] += 1; return
    motu, cog = m["motu"], m["cog"]
    if not want(cog): return
    if motu.endswith("_unassigned"): n["unassigned bin"] += 1; return
    if sid in blocked: n["blocklisted"] += 1; return
    s = "".join(seq)
    if len(s) < a.min_len: n["short"] += 1; return
    if a.max_len and len(s) > a.max_len: n["long"] += 1; return
    if not a.per_motu: write(sid, s, motu, cog); return
    key = (motu, cog)
    lst = chosen[key]
    if not lst: order.append(key)
    lst.append((0 if sid in from_rep else 1, len(lst), sid, s))
    if len(lst) > a.per_motu:
        lst.sort(); lst.pop(); n["over per-mOTU cap"] += 1

sid, seq = None, []
with gzip.open(P("db.fna.gz"), "rt") as f:
    for line in f:
        if line.startswith(">"):
            if sid is not None: consider(sid, seq)
            sid, seq = line[1:].split()[0], []
        else:
            seq.append(line.strip())
if sid is not None: consider(sid, seq)
if a.concat:
    by_motu = collections.defaultdict(list)
    for motu, cog in order: by_motu[motu].append(cog)
    for motu, cogs in by_motu.items():
        parts = [(s_id, s) for c in sorted(cogs) for _, _, s_id, s in sorted(chosen[motu, c])]
        joined = "N".join(s for _, s in parts)
        kept.append((motu, motu, "+".join(sorted(cogs)), len(joined)))
        fo.write(">" + motu + "\n")
        for i in range(0, len(joined), 80): fo.write(joined[i:i + 80] + "\n")
        n["joined sequences"] += len(parts)
        n["joined from representative"] += sum(s_id in from_rep for s_id, _ in parts)
else:
    for key in order:
        for _, _, s_id, s in sorted(chosen[key]): write(s_id, s, *key)
fo.close()

rows = [row(m) for _, m, _, _ in kept]
cols = ["domain", "group", "phylum", "class", "order", "family", "genus", "species", "lineage_rank", "lineage", "organism"]
with open(a.out + ".tax.tsv", "w") as f:
    f.write("id\t" + "\t".join(cols) + "\n")
    for (i, *_), r in zip(kept, rows): f.write(i + "\t" + "\t".join(x.replace("\t", " ") for x in r) + "\n")
codes = []
for j, deep in ((0, 1), (1, 2), (3, 4), (4, 5)):  # domain, group, class, order
    ix, col = {}, []
    for r in rows:
        key = ";".join(r[:deep]) if r[j] else None
        col.append(ix.setdefault(key, len(ix)) if key else -1)
    codes.append(col)
with open(a.out + ".anchor", "w") as f:
    for r in zip(*codes): f.write(" ".join(map(str, r)) + "\n")
with open(a.out + ".source.tsv", "w") as f:
    f.write("id\tmotu\tcog\tlength\tsource_genome\tfrom_representative\n")
    for sid_, m, c, L in kept: f.write(f"{sid_}\t{m}\t{c}\t{L}\t{origin.get(sid_, '')}\t{int(sid_ in from_rep)}\n")

lens, motus_per = collections.defaultdict(list), collections.defaultdict(set)
for _, m, c, L in kept: lens[c].append(L); motus_per[c].add(m)
def q(v):
    v = sorted(v); k = len(v)
    return [v[0], v[k // 4], v[k // 2], v[3 * k // 4], v[-1]]
with open(a.out + ".stats.tsv", "w") as f:
    f.write("cog\tsequences\tmotus\tlen_min\tlen_q1\tlen_median\tlen_q3\tlen_max\n")
    for c in sorted(lens, key=lambda c: (c.count("+"), c)): f.write(f"{c}\t{len(lens[c])}\t{len(motus_per[c])}\t" + "\t".join(map(str, q(lens[c]))) + "\n")

motus_kept = {m for _, m, _, _ in kept}
print(f"{n['read']} sequences read; kept {len(kept)} ({a.cog}, per mOTU {a.per_motu or 'all'}, min length {a.min_len}"
      f"{', max ' + str(a.max_len) if a.max_len else ''}); skipped " + ", ".join(f"{v} {k}" for k, v in n.items() if k != "read"))
n_rep = n["joined from representative"] if a.concat else sum(1 for s, *_ in kept if s in from_rep)
print(f"{len(motus_kept)} mOTUs, {n_rep} sequences from the representative genome, "
      f"{len(motus_kept - set(tax))} mOTUs without a taxonomy row; {len(set(r[0] for r in rows if r[0]))} domains, "
      f"{len(set(r[2] for r in rows if r[2]))} phyla, {len(set(r[6] for r in rows if r[6]))} genera, "
      f"{len(set(r[7] for r in rows if r[7]))} named species")
print(f"-> {a.out}.fasta / .tax.tsv / .anchor / .source.tsv / .stats.tsv")
