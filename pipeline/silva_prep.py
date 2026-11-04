#!/usr/bin/env python3
"""SILVA SSU fasta -> inputs for the globe.

  OUT.fsb       k-mer counts (float32, fastsne binary input; rows in fasta order)
  OUT.tax.tsv   id, domain, group, phylum, class, order, family, genus, species, lineage_rank, lineage, organism
  OUT.anchor    group ids for fastsne --anchor-file: domain, group, class, order (-1 = rank missing)

Ranks come from SILVA's rank file (tax_slv_ssu_*.txt), not from positions in the path: eukaryote paths stop at
whatever node SILVA classified to (only ~13% reach genus). The last header field is the source organism, which for
bacteria and archaea is often the host or sample (a Pseudomonas 16S named "Rhizoctonia solani"), so a binomial is
kept as the species only when its genus matches SILVA's genus; eukaryotes without a SILVA genus take genus and
species from the binomial.

usage: silva_prep.py SILVA_..._trunc.fasta.gz tax_slv_ssu_NNN.txt.gz OUT [--k 5]
"""
import argparse, gzip, re, struct
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("fasta"); ap.add_argument("ranks"); ap.add_argument("out")
ap.add_argument("--k", type=int, default=5)
a = ap.parse_args()

rank = {}
with gzip.open(a.ranks, "rt") as f:
    for line in f:
        p, _, r, *_ = line.rstrip("\n").split("\t")
        rank[p] = r

BAD = re.compile(r"^(uncultured|unidentified|unclassified|metagenome|environmental|bacterium|archaeon|eukaryote|"
                 r"marine|soil|freshwater|human|mouse|gut)$", re.I)
def binomial(org):
    w = org.replace("[", "").replace("]", "").split()
    if w and w[0] == "Candidatus": w = w[1:]
    if len(w) < 2 or not re.fullmatch(r"[A-Z][a-z]+", w[0]) or not re.fullmatch(r"[a-z][a-z-]+", w[1]): return None
    if BAD.match(w[0]) or w[1] in ("sp", "bacterium", "archaeon", "endosymbiont", "symbiont", "clone", "isolate"): return None
    return w[0] + " " + w[1]

def taxonomy(header_tax):
    t = header_tax.split(";")
    path, org = t[:-1], t[-1]
    R = {}
    for d in range(1, len(path) + 1):
        R.setdefault(rank.get(";".join(path[:d]) + ";"), path[d - 1])
    dom = path[0] if path else ""
    group = (path[1] if len(path) > 1 else "") if dom == "Eukaryota" else R.get("phylum", "")
    clean = lambda x: "" if not x or x == "Incertae Sedis" else x
    g, sp = clean(R.get("genus")), binomial(org)
    if dom != "Eukaryota":
        gw = [x.split()[0] for x in g.replace("Candidatus ", "").split("-") if x.split()]
        if sp and sp.split()[0] not in gw: sp = None
    elif not g and sp: g = sp.split()[0]
    elif g and sp and sp.split()[0] != g.split()[0]: sp = None
    lin_rank, lin = "", ""
    for rk in ("family", "order", "class", "phylum"):
        if clean(R.get(rk)): lin_rank, lin = rk, R[rk]; break
    return [dom, group, clean(R.get("phylum")), clean(R.get("class")), clean(R.get("order")),
            clean(R.get("family")), g, sp or "", lin_rank, lin, org]

lut = np.full(256, 255, np.uint8)
for i, c in enumerate(b"ACGU"): lut[c] = i; lut[c + 32] = i
lut[ord("T")] = lut[ord("t")] = 3
k = a.k
pw = (4 ** np.arange(k - 1, -1, -1)).astype(np.int64)
def counts(s):
    c = lut[np.frombuffer(s, np.uint8)]
    bad = c > 3
    c = np.where(bad, 0, c).astype(np.int64)
    if len(c) < k: return np.zeros(4 ** k, np.float32)
    w = np.lib.stride_tricks.sliding_window_view(c, k) @ pw
    ok = ~np.lib.stride_tricks.sliding_window_view(bad, k).any(1)  # skip k-mers spanning ambiguity codes
    return np.bincount(w[ok], minlength=4 ** k).astype(np.float32)

ids, rows = [], []
fsb = open(a.out + ".fsb", "wb")
fsb.write(b"FSB1" + struct.pack("<iq", 4 ** k, 0))  # N patched at the end
def flush(seq):
    fsb.write(counts(b"".join(seq)).tobytes())
seq = None
with gzip.open(a.fasta, "rb") as f:
    for line in f:
        if line[:1] == b">":
            if seq is not None: flush(seq)
            sid, _, tx = line[1:].decode().rstrip("\n").partition(" ")
            ids.append(sid); rows.append(taxonomy(tx)); seq = []
        else:
            seq.append(line.strip())
if seq is not None: flush(seq)
fsb.write("\n".join(ids).encode())
fsb.seek(8); fsb.write(struct.pack("<q", len(ids))); fsb.close()

cols = ["domain", "group", "phylum", "class", "order", "family", "genus", "species", "lineage_rank", "lineage", "organism"]
with open(a.out + ".tax.tsv", "w") as f:
    f.write("id\t" + "\t".join(cols) + "\n")
    for i, r in zip(ids, rows): f.write(i + "\t" + "\t".join(x.replace("\t", " ") for x in r) + "\n")
codes = []
for j, deep in ((0, 1), (1, 2), (3, 4), (4, 5)):  # domain, group, class, order; nested so names repeat safely
    ix, col = {}, []
    for r in rows:
        key = ";".join(r[:deep]) if r[j] else None
        col.append(ix.setdefault(key, len(ix)) if key else -1)
    codes.append(col)
with open(a.out + ".anchor", "w") as f:
    for row in zip(*codes): f.write(" ".join(map(str, row)) + "\n")
print(f"{len(ids)} sequences, {4 ** k} {k}-mer columns -> {a.out}.fsb / .tax.tsv / .anchor")
