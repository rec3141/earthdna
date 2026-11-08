#!/usr/bin/env python3
"""The shared fragment set for the DNA-model globes (microscape.app/earthdna/models): one window per GTDB genus, so
every model (5-mers, Evo 2 7B / 20B, DNABERT) embeds the same sequences and the viewer can toggle between their maps.

For each genus: the species representative with the lowest CRC32(accession) among those with a fragment in the
globe's reference FASTA (gtdb_prep.py's .ref.fasta.gz, so every window also sits on the GTDB globe), then that
genome's lowest-CRC32 fragment, cut to its central --len bp (windows with more than 1% N are skipped).

  OUT.fa        windows, ids "<fragment id>"   (fragment id = <accession>_<contig>_<offset>)
  OUT.tax.tsv   id + the genome's taxonomy columns (pack_globe.py / fastsne anchor input)

usage: models_set.py GENOMES.tsv REF.fasta.gz OUT [--len 4000]
"""
import argparse, gzip, zlib
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("genomes"); ap.add_argument("fasta"); ap.add_argument("out"); ap.add_argument("--len", type=int, default=4000)
a = ap.parse_args()
crc = lambda s: zlib.crc32(s.encode())

g = pd.read_csv(a.genomes, sep="\t")
acc_of = lambda fid: fid.rsplit("_", 2)[0]
# genome -> its lowest-hash fragment id (one pass over the FASTA ids)
best = {}
with gzip.open(a.fasta, "rt") as f:
    for l in f:
        if l[0] == ">":
            fid = l[1:].split()[0]; acc = acc_of(fid); h = crc(fid)
            if acc not in best or h < best[acc][0]: best[acc] = (h, fid)
g = g[g.accession.isin(best)].copy()
g["h"] = g.accession.map(crc)
pick = g.sort_values("h").groupby("genus", sort=False).head(1)
want = {best[acc][1]: acc for acc in pick.accession}
print(f"{g.genus.nunique():,} genera with a reference fragment; {len(want):,} windows to cut")

kept, seq, fid = [], [], None
def flush(out):
    if fid in want:
        s = "".join(seq); m = (len(s) - a.len) // 2; w = s[max(m, 0):max(m, 0) + a.len]
        if len(w) == a.len and w.upper().count("N") <= 0.01 * a.len:
            out.write(f">{fid}\n{w}\n"); kept.append(fid)
with gzip.open(a.fasta, "rt") as f, open(f"{a.out}.fa", "w") as out:
    for l in f:
        if l[0] == ">":
            flush(out); fid = l[1:].split()[0]; seq = []
        elif fid in want: seq.append(l.strip())
    flush(out)
tax = pick.set_index("accession").drop(columns=["contigs", "bp", "contigs_ge_frag", "fragments", "ref_fragments", "h"])
t = tax.loc[[want[k] for k in kept]].reset_index(drop=True); t.insert(0, "id", kept)
t.to_csv(f"{a.out}.tax.tsv", sep="\t", index=False)
print(f"{len(kept):,} windows of {a.len} bp -> {a.out}.fa, {a.out}.tax.tsv")
