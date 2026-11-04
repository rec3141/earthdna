#!/usr/bin/env python3
"""Cut an amplicon region out of every aligned reference that spans it, for region-matched ASV placement.

The region's alignment columns come from references where both primers match (in-silico PCR on the ungapped
sequence, IUPAC codes allowed); the commonest start and end columns then cut every sequence that covers them,
including the many whose primer sites have mismatches or lie outside a trimmed read. Primers are excluded.

usage: region_refs.py SILVA_..._full_align_trunc.fasta.gz FWD REV out.fasta [--sample 200000]
  e.g. V4 515F-Y / 806R-B:  GTGYCAGCMGCCGCGGTAA  GGACTACNVGGGTWTCTAAT
"""
import argparse, collections, gzip, re
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("aln"); ap.add_argument("fwd"); ap.add_argument("rev"); ap.add_argument("out")
ap.add_argument("--sample", type=int, default=200000, help="sequences used to find the region's columns")
ap.add_argument("--min-len", type=int, default=100)
a = ap.parse_args()

IUPAC = dict(A="A", C="C", G="G", T="T", U="T", R="[AG]", Y="[CT]", S="[CG]", W="[AT]", K="[GT]", M="[AC]",
             B="[CGT]", D="[AGT]", H="[ACT]", V="[ACG]", N="[ACGT]")
rc = lambda p: p[::-1].translate(str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN"))
fwd = re.compile("".join(IUPAC[c] for c in a.fwd.upper()))
rev = re.compile("".join(IUPAC[c] for c in rc(a.rev.upper())))

def records():
    name, seq = None, []
    with gzip.open(a.aln, "rt") as f:
        for line in f:
            if line[0] == ">":
                if name: yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
    if name: yield name, "".join(seq)

def bases(aln):
    arr = np.frombuffer(aln.encode(), np.uint8)
    cols = np.where(((arr >= 65) & (arr <= 90)) | ((arr >= 97) & (arr <= 122)))[0]
    return cols, aln.upper().replace("U", "T")

# pass 1: alignment columns of the first base after the forward primer and the last base before the reverse primer
starts, ends = collections.Counter(), collections.Counter()
for n, (name, aln) in enumerate(records()):
    if n >= a.sample: break
    cols, up = bases(aln)
    s = "".join(up[c] for c in cols)
    m1 = fwd.search(s)
    if not m1: continue
    m2 = rev.search(s, m1.end())
    if not m2 or m2.start() <= m1.end(): continue
    starts[cols[m1.end()]] += 1; ends[cols[m2.start() - 1]] += 1
if not starts: raise SystemExit("no sequence in the sample matched both primers")
c0, c1 = starts.most_common(1)[0][0], ends.most_common(1)[0][0]
print(f"region: alignment columns {c0}-{c1} ({starts[c0]} / {ends[c1]} primer-matched sequences agree)")

# pass 2: every sequence with bases on both sides of the region
kept = 0; lens = []
with open(a.out, "w") as out:
    for name, aln in records():
        cols, up = bases(aln)
        if len(cols) == 0 or cols[0] > c0 or cols[-1] < c1: continue
        reg = "".join(ch for ch in up[c0:c1 + 1] if ch.isalpha())
        if len(reg) < a.min_len: continue
        out.write(f">{name}\n{reg}\n"); kept += 1; lens.append(len(reg))
print(f"{kept} references span the region, median length {int(np.median(lens))}")
