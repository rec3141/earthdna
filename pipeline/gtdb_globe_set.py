#!/usr/bin/env python3
"""Reference + a projected subsample -> one point set for pack_globe.py.

The globe viewer holds every point in memory as base64 text and draws them in one call, so the full projected set
(tens of millions of fragments) stays on disk and the globe gets the embedded reference plus EXTRA projected
non-reference fragments, chosen by a hash of the id (the same fragments in every layout). Projected rows take their
taxonomy from the genome table (accession = the first two '_' fields of the id) and conflict flag 0 (conflicts are
computed on the reference kNN graph only).

  OUT.tax.tsv, OUT.conflicts.tsv, OUT.NAME.tsv (id lon lat) per layout

With a reference bigger than the globe can draw, --ref-shown F keeps only the reference fragments whose id hash
(CRC32, as gtdb_prep.py chooses them) is below F * 2^32, e.g. 0.065 shows the 6.5% subset of a 30% reference.

usage: gtdb_globe_set.py REF.tax.tsv GENOMES.tsv REF.conflicts.tsv OUT --extra 1500000 [--ref-shown 0.065] \
           --layout similarity=ref.sphere.tsv:all.sphere.tsv [--layout sculpt=ref.sculpt.tsv:all.sculpt.tsv]
"""
import argparse, zlib

ap = argparse.ArgumentParser()
ap.add_argument("reftax"); ap.add_argument("genomes"); ap.add_argument("conflicts"); ap.add_argument("out")
ap.add_argument("--extra", type=int, default=1500000); ap.add_argument("--layout", action="append", required=True)
ap.add_argument("--ref-shown", type=float, default=1.0, help="show only reference ids with CRC32 below this fraction")
a = ap.parse_args()

gen, total = {}, 0
with open(a.genomes) as f:
    h = next(f).rstrip("\n").split("\t")
    t0 = h.index("domain")
    for line in f:
        r = line.rstrip("\n").split("\t")
        gen[r[0]] = r[t0:]; total += int(r[h.index("fragments")])
refids = set()
with open(a.reftax) as f, open(a.out + ".tax.tsv", "w") as o:
    o.write(next(f))
    RT = int(a.ref_shown * 2 ** 32)
    shown = lambda i: a.ref_shown >= 1 or zlib.crc32(i.encode()) < RT
    for line in f:
        i = line.split("\t", 1)[0]
        if shown(i): refids.add(i); o.write(line)
T = int(min(1.0, a.extra / max(1, total - len(refids))) * 2 ** 32)
pick = lambda i: i not in refids and zlib.crc32(b"globe:" + i.encode()) < T

extra = []
for k, spec in enumerate(a.layout):
    name, paths = spec.split("=", 1)
    refp, allp = paths.split(":")
    with open(a.out + f".{name}.tsv", "w") as o:
        o.write("id\tlon\tlat\n")
        with open(refp) as f:
            next(f)
            for line in f:
                r = line.rstrip("\n").split("\t")[:3]
                if r[0] in refids: o.write("\t".join(r) + "\n")
        n = 0
        with open(allp) as f:
            next(f)
            for line in f:
                i, lo, la = line.split("\t", 3)[:3]
                if pick(i):
                    o.write(f"{i}\t{lo}\t{la.strip()}\n"); n += 1
                    if k == 0: extra.append(i)
    print(f"{name}: {len(refids)} reference + {n} projected")
with open(a.out + ".tax.tsv", "a") as o:
    for i in extra:
        acc = "_".join(i.split("_", 2)[:2])
        o.write(i + "\t" + "\t".join(gen[acc]) + "\n")
with open(a.out + ".conflicts.tsv", "w") as o, open(a.conflicts) as f:
    o.write(next(f))
    for line in f:
        if line.split("\t", 1)[0] in refids: o.write(line)
    for i in extra: o.write(f"{i}\t0\t\t\n")
print(f"{len(refids) + len(extra)} points -> {a.out}.tax.tsv / .conflicts.tsv")
