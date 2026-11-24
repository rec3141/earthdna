#!/usr/bin/env python3
"""usage: gtdb_all_tax.py SHARDS.txt GENOMES.tsv OUT.tax.tsv

id, domain, group for every row of the shards, in list order (what conflicts.py reads)."""
import struct, sys
shards, genomes, out = sys.argv[1:4]
lines = open(genomes).read().splitlines()
h = lines[0].split("\t"); d0, g0 = h.index("domain"), h.index("group")
tax = {l.split("\t")[0]: l.split("\t") for l in lines[1:]}
with open(out, "w") as o:
    o.write("id\tdomain\tgroup\n")
    for p in open(shards).read().split():
        f = open(p, "rb"); vb = 2 if f.read(4) == b"FSB2" else 4; d, rows = struct.unpack("<iq", f.read(12))
        f.seek(16 + vb * d * rows)
        for fid in f.read().decode().split("\n"):
            r = tax["_".join(fid.split("_")[:2])]
            o.write(f"{fid}\t{r[d0]}\t{r[g0]}\n")
