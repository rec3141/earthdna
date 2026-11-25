#!/usr/bin/env python3
"""usage: gtdb_all_tax.py SHARDS.txt GENOMES.tsv OUT.tax.tsv [--full]

id, domain, group for every row of the shards, in list order (what conflicts.py reads). With --full, every taxonomy
column of GENOMES.tsv from domain on (gtdb_reselect.py's REF.tax.tsv layout), for gtdb_globe_set.py when every
fragment is in the embedding."""
import struct, sys
full = "--full" in sys.argv[1:]
shards, genomes, out = [a for a in sys.argv[1:] if a != "--full"][:3]
lines = open(genomes).read().splitlines()
h = lines[0].split("\t"); d0, g0 = h.index("domain"), h.index("group")
tax = {l.split("\t")[0]: l.split("\t") for l in lines[1:]}
with open(out, "w") as o:
    o.write("id\t" + "\t".join(h[d0:] if full else ["domain", "group"]) + "\n")
    for p in open(shards).read().split():
        f = open(p, "rb"); vb = 2 if f.read(4) == b"FSB2" else 4; d, rows = struct.unpack("<iq", f.read(12))
        f.seek(16 + vb * d * rows)
        for fid in f.read().decode().split("\n"):
            r = tax["_".join(fid.split("_")[:2])]
            o.write(fid + "\t" + "\t".join(r[d0:] if full else [r[d0], r[g0]]) + "\n")
