#!/usr/bin/env python3
"""usage: gtdb_all_anchor.py SHARDS.txt GENOMES.tsv OUT.anchor

Anchor groups (domain, group, class, order; gtdb_reselect.py's codes) for every row of the shards, in list order."""
import struct, sys
shards, genomes, out = sys.argv[1:4]
lines = open(genomes).read().splitlines()
c0 = lines[0].split("\t").index("domain")
tax = {l.split("\t")[0]: l.split("\t")[c0:] for l in lines[1:]}
ix = [{} for _ in range(4)]
n = 0
with open(out, "w") as o:
    for p in open(shards).read().split():
        f = open(p, "rb"); vb = 2 if f.read(4) == b"FSB2" else 4; d, rows = struct.unpack("<iq", f.read(12))
        f.seek(16 + vb * d * rows)
        for fid in f.read().decode().split("\n"):
            r = tax["_".join(fid.split("_")[:2])]
            c = []
            for k, (j, deep) in enumerate(((0, 1), (1, 2), (3, 4), (4, 5))):
                c.append(ix[k].setdefault(";".join(r[:deep]), len(ix[k])) if r[j] else -1)
            o.write("%d %d %d %d\n" % tuple(c)); n += 1
        print(p, rows, n, flush=True)
