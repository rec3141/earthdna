#!/usr/bin/env python3
"""A new reference set for a GTDB globe, taken from gtdb_prep.py's all-fragment shards (no re-reading of genomes).

Fragments are chosen as gtdb_prep.py chooses them, by CRC32 of the fragment id below --ref * 2^32, so a bigger
fraction contains every smaller one (and the earlier OUT.ref.fasta.gz still holds sequences for a subset). Every
genome then gets at least --min-per-genome reference fragments (all of them if it has fewer): the lowest-hash
fragments it has, so small genomes that a hash fraction misses still shape the map.

Writes OUT.ref.fsb, OUT.ref.tax.tsv and OUT.ref.anchor in gtdb_prep.py's formats (taxonomy from OUT_OLD.genomes.tsv).

usage: gtdb_reselect.py ALL_FSB_LIST GENOMES_TSV OUT [--ref 0.30] [--min-per-genome 5]
"""
import argparse, struct, zlib
from collections import defaultdict
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("shards"); ap.add_argument("genomes"); ap.add_argument("out")
ap.add_argument("--ref", type=float, default=0.30); ap.add_argument("--min-per-genome", type=int, default=5)
a = ap.parse_args()
T = int(a.ref * 2 ** 32)

lines = open(a.genomes).read().splitlines()
gcols = lines[0].split("\t")
c0 = gcols.index("domain")
tax = {l.split("\t")[0]: l.split("\t")[c0:] for l in lines[1:]}
cols = gcols[c0:]
acc_of = lambda fid: "_".join(fid.split("_")[:2])          # GCA_000008085.1_<contig>_<start>

paths = [p for p in open(a.shards).read().split() if p]
def header(p):
    f = open(p, "rb"); f.read(4); d, n = struct.unpack("<iq", f.read(12)); return f, d, n
def ids_of(p):
    f, d, n = header(p); f.seek(16 + 4 * d * n); return f.read().decode().split("\n")

# pass 1: hashes, and per genome the fragments to add when the hash fraction gives it too few
pick = []                       # per shard: boolean mask
per_genome = defaultdict(list)  # acc -> [(hash, shard, row)] of fragments not picked by the fraction
count = defaultdict(int)
for s, p in enumerate(paths):
    ids = ids_of(p)
    h = np.array([zlib.crc32(i.encode()) for i in ids], np.uint64)
    m = h < T
    pick.append(m)
    for r, (fid, hv, ok) in enumerate(zip(ids, h, m)):
        acc = acc_of(fid)
        if ok: count[acc] += 1
        else:
            c = per_genome[acc]; c.append((int(hv), s, r))
            if len(c) > 4 * a.min_per_genome: c.sort(); del c[a.min_per_genome:]   # only the lowest hashes matter
    print(f"shard {s}: {m.sum():,} of {len(ids):,} by hash", flush=True)
extra = 0
for acc, cand in per_genome.items():
    need = a.min_per_genome - count[acc]
    if need <= 0: continue
    for _, s, r in sorted(cand)[:need]: pick[s][r] = True; extra += 1

# pass 2: write the chosen rows, their taxonomy and anchor groups
D = header(paths[0])[1]
out = open(a.out + ".ref.fsb", "wb"); out.write(b"FSB1" + struct.pack("<iq", D, 0))
taxf = open(a.out + ".ref.tax.tsv", "w"); taxf.write("id\t" + "\t".join(cols) + "\n")
refids, rows = [], []
for s, p in enumerate(paths):
    f, d, n = header(p)
    m = pick[s]
    ids = ids_of(p)
    f.seek(16)
    for start in range(0, n, 200000):                     # stream rows in blocks
        stop = min(n, start + 200000)
        X = np.frombuffer(f.read(4 * d * (stop - start)), np.float32).reshape(-1, d)
        out.write(X[m[start:stop]].tobytes())
    for i in np.flatnonzero(m):
        r = tax[acc_of(ids[i])]
        refids.append(ids[i]); rows.append(r)
        taxf.write(ids[i] + "\t" + "\t".join(r) + "\n")
out.write("\n".join(refids).encode()); out.seek(8); out.write(struct.pack("<q", len(refids))); out.close()
taxf.close()
codes = []
for j, deep in ((0, 1), (1, 2), (3, 4), (4, 5)):  # domain, group, class, order, as gtdb_prep.py
    ix, c = {}, []
    for r in rows:
        key = ";".join(r[:deep]) if r[j] else None
        c.append(ix.setdefault(key, len(ix)) if key else -1)
    codes.append(c)
with open(a.out + ".ref.anchor", "w") as f:
    for row in zip(*codes): f.write(" ".join(map(str, row)) + "\n")
print(f"{len(refids):,} reference fragments ({len(refids) - extra:,} by hash < {a.ref}, {extra:,} added so every genome "
      f"has >= {a.min_per_genome}) -> {a.out}.ref.fsb")
