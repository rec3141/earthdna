#!/usr/bin/env python3
"""GTDB species-representative genomes -> every 10 kb genome fragment, as inputs for a sample + project globe.

Each genome is cut into non-overlapping FRAG-bp windows from the start of every contig at least FRAG bp long; windows
with an ambiguity run (2+ consecutive non-ACGT bases) or more than 10 ambiguous bases are dropped. Every remaining
window is kept. Features are canonical k-mer counts: a k-mer and its reverse complement share one column (512 for
k = 5), summing both strands, so the profile does not depend on which strand a fragment came from; the values equal
kmers.py --rc restricted to the canonical columns. A fraction --ref of the fragments (chosen by a hash of the
fragment id, so the choice does not depend on worker count or order) is also written as the reference set that
fastsne embeds; everything is then projected onto that map.

  OUT.all.NNN.fsb   all fragments, sharded (float32 counts, fastsne binary input), listed in OUT.all_fsb.txt
  OUT.ref.fsb       the reference fragments
  OUT.ref.tax.tsv   id, domain, group, phylum, class, order, family, genus, species, lineage_rank, lineage, organism
                    for the reference (silva_prep.py's layout: group = GTDB phylum, organism = genome accession)
  OUT.ref.anchor    group ids for fastsne --anchor-file, one row per reference row: domain, group, class, order
  OUT.ref.fasta.gz  the reference fragments' sequences
  OUT.genomes.tsv   per genome: accession, contigs, bp, contigs >= FRAG, fragments, reference fragments, then the
                    taxonomy columns above (the taxonomy of any fragment, from the accession in its id)

Fragment ids are <accession>_<contig>_<start> (0-based start). Genomes are read in parallel; only counts (and the
reference sequences) leave the workers.

usage: gtdb_prep.py GTDB_RELEASE_DIR OUT [--ref 0.065] [--frag 10000] [--k 5] [--workers 16] [--shard 2000000]
                    [--shard-dir DIR] [--limit N]
  GTDB_RELEASE_DIR holds skani/genome_paths.tsv, skani/database/... and taxonomy/gtdb_taxonomy.tsv
"""
import argparse, gzip, os, struct, time, zlib
from multiprocessing import Pool
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("gtdb"); ap.add_argument("out")
ap.add_argument("--ref", type=float, default=0.065); ap.add_argument("--frag", type=int, default=10000)
ap.add_argument("--k", type=int, default=5); ap.add_argument("--workers", type=int, default=os.cpu_count())
ap.add_argument("--shard", type=int, default=2000000, help="fragments per all-fragment .fsb shard")
ap.add_argument("--shard-dir", help="write the all-fragment shards here instead of next to OUT (they are big)")
ap.add_argument("--seed", type=int, default=1)
ap.add_argument("--limit", type=int, default=0, help="only this many genomes, chosen at random (pilots)")
a = ap.parse_args()

k, F = a.k, a.frag
lut = np.full(256, 255, np.uint8)
for i, c in enumerate(b"ACGT"): lut[c] = i; lut[c + 32] = i
rc = np.array([sum((3 - ((w >> (2 * j)) & 3)) * 4 ** (k - 1 - j) for j in range(k)) for w in range(4 ** k)])
canon = np.unique(np.minimum(np.arange(4 ** k), rc))  # canonical k-mers in index order
col = np.zeros(4 ** k, np.int64); col[canon] = np.arange(len(canon)); col[rc[canon]] = np.arange(len(canon))
D = len(canon)
REF_T = int(a.ref * 2 ** 32)

def shred(job):
    acc, path = job
    with gzip.open(path, "rb") as f: data = f.read()
    ids, Xs, fa, ncontig, bp, long, nref = [], [], [], 0, 0, 0, 0
    for rec in data.split(b">")[1:]:
        head, _, body = rec.partition(b"\n")
        s = body.replace(b"\n", b"").replace(b"\r", b"")
        ncontig += 1; bp += len(s)
        if len(s) < F: continue
        long += 1
        name = head.split()[0].decode()
        nw = len(s) // F
        code = lut[np.frombuffer(s, np.uint8)[:nw * F]]
        bad = code > 3
        st = np.arange(nw) * F
        cb = np.concatenate([[0], np.cumsum(bad)])
        cr = np.concatenate([[0], np.cumsum(bad[:-1] & bad[1:])])
        okw = (cb[st + F] - cb[st] <= 10) & (cr[st + F - 1] - cr[st] == 0)
        if not okw.any(): continue
        # k-mer codes inside each kept window, skipping those that span ambiguity codes
        C = np.where(bad, 0, code).astype(np.uint16).reshape(nw, F)[okw]
        B = bad.reshape(nw, F)[okw]
        m, L = len(C), F - k + 1
        w = np.zeros((m, L), np.uint16); b = np.zeros((m, L), bool)
        for j in range(k): w = w * 4 + C[:, j:j + L]; b |= B[:, j:j + L]
        idx = col[w] + (np.arange(m) * D)[:, None]
        n = np.bincount(idx[~b], minlength=m * D).reshape(m, D)
        for wi in np.flatnonzero(okw):
            fid = f"{acc}_{name}_{wi * F}"
            ids.append(fid)
            if zlib.crc32(fid.encode()) < REF_T:
                nref += 1
                fa.append(b">" + fid.encode() + b"\n" + s[wi * F:(wi + 1) * F].upper() + b"\n")
        Xs.append(n.astype(np.float32))
    X = np.concatenate(Xs) if Xs else np.zeros((0, D), np.float32)
    isref = np.array([zlib.crc32(i.encode()) < REF_T for i in ids], bool)
    fz = gzip.compress(b"".join(fa), 6) if fa else b""  # concatenated gzip members are one valid gzip file
    return acc, ids, X, isref, fz, (ncontig, bp, long, len(ids), nref)

class Shards:
    def __init__(self, prefix, size):
        self.prefix, self.size, self.paths, self.f, self.ids = prefix, size, [], None, []
    def _close(self):
        if self.f:
            self.f.write("\n".join(self.ids).encode())
            self.f.seek(8); self.f.write(struct.pack("<q", len(self.ids))); self.f.close()
            self.f, self.ids = None, []
    def write(self, ids, X):
        if self.f is None or len(self.ids) >= self.size:
            self._close()
            p = f"{self.prefix}.{len(self.paths):03d}.fsb"
            self.paths.append(os.path.abspath(p))
            self.f = open(p, "wb"); self.f.write(b"FSB1" + struct.pack("<iq", D, 0))
        self.f.write(X.tobytes()); self.ids += ids
    def close(self): self._close()

if __name__ == "__main__":
    t0 = time.time()
    tax = {}
    for line in open(os.path.join(a.gtdb, "taxonomy", "gtdb_taxonomy.tsv")):
        g, lin = line.rstrip("\n").split("\t")
        r = dict(x.split("__", 1) for x in lin.split(";"))
        tax[g[3:] if g[:3] in ("RS_", "GB_") else g] = r
    jobs = []
    for line in open(os.path.join(a.gtdb, "skani", "genome_paths.tsv")):
        fn, d = line.split()[:2]
        acc = fn.replace("_genomic.fna.gz", "")
        if acc in tax: jobs.append((acc, os.path.join(a.gtdb, "skani", d, fn)))
    jobs.sort()
    if a.limit:
        keep = np.random.default_rng(a.seed).choice(len(jobs), min(a.limit, len(jobs)), replace=False)
        jobs = [jobs[i] for i in sorted(keep)]
    print(f"{len(jobs)} genomes with taxonomy ({len(tax)} in taxonomy), {F} bp fragments, {D} canonical {k}-mers, "
          f"reference fraction {a.ref}", flush=True)

    def taxrow(acc):
        r = tax[acc]
        ph = r.get("p", "")
        rk = next(((n, r[c]) for n, c in (("family", "f"), ("order", "o"), ("class", "c"), ("phylum", "p")) if r.get(c)),
                  ("", ""))
        return [r.get("d", ""), ph, ph, r.get("c", ""), r.get("o", ""), r.get("f", ""), r.get("g", ""), r.get("s", ""),
                rk[0], rk[1], acc]

    cols = ["domain", "group", "phylum", "class", "order", "family", "genus", "species", "lineage_rank", "lineage", "organism"]
    sp = os.path.join(a.shard_dir, os.path.basename(a.out)) if a.shard_dir else a.out
    allf = Shards(sp + ".all", a.shard)
    ref = open(a.out + ".ref.fsb", "wb"); ref.write(b"FSB1" + struct.pack("<iq", D, 0))
    taxf = open(a.out + ".ref.tax.tsv", "w"); taxf.write("id\t" + "\t".join(cols) + "\n")
    gen = open(a.out + ".genomes.tsv", "w")
    gen.write("accession\tcontigs\tbp\tcontigs_ge_frag\tfragments\tref_fragments\t" + "\t".join(cols) + "\n")
    fa = open(a.out + ".ref.fasta.gz", "wb")
    refids, rows, tot = [], [], np.zeros(5, np.int64)
    with Pool(a.workers) as pool:
        for n, (acc, ids, X, isref, fz, st) in enumerate(pool.imap(shred, jobs, chunksize=2), 1):
            r = taxrow(acc)
            if ids: allf.write(ids, X)
            if isref.any():
                ref.write(X[isref].tobytes()); fa.write(fz)
                for i in np.flatnonzero(isref):
                    taxf.write(ids[i] + "\t" + "\t".join(x.replace("\t", " ") for x in r) + "\n")
                    refids.append(ids[i]); rows.append(r)
            gen.write(acc + "\t" + "\t".join(map(str, st)) + "\t" + "\t".join(r) + "\n")
            tot += st
            if n % 5000 == 0 or n == len(jobs):
                print(f"[{time.time() - t0:7.0f}s] {n} genomes, {tot[1] / 1e9:.1f} Gbp, {tot[3]} fragments, "
                      f"{tot[4]} reference", flush=True)
    allf.close()
    ref.write("\n".join(refids).encode())
    ref.seek(8); ref.write(struct.pack("<q", len(refids))); ref.close()
    taxf.close(); gen.close(); fa.close()
    open(a.out + ".all_fsb.txt", "w").write("\n".join(allf.paths) + "\n")
    codes = []
    for j, deep in ((0, 1), (1, 2), (3, 4), (4, 5)):  # domain, group, class, order; nested so names repeat safely
        ix, c = {}, []
        for r in rows:
            key = ";".join(r[:deep]) if r[j] else None
            c.append(ix.setdefault(key, len(ix)) if key else -1)
        codes.append(c)
    with open(a.out + ".ref.anchor", "w") as f:
        for row in zip(*codes): f.write(" ".join(map(str, row)) + "\n")
    print(f"{tot[3]} fragments in {len(allf.paths)} shards, {len(refids)} reference, from {len(jobs)} genomes "
          f"({time.time() - t0:.0f}s)")
