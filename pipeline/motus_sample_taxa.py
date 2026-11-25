#!/usr/bin/env python3
"""Which taxa each metagenome contains, for searching the mOTUs sample globes by any taxon (not only the one that
dominates), and how abundant its dominant taxon is, for the tooltip.

  OUT.contains.npz  the samples of SAMPLES.tax.tsv (array "samples") and, for each rank (phylum .. species), the
                    taxa present in at least one sample at >= MIN relative abundance ("<rank>_names") with the
                    samples carrying each, as CSR over sample rows ("<rank>_indptr", "<rank>_indices"), and the
                    taxon's summed relative abundance in each of them ("<rank>_data")
  OUT.dompct.tsv    id, percent: the relative abundance of the sample's dominant species as SAMPLES.tax.tsv names it
                    (or of its dominant genus when no species is named), as a percent of all its inserts
  OUT.rankpct.tsv   id, then for each rank (domain .. species) the percent of the sample's inserts in the taxon
                    SAMPLES.tax.tsv names at that rank (blank where it names none), for the tooltip of a drilled-down
                    globe, whose points are coloured by their taxon one rank below the picked one
  OUT.info.tsv      (rewritten when present) gains a Dominant column: "34% Escherichia coli"

Relative abundance is over all of a sample's inserts, the unassigned bin included, as motus profile reports it.

usage: motus_sample_taxa.py RAW_DIR MOTU.tax.tsv SAMPLES.tax.tsv OUT [--min 0.001]
"""
import argparse, os
import numpy as np
import scipy.sparse as sp

ap = argparse.ArgumentParser()
ap.add_argument("raw"); ap.add_argument("motu_tax"); ap.add_argument("samples_tax"); ap.add_argument("out")
ap.add_argument("--min", type=float, default=0.001)
a = ap.parse_args()

X = sp.load_npz(f"{a.raw}/motus41.profiles_sparse.npz").tocsc().astype(np.float64)
motus = np.loadtxt(f"{a.raw}/motus41.profiles_motu_index.csv.gz", dtype=str)
smp = np.loadtxt(f"{a.raw}/motus41.profiles_sample_index.csv.gz", dtype=str)
RA = (X @ sp.diags(1 / np.asarray(X.sum(0)).ravel())).tocsr()
mt = [l.rstrip("\n").split("\t") for l in open(a.motu_tax)]
mh, mt = mt[0], mt[1:]
st = [l.rstrip("\n").split("\t") for l in open(a.samples_tax)]
sh, st = st[0], st[1:]
row = {m: i for i, m in enumerate(motus)}
col = {s: j for j, s in enumerate(smp)}
ids = [r[0] for r in st]
A = RA[np.array([row[r[0]] for r in mt])][:, np.array([col[s] for s in ids])].T.tocsr()   # sample x mOTU

RANKS = ["phylum", "class", "order", "family", "genus", "species"]
out = {"samples": np.array(ids), "min": np.array(a.min)}
dom = {}
for rk in ["domain"]:                                              # shares only: no search by domain
    lab = np.array([r[mh.index(rk)] for r in mt])
    names, code = np.unique(lab, return_inverse=True)
    ok = lab != ""
    M = sp.csr_matrix((np.ones(ok.sum()), (np.flatnonzero(ok), code[ok])), shape=(len(mt), len(names)))
    dom[rk] = ((A @ M).tocsc(), {n: k for k, n in enumerate(names)})
for rk in RANKS:
    lab = np.array([r[mh.index(rk)] for r in mt])
    names, code = np.unique(lab, return_inverse=True)
    ok = lab != ""
    M = sp.csr_matrix((np.ones(ok.sum()), (np.flatnonzero(ok), code[ok])), shape=(len(mt), len(names)))
    T = (A @ M).tocsc()                                            # sample x taxon
    dom[rk] = (T.copy(), {n: k for k, n in enumerate(names)})
    T.data[T.data < a.min] = 0
    T.eliminate_zeros()
    keep = np.flatnonzero((np.diff(T.indptr) > 0) & (names != ""))
    T = T[:, keep].tocsc()
    out[f"{rk}_names"] = names[keep]
    out[f"{rk}_indptr"] = T.indptr.astype(np.int64)
    out[f"{rk}_indices"] = T.indices.astype(np.int32)
    out[f"{rk}_data"] = T.data.astype(np.float32)                 # relative abundance, aligned with indices
    print(f"{rk}: {len(keep)} taxa present at >= {a.min:.1%} somewhere, {T.nnz} sample-taxon pairs")
np.savez_compressed(a.out + ".contains.npz", **out)

# the dominant species' (else genus') share, matching the name samples.tax.tsv shows
si, gi = sh.index("species"), sh.index("genus")
pct = []
for k, r in enumerate(st):
    rk, name = ("species", r[si]) if r[si] else ("genus", r[gi])
    T, ix = dom[rk]
    pct.append(100 * T[k, ix[name]] if name in ix else float("nan"))
with open(a.out + ".dompct.tsv", "w") as f:
    f.write("id\tpercent\n")
    f.writelines(f"{s}\t{p:.1f}\n" for s, p in zip(ids, pct) if p == p)
ALL = ["domain"] + RANKS
cols = [sh.index(rk) for rk in ALL]
rows_of = {rk: dom[rk][0].tocsr() for rk in ALL}
with open(a.out + ".rankpct.tsv", "w") as f:
    f.write("id\t" + "\t".join(ALL) + "\n")
    for k, r in enumerate(st):
        v = []
        for rk, c in zip(ALL, cols):
            ix, R = dom[rk][1], rows_of[rk]
            v.append(f"{100 * R[k, ix[r[c]]]:.2f}" if r[c] and r[c] in ix else "")
        f.write(r[0] + "\t" + "\t".join(v) + "\n")
info = a.out + ".info.tsv"
if os.path.exists(info):
    lines = [l.rstrip("\n").split("\t") for l in open(info)]
    if "Dominant" not in lines[0]:
        by = {s: (p, r[si] or (r[gi] + " sp." if r[gi] else "")) for s, p, r in zip(ids, pct, st)}
        with open(info, "w") as f:
            f.write("\t".join(lines[0] + ["Dominant"]) + "\n")
            for l in lines[1:]:
                p, n = by.get(l[0], (float("nan"), ""))
                f.write("\t".join(l + [f"{p:.0f}% {n}" if p == p and n else ""]) + "\n")
print(f"dominant taxon share: median {np.nanmedian(pct):.0f}% over {len(ids)} samples")
