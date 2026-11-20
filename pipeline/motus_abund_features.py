#!/usr/bin/env python3
"""mOTUs 4.1 sample x mOTU profiles -> fastsne inputs for the two abundance views of the EarthDNA mOTUs globes.

  OUT_MOTU.fsb       one row per mOTU of TAX (that is profiled), in TAX order: the mOTU's relative abundance
                     across all samples (t(sample matrix)), Hellinger (square root), scaled to unit length so that
                     the pattern of occurrence matters and not overall abundance, then a 128-d sparse SVD
  OUT_MOTU.tax.tsv   TAX restricted to the same mOTUs (the shared point set of the three mOTU layouts)
  OUT_SAMPLE.fsb     one row per sample: its mOTU composition (relative abundance over assigned mOTUs, the
                     mOTUv4.0_unassigned bin left out), Hellinger, then a 128-d sparse SVD
  OUT_SAMPLE.tax.tsv EarthDNA taxonomy columns filled with each sample's most abundant taxon at each rank
                     (summed relative abundance of the mOTUs carrying that label); organism = sample id
  OUT_SAMPLE.study.tsv  sample -> study
  OUT_SAMPLE.traits.tsv  per-sample traits for pack_globe.py --traits

Relative abundance divides by the sample's total insert count including the unassigned bin, as motus profile does.
SVD scores are signed: run fastsne on these with --transform none.

usage: motus_abund_features.py RAW_DIR TAX.tsv OUT_MOTU OUT_SAMPLE [--dims 128]
  RAW_DIR  holds motus41.profiles_sparse.npz, its two index files and mOTUs4.1.samples.tsv.gz
  TAX.tsv  motus_prep.py's taxonomy table (one row per mOTU)
"""
import argparse, csv, gzip, struct, time
import numpy as np
import scipy.sparse as sp
from sklearn.decomposition import TruncatedSVD

ap = argparse.ArgumentParser()
ap.add_argument("raw"); ap.add_argument("tax"); ap.add_argument("out_motu"); ap.add_argument("out_sample")
ap.add_argument("--dims", type=int, default=128)
a = ap.parse_args()
t0 = time.time()
log = lambda m: print(f"[{time.time() - t0:6.1f}s] {m}", flush=True)

X = sp.load_npz(f"{a.raw}/motus41.profiles_sparse.npz").tocsc().astype(np.float32)   # mOTU x sample counts
motus = np.loadtxt(f"{a.raw}/motus41.profiles_motu_index.csv.gz", dtype=str)
samples = np.loadtxt(f"{a.raw}/motus41.profiles_sample_index.csv.gz", dtype=str)
assert X.shape == (len(motus), len(samples)), (X.shape, len(motus), len(samples))
log(f"profiles {X.shape[0]} mOTUs x {X.shape[1]} samples, {X.nnz} non-zero")

total = np.asarray(X.sum(0)).ravel()
RA = (X @ sp.diags(1 / total)).tocsr()                      # relative abundance, columns sum to 1
una = np.flatnonzero(motus == "mOTUv4.0_unassigned")
una_frac = RA[una].toarray().ravel() if len(una) else np.zeros(len(samples))

with open(a.tax) as f:
    header = f.readline().rstrip("\n").split("\t")
    tax = [l.rstrip("\n").split("\t") for l in f]
row_of = {m: i for i, m in enumerate(motus)}
keep = [r for r in tax if r[0] in row_of]
log(f"taxonomy {len(tax)} mOTUs, {len(keep)} profiled (kept), {len(tax) - len(keep)} never detected (dropped)")
rows = np.array([row_of[r[0]] for r in keep])

def write_fsb(path, M, ids):
    M = np.ascontiguousarray(M, dtype=np.float32)
    assert np.isfinite(M).all() and M.shape[0] == len(ids)
    with open(path, "wb") as o:
        o.write(b"FSB1" + struct.pack("<iq", M.shape[1], M.shape[0])); o.write(M.tobytes()); o.write("\n".join(ids).encode())

def svd(M, what):
    s = TruncatedSVD(a.dims, algorithm="randomized", n_iter=7, random_state=0)
    Y = s.fit_transform(M)
    log(f"{what}: SVD {M.shape} -> {a.dims} d, explained variance {s.explained_variance_ratio_.sum():.3f}")
    return Y

# mOTU view: rows of t(sample matrix)
H = RA[rows].sqrt()
norm = np.sqrt(np.asarray(H.multiply(H).sum(1)).ravel())
H = sp.diags(1 / norm) @ H
write_fsb(a.out_motu + ".fsb", svd(H, "mOTUs"), [r[0] for r in keep])
with open(a.out_motu + ".tax.tsv", "w") as f:
    f.write("\t".join(header) + "\n"); f.writelines("\t".join(r) + "\n" for r in keep)

# sample view: rows of the sample matrix over assigned mOTUs
A = RA.T.tocsr()[:, rows]                                   # sample x assigned mOTU (in TAX order)
asum = np.asarray(A.sum(1)).ravel()
empty = asum == 0
log(f"samples with no assigned mOTU (dropped): {int(empty.sum())}")
S = (sp.diags(1 / np.where(empty, 1, asum)) @ A).tocsr()[~empty]
smp = samples[~empty]
write_fsb(a.out_sample + ".fsb", svd(S.sqrt(), "samples"), list(smp))

# dominant taxon per rank: summed relative abundance of the mOTUs carrying each label
C = {h: i for i, h in enumerate(header)}
RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
dom = {}
for rk in RANKS:
    lab = np.array([r[C[rk]] for r in keep])
    names, code = np.unique(lab, return_inverse=True)
    ok = lab != ""
    M = sp.csr_matrix((np.ones(ok.sum(), np.float32), (np.flatnonzero(ok), code[ok])), shape=(len(keep), len(names)))
    T = (S @ M).tocsr()
    best = np.asarray(T.argmax(1)).ravel()
    has = np.asarray(T.max(1).todense()).ravel() > 0
    dom[rk] = np.where(has, names[best], "")
fam_ok = dom["family"] != ""
with open(a.out_sample + ".tax.tsv", "w") as f:
    f.write("\t".join(header) + "\n")
    for i, s in enumerate(smp):
        lr, ln = next(((rk, dom[rk][i]) for rk in ("family", "order", "class", "phylum") if dom[rk][i]), ("", ""))
        rec = dict(id=s, group=dom["phylum"][i], lineage_rank=lr, lineage=ln, organism=s, **{rk: dom[rk][i] for rk in RANKS})
        f.write("\t".join(rec.get(h, "") for h in header) + "\n")

meta = {}
with gzip.open(f"{a.raw}/mOTUs4.1.samples.tsv.gz", "rt") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        meta[r["#SAMPLE"]] = r
missing = [s for s in smp if s not in meta]
log(f"samples missing from mOTUs4.1.samples.tsv.gz: {len(missing)} {missing[:3]}")
with open(a.out_sample + ".study.tsv", "w") as f:
    f.write("sample\tstudy\n"); f.writelines(f"{s}\t{meta.get(s, {}).get('STUDY', '')}\n" for s in smp)

arch = np.asarray((S @ sp.csr_matrix(np.array([[r[C["domain"]] == "Archaea"] for r in keep], np.float32))).todense()).ravel()
una_k, tot_k = una_frac[~empty], total[~empty]
traits = dict(added_4_1=np.array([meta.get(s, {}).get("ADDED") == "4.1" for s in smp]),
              mostly_unassigned=una_k > 0.5, archaea_1pct=arch >= 0.01, low_depth=tot_k < 1000)
with open(a.out_sample + ".traits.tsv", "w") as f:
    f.write("organism\t" + "\t".join(traits) + "\n")
    for i, s in enumerate(smp):
        f.write(s + "\t" + "\t".join("1" if v[i] else "0" for v in traits.values()) + "\n")
log("sample traits: " + ", ".join(f"{k} {int(v.sum())}" for k, v in traits.items()))
log(f"samples whose dominant family is known: {int(fam_ok.sum())} of {len(smp)}")
log(f"wrote {a.out_motu}.fsb ({len(keep)} rows), {a.out_sample}.fsb ({len(smp)} rows)")
