#!/usr/bin/env python3
"""mOTUs 4.1 taxonomic profiles -> per-mOTU prevalence, and a traits file for EarthDNA's pack_globe.py --traits.

Reads the sparse COO matrix (motus41.profiles_sparse.npz, mOTUs x samples, insert-scaled counts) with numpy alone
(allow_pickle=False; scipy is not needed) plus its two sorted index files. A sample's study is the part of its id
before the first "_" (e.g. ACIN21-1_SAMN05421555_METAG -> ACIN21-1).

  OUT.prevalence.tsv  motu, samples detected, studies detected, total count
  OUT.traits.tsv      header "organism" + traits, one row per mOTU, 0/1 per trait; keyed on the organism column
                      that motus_prep.py writes (the mOTU id):
                        profiled      detected in at least one of the 115,013 samples
                        widespread    detected in >= 1% of samples
                        multi_study   detected in >= 10 studies
                        rare          detected in 1-10 samples
                      with --genomes (db_mOTU/mOTUsv4.1.genomes.tsv.gz, its MOTU4 column), two more:
                        isolate       the mOTU holds at least one isolate or single-cell genome (_GENO_ ids)
                        mag_only      every genome of the mOTU is a MAG (_MAG_ ids)
                      mOTUs that the genome table names but no sample detects are written with profiled = 0.

usage: python3 -I motus_profile_traits.py profiles_sparse.npz motu_index.csv.gz sample_index.csv.gz OUT
                                          [--genomes mOTUsv4.1.genomes.tsv.gz]
"""
import argparse, gzip
import numpy as np

ap = argparse.ArgumentParser()
for x in ("npz", "motu_ix", "sample_ix", "out"): ap.add_argument(x)
ap.add_argument("--genomes")
A = ap.parse_args()
npz, motu_ix, sample_ix, out = A.npz, A.motu_ix, A.sample_ix, A.out
z = np.load(npz, allow_pickle=False)
if z["format"].item() != b"coo":
    raise SystemExit(f"expected a COO matrix, found {z['format'].item()!r}")
row, col, data = z["row"], z["col"], z["data"]
n_motu, n_sample = (int(x) for x in z["shape"])
with gzip.open(motu_ix, "rt") as f: motus = [l.strip() for l in f if l.strip()]
with gzip.open(sample_ix, "rt") as f: samples = [l.strip() for l in f if l.strip()]
if (len(motus), len(samples)) != (n_motu, n_sample):
    raise SystemExit(f"index lengths {len(motus)}, {len(samples)} do not match matrix shape {n_motu} x {n_sample}")
keep = data > 0
row, col, data = row[keep], col[keep], data[keep]

n_samples = np.bincount(row, minlength=n_motu)
total = np.bincount(row, weights=data.astype(np.float64), minlength=n_motu)
study_names, study_of_sample = np.unique([s.split("_", 1)[0] for s in samples], return_inverse=True)
pair = np.unique(row.astype(np.int64) * len(study_names) + study_of_sample[col])
n_studies = np.bincount(pair // len(study_names), minlength=n_motu)

wide = n_samples >= 0.01 * n_sample
with open(out + ".prevalence.tsv", "w") as f:
    f.write("motu\tsamples\tstudies\ttotal_count\n")
    for m, s, st, t in zip(motus, n_samples, n_studies, total):
        f.write(f"{m}\t{s}\t{st}\t{int(t)}\n")
iso, mag = set(), set()
if A.genomes:
    with gzip.open(A.genomes, "rt") as f:
        h = next(f).rstrip("\n").split("\t")
        ig, im = h.index("GENOME"), h.index("MOTU4")
        for line in f:
            r = line.rstrip("\n").split("\t")
            if r[im] == "no_mOTU": continue
            (mag if "_MAG_" in r[ig] else iso).add(r[im])
    extra = sorted((iso | mag) - set(motus))
    motus = motus + extra
    pad = lambda v: np.concatenate([v, np.zeros(len(extra), v.dtype)])
    n_samples, n_studies, wide = pad(n_samples), pad(n_studies), pad(wide)
traits = [("profiled", n_samples > 0), ("widespread", wide), ("multi_study", n_studies >= 10),
          ("rare", (n_samples > 0) & (n_samples <= 10))]
if A.genomes:
    traits += [("isolate", np.array([m in iso for m in motus])),
               ("mag_only", np.array([m in mag and m not in iso for m in motus]))]
with open(out + ".traits.tsv", "w") as f:
    f.write("organism\t" + "\t".join(n for n, _ in traits) + "\n")
    for i, m in enumerate(motus):
        if m.endswith("_unassigned"): continue
        f.write(m + "\t" + "\t".join(str(int(v[i])) for _, v in traits) + "\n")
print(f"{n_motu} mOTUs x {n_sample} samples ({len(study_names)} studies), {len(data)} non-zero counts")
print(f"  {len(motus) - n_motu} mOTUs from the genome table not in the profiles" if A.genomes else "")
print("  " + ", ".join(f"{n} {int(v.sum())}" for n, v in traits))
print(f"  -> {out}.prevalence.tsv, {out}.traits.tsv")
