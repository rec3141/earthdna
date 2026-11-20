#!/usr/bin/env python3
"""motus-db.org sample environments -> environment label tables for the two mOTUs globes (pack_globe.py --labels).

  OUT_SAMPLE.env.tsv  per sample (the rows of SAMPLES.tax.tsv): context domain, realm, ecosystem type, ecosystem,
                      compartment, sample source, as motus-db.org curates them
  OUT_SAMPLE.coords.tsv  per sample: latitude, longitude (empty when the site has none)
  OUT_MOTU.env.tsv    per mOTU (the rows of MOTU.tax.tsv): the environment holding most of its abundance, chosen top
                      down so the levels nest: the context domain with the largest summed relative abundance over the
                      mOTU's samples, then the largest realm inside that domain, and so on to the sample source

"NA" and empty values are left blank. A name that sits under two parents at one level (3 sample sources) gets the
parent in brackets, so each name has one place in the hierarchy.

usage: motus_env_labels.py SITE_SAMPLES.tsv RAW_DIR SAMPLES.tax.tsv MOTU.tax.tsv OUT_SAMPLE OUT_MOTU
  SITE_SAMPLES.tsv  fetch_sample_table.py's dump of the motus-db.org Samples table
  RAW_DIR           holds motus41.profiles_sparse.npz and its two index files
"""
import collections, csv, sys
import numpy as np
import scipy.sparse as sp

site, raw, s_tax, m_tax, out_s, out_m = sys.argv[1:7]
L = ["context_domain", "realm", "ecosystem_type", "ecosystem_name", "ecosystem_compartment", "sample_source"]
HEAD = ["context domain", "realm", "ecosystem type", "ecosystem", "compartment", "sample source"]
blank = lambda v: "" if v in ("", "NA", "None", "nan") else v

rows = {r["sushisample"]: r for r in csv.DictReader(open(site), delimiter="\t")}
samples = [l.split("\t", 1)[0] for l in list(open(s_tax))[1:]]
missing = [s for s in samples if s not in rows]
assert not missing, f"{len(missing)} samples not in the site table, e.g. {missing[:3]}"
path = [[blank(rows[s][c]) for c in L] for s in samples]
for p in path:  # a blank level blanks everything below it, so the levels nest
    for k in range(1, len(L)):
        if not p[k - 1]: p[k] = ""
parents = [collections.defaultdict(set) for _ in L]
for p in path:
    for k in range(1, len(L)):
        if p[k]: parents[k][p[k]].add(p[k - 1])
for p in path:
    for k in range(1, len(L)):
        if p[k] and len(parents[k][p[k]]) > 1: p[k] = f"{p[k]} ({p[k - 1]})"
with open(out_s + ".env.tsv", "w") as f:
    f.write("id\t" + "\t".join(HEAD) + "\n")
    f.writelines(s + "\t" + "\t".join(p) + "\n" for s, p in zip(samples, path))
with open(out_s + ".coords.tsv", "w") as f:
    f.write("id\tlat\tlon\n")
    for s in samples:
        la, lo = blank(rows[s]["latitude_degN"]), blank(rows[s]["longitude_degE"])
        f.write(f"{s}\t{la}\t{lo}\n")
print(f"samples: {len(samples)}, with a context domain {sum(1 for p in path if p[0])}, "
      f"with coordinates {sum(1 for s in samples if blank(rows[s]['latitude_degN']))}")

# mOTUs: summed relative abundance over each sample's environment path, chosen top down
X = sp.load_npz(f"{raw}/motus41.profiles_sparse.npz").tocsc().astype(np.float64)
motus = np.loadtxt(f"{raw}/motus41.profiles_motu_index.csv.gz", dtype=str)
allsmp = np.loadtxt(f"{raw}/motus41.profiles_sample_index.csv.gz", dtype=str)
RA = (X @ sp.diags(1 / np.asarray(X.sum(0)).ravel())).tocsr()
col = {s: j for j, s in enumerate(allsmp)}
cols = np.array([col[s] for s in samples])
mrow = {m: i for i, m in enumerate(motus)}
mids = [l.split("\t", 1)[0] for l in list(open(m_tax))[1:]]
A = RA[np.array([mrow[m] for m in mids])][:, cols].tocsr()      # mOTU x labelled sample
choice = [[""] * len(L) for _ in mids]
prefix_of_motu = np.zeros(len(mids), np.int64)                   # code of the chosen path so far (level k-1)
prev_prefixes = {(): 0}
for k in range(len(L)):
    pre = [tuple(p[:k + 1]) if all(p[:k + 1]) else None for p in path]
    keys = sorted({x for x in pre if x is not None})
    ix = {x: i for i, x in enumerate(keys)}
    ok = np.array([x is not None for x in pre])
    I = sp.csr_matrix((np.ones(ok.sum()), (np.flatnonzero(ok), [ix[x] for x in pre if x is not None])), shape=(len(samples), len(keys)))
    T = (A @ I).toarray()                                         # mOTU x path prefix of length k+1
    parent = np.array([prev_prefixes.get(x[:k], -1) for x in keys])
    T[parent[None, :] != prefix_of_motu[:, None]] = 0             # only children of the path chosen so far
    best, has = T.argmax(1), T.max(1) > 0
    for i in np.flatnonzero(has): choice[i][k] = keys[best[i]][k]
    prefix_of_motu = np.where(has, best, -2)
    prev_prefixes = ix
    print(f"mOTUs labelled at {HEAD[k]}: {int(has.sum())} of {len(mids)}, {len(set(keys[b] for b in best[has]))} categories")
with open(out_m + ".env.tsv", "w") as f:
    f.write("id\t" + "\t".join(HEAD) + "\n")
    f.writelines(m + "\t" + "\t".join(c) + "\n" for m, c in zip(mids, choice))
