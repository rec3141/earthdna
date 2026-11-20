#!/usr/bin/env python3
"""motus-db.org sample environments -> environment label tables for the two mOTUs globes (pack_globe.py --labels).

  OUT_SAMPLE.env.tsv  per sample (the rows of SAMPLES.tax.tsv): context domain, realm, ecosystem type, ecosystem,
                      compartment, sample source, as motus-db.org curates them
  OUT_SAMPLE.coords.tsv  per sample: latitude, longitude (empty when the site has none)
  OUT_SAMPLE.geo.tsv  the samples that have coordinates, as a map layout (id, lon, lat). Samples that share exact
                      coordinates are scattered around them, each by a fixed amount, within 0.02 * sqrt(n) degrees
                      for n samples at that spot: a lone sample stays within about 2 km, 1,000 at one spot spread
                      over about 0.6 degrees
  OUT_SAMPLE.geo.tax.tsv  SAMPLES.tax.tsv restricted to those samples
  OUT_SAMPLE.info.tsv  per sample, for the viewer's tooltip while the environment labels are on: sample, study,
                      BioSample, lat, lon
  OUT_MOTU.info.tsv   per mOTU, the same: mOTU, its GTDB species (or genus), the number of samples and studies it
                      is detected in
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
def coords(r):
    """the site's latitude, longitude as given, or ("", "") when either is missing or out of range (27 DRAG22-1
    samples carry latitudes near -170 and longitudes near -350)"""
    la, lo = blank(r["latitude_degN"]), blank(r["longitude_degE"])
    try:
        return (la, lo) if la and lo and abs(float(la)) <= 90 and abs(float(lo)) <= 180 else ("", "")
    except ValueError:
        return "", ""

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
        la, lo = coords(rows[s])
        f.write(f"{s}\t{la}\t{lo}\n")
import zlib
s_head, s_rows = open(s_tax).readline(), {l.split("\t", 1)[0]: l for l in list(open(s_tax))[1:]}
at_spot = collections.Counter(coords(rows[s]) for s in samples)
with open(out_s + ".geo.tsv", "w") as g, open(out_s + ".geo.tax.tsv", "w") as t:
    g.write("id\tlon\tlat\n"); t.write(s_head)
    for s in samples:
        la, lo = coords(rows[s])
        if not la: continue
        h, rmax = zlib.crc32(s.encode()), 0.02 * np.sqrt(at_spot[(la, lo)])
        ang, rad = (h & 0xffff) / 65536 * 2 * np.pi, rmax * np.sqrt((h >> 16) / 65536)
        la2 = max(-89.9, min(89.9, float(la) + rad * np.sin(ang)))
        lo2 = (float(lo) + rad * np.cos(ang) / max(0.05, np.cos(np.radians(la2))) + 180) % 360 - 180
        g.write(f"{s}\t{lo2:.6f}\t{la2:.6f}\n"); t.write(s_rows[s])
with open(out_s + ".info.tsv", "w") as f:
    f.write("id\tSample\tStudy\tBioSample\tlat\tlon\n")
    for s in samples:
        r = rows[s]
        f.write("\t".join([s, s, r["sushistudy"], r["biosample"], *coords(r)]) + "\n")
print(f"samples: {len(samples)}, with a context domain {sum(1 for p in path if p[0])}, "
      f"with coordinates {sum(1 for s in samples if coords(rows[s])[0])}")

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
D = X.tocsr()[np.array([mrow[m] for m in mids])]                 # every profiled sample, for the detection counts
study_of = np.array([s.split("_", 1)[0] for s in allsmp])
tax = {l.split("\t")[0]: l.rstrip("\n").split("\t") for l in list(open(m_tax))[1:]}
th = open(m_tax).readline().rstrip("\n").split("\t")
gi, si = th.index("genus"), th.index("species")
with open(out_m + ".info.tsv", "w") as f:
    f.write("id\tmOTU\ttaxon\tSamples\tStudies\n")
    for i, m in enumerate(mids):
        js = D.indices[D.indptr[i]:D.indptr[i + 1]]
        t = tax[m]
        name = t[si] or (t[gi] + " sp." if t[gi] else "")
        f.write(f"{m}\t{m}\t{name}\t{len(js)}\t{len(set(study_of[js]))}\n")
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
