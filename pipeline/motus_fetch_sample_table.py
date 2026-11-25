#!/usr/bin/env python3
"""motus-db.org sample table (the site's Samples page) -> one TSV: per sample its study, BioSample, BioProject, MAG
and genome counts, the six-level environment ontology and coordinates.

The page loads rows from POST https://motus-db.org/api/sample-cols a page at a time; its export buttons are
switched off, and the Zenodo samples table carries none of these columns. This asks for the same pages the
page does, in id order, one request at a time with a pause between them, and stops at the first empty page.
Rows are deduplicated by id (pages overlap by one row).

usage: python3 -I motus_fetch_sample_table.py OUT.tsv [--page 10000] [--pause 1.5]
"""
import argparse, json, sys, time, urllib.request

ap = argparse.ArgumentParser()
ap.add_argument("out"); ap.add_argument("--page", type=int, default=10000); ap.add_argument("--pause", type=float, default=1.5)
a = ap.parse_args()
COLS = ["id", "sushisample", "sushistudy", "biosample", "bioproject", "num_mags", "num_genome", "context_domain", "realm",
        "ecosystem_type", "ecosystem_name", "ecosystem_compartment", "sample_source", "latitude_degN", "longitude_degE",
        "environment", "ecosystem", "public_sample_link"]

def page(start):
    body = json.dumps(dict(repo_id="motus", global_search_value="", search_value_pair=[], start=start, length=a.page,
                           window=0, fetchall=0, order=[], order_cols=[])).encode()
    req = urllib.request.Request("https://motus-db.org/api/sample-cols", data=body, method="POST",
                                 headers={"Content-Type": "application/json", "Accept": "application/json",
                                          "User-Agent": "Mozilla/5.0 (research data client)"})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.load(r)
        except Exception as e:
            print(f"start {start}: {e}; attempt {attempt + 1}", file=sys.stderr, flush=True)
            time.sleep(30 * (attempt + 1))
    raise SystemExit(f"start {start}: failed 3 times; rerun to start over")

seen, total, start, extra = set(), None, 0, set()
clean = lambda v: "" if v is None else str(v).replace("\t", " ").replace("\n", " ")
with open(a.out, "w") as o:
    o.write("\t".join(COLS) + "\n")
    while True:
        d = page(start)
        total = total or d.get("total_count")
        rows = d.get("data") or []
        new = [r for r in rows if r.get("id") not in seen]
        if not new: break
        for r in new:
            seen.add(r["id"]); extra |= set(r) - set(COLS)
            o.write("\t".join(clean(r.get(c)) for c in COLS) + "\n")
        print(f"start {start}: {len(rows)} rows, {len(new)} new, {len(seen)} of {total}", flush=True)
        start += a.page
        time.sleep(a.pause)
print(f"done: {len(seen)} rows (site total_count {total}); fields not written: {sorted(extra) or 'none'}", flush=True)
