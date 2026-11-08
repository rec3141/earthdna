#!/usr/bin/env python3
"""GTDB genome metadata -> per-genome traits for pack_globe.py --traits (keyed on the organism column = accession).

  cultured       an isolate genome: NCBI genome category is not "derived from metagenome / single cell /
                 environmental sample" (MAGs and SAGs are not cultured)
  named          the GTDB species has a Latin epithet (incl. GTDB's split names such as coli_A), not a placeholder
                 such as sp910578175
  type_material  NCBI type material: the type strain of a species, a subspecies or a heterotypic synonym

Only GTDB species representatives are written (one per species, as on the globe).

usage: gtdb_traits.py bac120_metadata_r226.tsv.gz ar53_metadata_r226.tsv.gz out.tsv
"""
import gzip, re, sys

*meta, out = sys.argv[1:]
epithet = re.compile(r"^[a-z]+(_[A-Z]+)?$")
n = [0, 0, 0, 0]
with open(out, "w") as o:
    o.write("organism\tcultured\tnamed\ttype_material\n")
    for path in meta:
        with gzip.open(path, "rt") as f:
            h = next(f).rstrip("\n").split("\t")
            c = {k: h.index(k) for k in ("accession", "gtdb_representative", "gtdb_taxonomy", "ncbi_genome_category",
                                         "gtdb_type_designation_ncbi_taxa")}
            for line in f:
                r = line.rstrip("\n").split("\t")
                if r[c["gtdb_representative"]] != "t": continue
                acc = r[c["accession"]].split("_", 1)[1]                       # RS_GCF_... / GB_GCA_... -> GCF_...
                sp = r[c["gtdb_taxonomy"]].rsplit("s__", 1)[-1].split(" ")
                cult = r[c["ncbi_genome_category"]] in ("none", "")
                named = len(sp) > 1 and bool(epithet.match(sp[1]))
                typ = r[c["gtdb_type_designation_ncbi_taxa"]].startswith("type strain")
                o.write(f"{acc}\t{int(cult)}\t{int(named)}\t{int(typ)}\n")
                n[0] += 1; n[1] += cult; n[2] += named; n[3] += typ
print(f"{n[0]} species representatives: {n[1]} cultured, {n[2]} named, {n[3]} type material -> {out}")
