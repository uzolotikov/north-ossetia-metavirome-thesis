#!/usr/bin/env python3
"""Novelty and lysogeny markers of the 88 viral vOTU.

Novelty: average amino-acid identity (AAI) of each vOTU to its closest genome in
the CheckV database (aai_id in CheckV completeness.tsv of the catalogue input).
Lysogeny: geNomad gene annotations containing integrase, recombinase, excisionase,
repressor, antirepressor or the CII lysogeny regulator (PF06892); listed separately for vOTU not already found embedded
as prophages.

Output: report/votu_novelty_lysogeny.tsv and a printed summary.
"""
import csv, io, re, sys
from collections import Counter
from statistics import median

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
K = "/NatureUsers/nsinyavskiy/metaviroms_kgt"
W = "/NatureUsers/glitch/mv_review"
KEYS = ("integrase", "recombinase", "excisionase", "antirepressor", "repressor", "regulatory protein cii")

orig = {}
for line in open(K + "/votu/no_NC/catalog/votu_no_NC_95_85/representatives_ge5kb.fna"):
    if line.startswith(">"):
        p = line[1:].split(); orig[p[0]] = p[1].split("=", 1)[1]
comp = {r["contig_id"]: r for r in csv.DictReader(
    open(K + "/checkv/filtered_no_NC/strict_supported_pass2/completeness.tsv"), delimiter="\t")}
def checkv(v):
    o = orig[v]
    for k in (o, re.sub(r"_\d+$", "", o)):
        if k in comp:
            return comp[k]
    return None

rows = [l.rstrip("\n").split("\t") for l in open(W + "/report/votu_origin.tsv", encoding="utf-8") if not l.startswith("#")]
ix = {k: i for i, k in enumerate(rows[0])}
origin = {r[0]: r for r in rows[1:]}
viral = sorted(v for v, r in origin.items() if r[ix["category"]] != "0_not_virus")

lyso = {}
for r in csv.DictReader(open(W + "/genomad/votu97/representatives_ge5kb_annotate/representatives_ge5kb_genes.tsv"), delimiter="\t"):
    v = re.match(r"(vOTU_\d+)", r["gene"]).group(1)
    d = (r["annotation_description"] or "").lower()
    hit = [k for k in KEYS if k in d and not (k == "repressor" and "antirepressor" in d)]
    if hit:
        lyso.setdefault(v, set()).update(hit)

out = []
for v in viral:
    c = checkv(v)
    out.append((v, origin[v][ix["taxon"]], c["aai_id"], c["aai_af"], c["aai_confidence"], c["aai_top_hit"],
                origin[v][ix["category"]], "yes" if origin[v][ix["host_contig"]] else "no",
                ",".join(sorted(lyso.get(v, ())))))
with open(W + "/report/votu_novelty_lysogeny.tsv", "w", encoding="utf-8") as fh:
    fh.write("# AAI to the closest genome of the CheckV database; lysogeny-associated gene annotations by geNomad\n")
    fh.write("votu\ttaxon\taai_id\taai_af\taai_confidence\taai_top_hit\torigin_category\tembedded_prophage\tlysogeny_genes\n")
    for r in out:
        fh.write("\t".join(r) + "\n")

aai = sorted(float(r[2]) for r in out)
fam = Counter(r[1] for r in out)
emb = {r[0] for r in out if r[7] == "yes"}
lys_extra = {r[0] for r in out if r[8] and r[0] not in emb}
print("вирусных vOTU:", len(out))
print("таксоны:", dict(fam))
print("AAI: min %.1f, медиана %.1f, max %.1f; >=90: %d; 70-90: %d; <70: %d" % (
    aai[0], median(aai), aai[-1], sum(x >= 90 for x in aai), sum(70 <= x < 90 for x in aai), sum(x < 70 for x in aai)))
print("профаги (встроены):", len(emb), "| с генами лизогении, не встроены:", len(lys_extra), "| итого умеренных:", len(emb | lys_extra))
print("NOVELTY_DONE")
