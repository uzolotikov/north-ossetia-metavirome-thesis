#!/usr/bin/env python3
"""Origin of every vOTU: associated with the MV-specific contaminant or not.

Categories, assigned in priority order:
  0 not_virus        geNomad calls it a plasmid or not viral
  1 provirus_MV      embedded in a contig of an MV-specific bacterium;
                     split by host assembly: MV-4128 (Sphingomonas genome) / MV-4122 community
  2 crispr_MV        target of CRISPR spacers of an MV-specific bacterium
  3 follows_sph      Tib2/Biragzang RPM ratio inside the band of the Sphingomonas
                     contaminant (0.339-0.693), no reliable cellular confirmation
  4 native           reliably confirmed by the cellular fraction of its own site
  5 unknown          none of the above

Outputs: report/votu_origin.tsv and a printed summary.
"""
import csv
import io
import sys
from collections import defaultdict, OrderedDict

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
K = "/NatureUsers/nsinyavskiy/metaviroms_kgt"
W = "/NatureUsers/glitch/mv_review"
REP = W + "/report"
T1, T2, BI = "MV-4122_S12", "MV-4125_S18", "MV-4128_S29"
MV = [T1, T2, BI]
NAME = {T1: "Tib1", T2: "Tib2", BI: "Biragzang"}
BAND = (0.339, 0.693)
CRISPR = {"vOTU_008", "vOTU_017"}

rpm, det, length = {}, {}, {}
with open(K + "/qc/votu97_comparison/votu97_all_libraries.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        rpm[(r["sample"], r["votu"])] = float(r["RPM"])
        det[(r["sample"], r["votu"])] = int(r["detected_75pct_1x"])
        length[r["votu"]] = int(r["length_bp"])
VOTUS = sorted(length)

flags = {}
with open(REP + "/votu97_catalog_flags.tsv") as fh:
    rows = [l.rstrip("\n").split("\t") for l in fh if not l.startswith("#")]
hdr = rows[0]
for r in rows[1:]:
    flags[r[0]] = dict(zip(hdr, r))

# full geNomad taxonomy
fam = {}
with open(W + "/genomad/votu97/representatives_ge5kb_summary/representatives_ge5kb_virus_summary.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        tax = [t for t in r["taxonomy"].split(";") if t]
        fam[r["seq_name"].split("|")[0]] = tax[-1] if tax else "NA"

# CheckV quality via original contig id
orig = {}
with open(K + "/votu/no_NC/catalog/votu_no_NC_95_85/representatives_ge5kb.fna") as fh:
    for line in fh:
        if line.startswith(">"):
            p = line[1:].split()
            orig[p[0]] = p[1].split("=", 1)[1]
cv = {}
with open(K + "/checkv/filtered_no_NC/strict_supported_pass2/quality_summary.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        cv[r["contig_id"]] = (r["checkv_quality"], r["completeness"])
def checkv(v):
    o = orig.get(v, "")
    for key in (o, o.rsplit("_", 1)[0] if o.endswith(("_1", "_2", "_3")) else o):
        if key in cv:
            return cv[key]
    return ("NA", "NA")

host = {}
with open(W + "/prophage/embedded.tsv") as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2:
            host[p[0]] = p[1]

out = []
for v in VOTUS:
    f = flags[v]
    g = f["genomad"]
    rel = f["cellular"] == "reliable"
    shared = det[(T2, v)] == 1 and det[(BI, v)] == 1 and rpm[(BI, v)] > 0
    ratio = rpm[(T2, v)] / rpm[(BI, v)] if shared else None
    if g != "virus":
        cat, sub = "0_not_virus", ("plasmid" if g == "plasmid" else "not_viral")
    elif v in host:
        cat = "1_provirus_MV"
        sub = "Sphingomonas_MV4128" if host[v].startswith("MV-4128") else "community_MV4122"
        sub += "_" + f["embedded"]
    elif v in CRISPR:
        cat, sub = "2_crispr_MV", ""
    elif ratio is not None and BAND[0] <= ratio <= BAND[1] and not rel:
        cat, sub = "3_follows_sph", ""
    elif rel:
        cat = "4_native"
        sub = "site_" + f["cellular_site"]
        if ratio is not None and BAND[0] <= ratio <= BAND[1]:
            sub += "_carried_with_contaminant"
    else:
        cat, sub = "5_unknown", ""
    q, comp = checkv(v)
    dets = "".join("+" if det[(m, v)] else "-" for m in MV)
    out.append(OrderedDict([
        ("votu", v), ("length", length[v]), ("category", cat), ("detail", sub),
        ("genomad", g), ("taxon", fam.get(v, "NA")), ("checkv_quality", q), ("checkv_completeness", comp),
        ("detected_T1_T2_Bi", dets), ("ratio_T2_Bi", "%.3f" % ratio if ratio is not None else ""),
        ("host_contig", host.get(v, "")),
        ("rpm_Tib1", "%.2f" % rpm[(T1, v)]), ("rpm_Tib2", "%.2f" % rpm[(T2, v)]), ("rpm_Biragzang", "%.2f" % rpm[(BI, v)]),
    ]))

with open(REP + "/votu_origin.tsv", "w", encoding="utf-8") as fh:
    fh.write("# origin of every vOTU; see 13_virus_origin.py docstring for category rules\n")
    fh.write("\t".join(out[0].keys()) + "\n")
    for r in out:
        fh.write("\t".join(str(x) for x in r.values()) + "\n")

cats = ["1_provirus_MV", "2_crispr_MV", "3_follows_sph", "4_native", "5_unknown", "0_not_virus"]
print("=== по числу ===")
for c in cats:
    vs = [r for r in out if r["category"] == c]
    print("%-15s %3d  %s" % (c, len(vs), " ".join(r["votu"][5:] for r in vs)))
print("=== подкатегории ===")
sc = defaultdict(list)
for r in out:
    if r["detail"]:
        sc[(r["category"], r["detail"])].append(r["votu"][5:])
for k in sorted(sc):
    print(k, len(sc[k]), " ".join(sc[k]))
print("=== обнаружено в образце, по категориям ===")
for m in MV:
    line = []
    for c in cats:
        line.append("%s=%d" % (c.split("_")[0], sum(1 for r in out if r["category"] == c and det[(m, r["votu"])])))
    print(NAME[m], sum(det[(m, v)] for v in VOTUS), " ".join(line))
print("=== доля RPM каталога по категориям, % ===")
for m in MV:
    tot = sum(rpm[(m, v)] for v in VOTUS)
    parts = []
    for c in cats:
        s = sum(rpm[(m, r["votu"])] for r in out if r["category"] == c)
        parts.append("%s=%.1f" % (c.split("_")[0], 100 * s / tot))
    print(NAME[m], "total_RPM=%.0f" % tot, " ".join(parts))
print("=== таксоны по категориям ===")
for c in cats:
    t = defaultdict(int)
    for r in out:
        if r["category"] == c:
            t[r["taxon"]] += 1
    print(c, dict(t))
print("=== CheckV по категориям ===")
for c in cats:
    t = defaultdict(int)
    for r in out:
        if r["category"] == c:
            t[r["checkv_quality"]] += 1
    print(c, dict(t))
print("=== длина, кб: медиана по категориям ===")
from statistics import median
for c in cats:
    L = [r["length"] for r in out if r["category"] == c]
    print(c, "%.1f" % (median(L) / 1000), "min %.1f max %.1f" % (min(L) / 1000, max(L) / 1000))
print("=== топ-10 по RPM в каждом образце ===")
for m in MV:
    top = sorted(VOTUS, key=lambda v: -rpm[(m, v)])[:10]
    cm = {r["votu"]: r["category"] for r in out}
    print(NAME[m], ", ".join("%s(%.0f,%s)" % (v[5:], rpm[(m, v)], cm[v].split("_")[0]) for v in top))
