#!/usr/bin/env python3
"""Flagged catalogue of all 97 vOTU.

Flags each vOTU by: geNomad verdict, embedded-prophage status, whether its
abundance follows the contaminant, and whether the cellular fraction of its own
site confirms it. Produces a keep / flag / exclude recommendation for every vOTU.

Output: report/votu97_catalog_flags.tsv
"""
import csv
import os
from collections import defaultdict

K = "/NatureUsers/nsinyavskiy/metaviroms_kgt"
W = "/NatureUsers/glitch/mv_review"
REP = f"{W}/report"
G = f"{W}/genomad/votu97/representatives_ge5kb_summary"

T1, T2, BI = "MV-4122_S12", "MV-4125_S18", "MV-4128_S29"
MV = [T1, T2, BI]
SHORT = {T1: "Tib1", T2: "Tib2", BI: "Biragzang"}
CELLS = {T1: ["4122-0-44_S10", "4122-0-22_S11"],
         T2: ["4125-0-44_S16", "4125-0-22_S17"],
         BI: ["4128-0-44_S27", "4128-0-22_S28", "4128_S30"]}

rpm, reads, det, length, clean = {}, {}, {}, {}, {}
with open(f"{K}/qc/votu97_comparison/votu97_all_libraries.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        k = (r["sample"], r["votu"])
        rpm[k] = float(r["RPM"]); reads[k] = int(r["mapped_reads_filtered"])
        det[k] = int(r["detected_75pct_1x"]); length[r["votu"]] = int(r["length_bp"])
        clean[r["sample"]] = int(r["clean_reads"])
VOTUS = sorted(length)

# geNomad verdicts
gv, gscore, gtax = {}, {}, {}
with open(f"{G}/representatives_ge5kb_virus_summary.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        v = r["seq_name"].split("|")[0]
        gv[v] = "virus"; gscore[v] = float(r["virus_score"])
        tax = [t for t in r["taxonomy"].split(";") if t]
        gtax[v] = tax[-1] if tax else "NA"
with open(f"{G}/representatives_ge5kb_plasmid_summary.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        v = r["seq_name"].split("|")[0]
        gv[v] = "plasmid"; gscore[v] = float(r["plasmid_score"])
        gtax[v] = "NA"

# embedded prophages
host = {}
with open(f"{W}/prophage/embedded.tsv") as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2:
            host[p[0]] = p[1]
ratios = defaultdict(dict)
with open(f"{W}/prophage/prophage_ratios.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        try:
            ratios[r["votu"]][r["library"]] = float(r["ratio_in_over_flank"])
        except ValueError:
            pass
emb = {}
for v, h in host.items():
    lib = h.split("__")[0]
    x = ratios.get(v, {}).get(lib)
    if x is not None:
        emb[v] = ("induced_prophage" if x >= 1.5 else "dormant_provirus", x, h)

# contaminant ratio, read-based (same basis as 10_votu_stats.py)
cls = {}
with open(f"{REP}/step6_contig_classes.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        cls[r["contig"]] = int(r["length"])
targ = {c for c, L in cls.items() if c.startswith("MV-4128_S29__") and L >= 50000}
cr = {}
for m in MV:
    tot = 0
    with open(f"{W}/step6/cov/{m}.idxstats.tsv") as fh:
        for line in fh:
            p = line.rstrip().split("\t")
            if len(p) >= 3 and p[0] in targ:
                tot += int(p[2])
    cr[m] = 1e6 * tot / clean[m]
R = cr[T2] / cr[BI]
lo, hi = R / 1.5, R * 1.5

# cellular confirmation
LEAK, CONFIRM, STRONG = 5e-4, 0.02, 0.1
conf = {}
for m in MV:
    for v in VOTUS:
        if rpm[(m, v)] <= 0:
            continue
        for lib in CELLS[m]:
            if reads[(lib, v)] < 10:
                continue
            ratio = rpm[(lib, v)] / rpm[(m, v)]
            if ratio > CONFIRM:
                conf[v] = ("reliable" if ratio >= STRONG else "confirmed", SHORT[m], lib, ratio)
            break

rows = []
for v in VOTUS:
    verdict = gv.get(v, "not_viral")
    e = emb.get(v)
    follows = ""
    if det[(T2, v)] == 1 and det[(BI, v)] == 1 and rpm[(BI, v)] > 0:
        r = rpm[(T2, v)] / rpm[(BI, v)]
        follows = "yes" if lo <= r <= hi else "no"
    c = conf.get(v)
    # recommendation
    if verdict == "plasmid":
        rec = "exclude_plasmid"
    elif verdict == "not_viral":
        rec = "exclude_not_confirmed_by_genomad"
    elif e and e[0] == "dormant_provirus":
        rec = "flag_dormant_provirus"
    elif c and c[0] == "reliable":
        rec = "keep_confirmed_by_cellular"
    else:
        rec = "keep"
    rows.append((v, length[v], verdict, gscore.get(v, float("nan")), gtax.get(v, "NA"),
                 e[0] if e else "", f"{e[1]:.2f}" if e else "", e[2] if e else "",
                 follows, c[0] if c else "", c[1] if c else "", f"{c[3]:.3f}" if c else "",
                 rpm[(T1, v)], rpm[(T2, v)], rpm[(BI, v)], rec))

with open(f"{REP}/votu97_catalog_flags.tsv", "w") as out:
    out.write("# flags for all 97 vOTU of the catalogue\n")
    out.write(f"# contaminant Tib2/Biragzang ratio = {R:.4f}; 'follows_contaminant' window {lo:.4f}..{hi:.4f}\n")
    out.write("# recommendation: exclude_plasmid | exclude_not_confirmed_by_genomad |\n")
    out.write("#                 flag_dormant_provirus | keep_confirmed_by_cellular | keep\n")
    out.write("votu\tlength\tgenomad\tgenomad_score\tgenomad_taxon\tembedded\t"
              "in_over_flank\thost_contig\tfollows_contaminant\tcellular\tcellular_site\t"
              "cellular_ratio\trpm_Tib1\trpm_Tib2\trpm_Biragzang\trecommendation\n")
    for r in rows:
        out.write("\t".join(str(x) if not isinstance(x, float) else f"{x:.2f}" for x in r) + "\n")
    out.write("#\n")
    cnt = defaultdict(int)
    for r in rows:
        cnt[r[-1]] += 1
    for k in sorted(cnt):
        out.write(f"# {k}: {cnt[k]}\n")
    out.write(f"# total: {len(rows)}\n")

print(open(f"{REP}/votu97_catalog_flags.tsv").read().rstrip().split("#\n")[-1])
print("CATALOG_FLAGS_DONE")
