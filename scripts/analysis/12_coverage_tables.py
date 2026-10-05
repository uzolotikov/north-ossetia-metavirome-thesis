#!/usr/bin/env python3
"""Coverage of the 97 vOTU across all 11 libraries.

Outputs:
  report/votu_coverage_by_library.tsv  per-library summary
  report/votu_coverage_matrix.tsv      full 97 x 11 matrix (breadth % / mean depth / RPM)
"""
import csv
from statistics import median

K = "/NatureUsers/nsinyavskiy/metaviroms_kgt"
W = "/NatureUsers/glitch/mv_review"
REP = f"{W}/report"

ORDER = ["MV-4122_S12", "MV-4125_S18", "MV-4128_S29",
         "4122-0-44_S10", "4122-0-22_S11",
         "4125-0-44_S16", "4125-0-22_S17",
         "4128-0-44_S27", "4128-0-22_S28", "4128_S30", "NC1_S64"]
LABEL = {"MV-4122_S12": "MV Tib1", "MV-4125_S18": "MV Tib2", "MV-4128_S29": "MV Biragzang",
         "4122-0-44_S10": "Tib1 cell 0.44", "4122-0-22_S11": "Tib1 cell 0.22",
         "4125-0-44_S16": "Tib2 cell 0.44", "4125-0-22_S17": "Tib2 cell 0.22",
         "4128-0-44_S27": "Bi cell 0.44", "4128-0-22_S28": "Bi cell 0.22",
         "4128_S30": "Bi whole water", "NC1_S64": "NC1 control"}

br, dp, rpm, rd, clean, length = {}, {}, {}, {}, {}, {}
with open(f"{K}/qc/votu97_comparison/votu97_all_libraries.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        k = (r["sample"], r["votu"])
        br[k] = float(r["breadth_1x_pct"]); dp[k] = float(r["mean_depth"])
        rpm[k] = float(r["RPM"]); rd[k] = int(r["mapped_reads_filtered"])
        clean[r["sample"]] = int(r["clean_reads"]); length[r["votu"]] = int(r["length_bp"])
VOTUS = sorted(length)

# share of library reads falling on whole-viral contigs of the common reference
viral_frac = {}
with open(f"{REP}/step6_read_fractions.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        tot = 0.0
        for k, v in r.items():
            if k.startswith("class:viral|") and not k.endswith("unmapped"):
                try:
                    tot += float(v)
                except (TypeError, ValueError):
                    pass
        viral_frac[r["library"]] = tot

with open(f"{REP}/votu_coverage_by_library.tsv", "w") as out:
    out.write("# coverage of the 97-vOTU catalogue per library\n")
    out.write("# detected = breadth >=75 % at >=1x (detection criterion of the catalogue)\n")
    out.write("# catalogue_read_pct = reads on the catalogue / clean reads\n")
    out.write("# viral_read_pct = share of library reads on whole-viral contigs of the step6 reference\n")
    out.write("library\tlabel\tclean_reads\tdetected_75pct\tbreadth_ge90\tbreadth_ge50\t"
              "zero_reads\tcatalogue_reads\tcatalogue_read_pct\tmedian_breadth_detected\t"
              "median_depth_detected\tmax_depth\tmax_depth_votu\tviral_read_pct\n")
    for s in ORDER:
        det = [v for v in VOTUS if br[(s, v)] >= 75]
        d90 = sum(1 for v in VOTUS if br[(s, v)] >= 90)
        d50 = sum(1 for v in VOTUS if br[(s, v)] >= 50)
        zero = sum(1 for v in VOTUS if rd[(s, v)] == 0)
        creads = sum(rd[(s, v)] for v in VOTUS)
        mb = median([br[(s, v)] for v in det]) if det else 0.0
        md = median([dp[(s, v)] for v in det]) if det else 0.0
        mx = max(VOTUS, key=lambda v: dp[(s, v)])
        out.write(f"{s}\t{LABEL[s]}\t{clean[s]}\t{len(det)}\t{d90}\t{d50}\t{zero}\t{creads}\t"
                  f"{100.0*creads/clean[s]:.3f}\t{mb:.1f}\t{md:.2f}\t{dp[(s,mx)]:.1f}\t{mx}\t"
                  f"{viral_frac.get(s, float('nan')):.2f}\n")

with open(f"{REP}/votu_coverage_matrix.tsv", "w") as out:
    out.write("# full coverage matrix: breadth_1x_pct / mean_depth / RPM per vOTU per library\n")
    out.write("votu\tlength")
    for s in ORDER:
        out.write(f"\tbreadth:{s}\tdepth:{s}\tRPM:{s}")
    out.write("\n")
    for v in VOTUS:
        out.write(f"{v}\t{length[v]}")
        for s in ORDER:
            out.write(f"\t{br[(s,v)]:.1f}\t{dp[(s,v)]:.2f}\t{rpm[(s,v)]:.2f}")
        out.write("\n")

print(open(f"{REP}/votu_coverage_by_library.tsv").read().rstrip())
print("\nCOVERAGE_TABLES_DONE")
