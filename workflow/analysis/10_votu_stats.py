#!/usr/bin/env python3
"""Сравнительная статистика каталога vOTU: отношение обилий контаминанта между фагомами,
vOTU, общие для Тиб-2 и Бирагзанга, подтверждение vOTU клеточной фракцией своей точки
и варианты матрицы детекции.

Выход (в report/):
  votu_contaminant_ratio.tsv   отношение обилий контаминанта между библиотеками
  votu_shared_ratios.tsv       24 vOTU, общие для Тиб-2 и Бирагзанга
  votu_confirmed_cellular.tsv  vOTU, подтверждённые клеточной фракцией своей точки
  votu_matrix_versions.tsv     варианты матрицы детекции, индексы Жаккара и Брея — Кёртиса
"""
import csv
from collections import defaultdict
from paths import OUT as K, ANALYSIS as W, REPORT as REP

T1, T2, BI = "MV-4122_S12", "MV-4125_S18", "MV-4128_S29"
MV = [T1, T2, BI]
SHORT = {T1: "Tib1", T2: "Tib2", BI: "Biragzang"}
CELLS = {T1: ["4122-0-44_S10", "4122-0-22_S11"],
         T2: ["4125-0-44_S16", "4125-0-22_S17"],
         BI: ["4128-0-44_S27", "4128-0-22_S28", "4128_S30"]}

# ---------------------------------------------------------------- входы
rpm, reads, det, length = {}, {}, {}, {}
with open(f"{K}/qc/votu97_comparison/votu97_all_libraries.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        k = (r["sample"], r["votu"])
        rpm[k] = float(r["RPM"])
        reads[k] = int(r["mapped_reads_filtered"])
        det[k] = int(r["detected_75pct_1x"])
        length[r["votu"]] = int(r["length_bp"])
VOTUS = sorted(length)

# встроенные профаги: контиг-хозяин каждой vOTU и отношение глубины «внутри / фланги»
embedded = {}
with open(f"{W}/prophage/embedded.tsv") as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2:
            embedded[p[0]] = p[1]
ratios = defaultdict(dict)
with open(f"{W}/prophage/prophage_ratios.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        try:
            ratios[r["votu"]][r["library"]] = float(r["ratio_in_over_flank"])
        except ValueError:
            pass
dormant, induced = [], []
for v, host in embedded.items():
    lib = host.split("__")[0]
    x = ratios.get(v, {}).get(lib)
    if x is None:
        continue
    (induced if x >= 1.5 else dormant).append(v)
dormant, induced = sorted(dormant), sorted(induced)

# ---------------------------------------------------------------- отношение обилий контаминанта
# Отношение вычисляется, а не берётся готовым. Важны два решения, оба записаны здесь.
#
# 1. Набор контигов. Рассуждение об общих vOTU касается контаминанта, преобладающего в MV-4128,
#    поэтому основа — контиги MV-4128 длиной >=50 кб. Если объединить контиги >=50 кб всех трёх
#    фагомов, примешивается собственное сообщество MV-4122, и отношение Тиб-1 / Бирагзанг теряет
#    смысл (0,64 вместо ~0,03).
#
# 2. Число прочтений, а не глубина. RPM vOTU в таблице покрытия считается по прочтениям, поэтому
#    и отношение для контаминанта должно считаться так же. Средняя глубина, взвешенная по длине,
#    из step6_contig_classes.tsv даёт для тех же контигов Тиб-2 / Бирагзанг = 0,34, потому что
#    meandepth samtools coverage и число прочтений idxstats на этом наборе расходятся. Смешение
#    отношения по глубине с отношениями vOTU по прочтениям сдвинуло бы число общих vOTU; отношение
#    по прочтениям воспроизводит доли прочтений контаминанта (68,6 % / 37,3 % прочтений MV-4128 /
#    MV-4125) и используется для окна «в пределах 1,5 раза».
cls = {}
with open(f"{W}/report/step6_contig_classes.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        cls[r["contig"]] = (int(r["length"]), r["origin"], r["class"])
targ = {c for c, (L, o, k) in cls.items() if c.startswith("MV-4128_S29__") and L >= 50000}


def mapped_reads(lib, contigs):
    tot = 0
    with open(f"{W}/step6/cov/{lib}.idxstats.tsv") as fh:
        for line in fh:
            p_ = line.rstrip().split("\t")
            if len(p_) >= 3 and p_[0] in contigs:
                tot += int(p_[2])
    return tot


clean = {}
with open(f"{K}/qc/votu97_comparison/votu97_all_libraries.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        clean[r["sample"]] = int(r["clean_reads"])

cont_rpm, cont_reads = {}, {}
for m in MV:
    cont_reads[m] = mapped_reads(m, targ)
    cont_rpm[m] = 1e6 * cont_reads[m] / clean[m]

# вариант по глубине, взвешенной по длине, — для сравнения
dep = defaultdict(float)
with open(f"{W}/report/step6_contig_classes.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["contig"] in targ:
            L = int(r["length"])
            for m in MV:
                dep[m] += float(r[f"depth:{m}"]) * L

with open(f"{REP}/votu_contaminant_ratio.tsv", "w") as out:
    out.write("# contaminant abundance ratio between the MV libraries\n")
    out.write(f"# basis: {len(targ)} MV-4128-derived contigs >=50 kb present in the step6 reference,\n")
    out.write(f"#        {sum(cls[c][0] for c in targ)} bp\n")
    out.write("library\tmapped_reads\tclean_reads\tread_pct_of_library\tRPM\n")
    for m in MV:
        out.write(f"{SHORT[m]}\t{cont_reads[m]}\t{clean[m]}\t"
                  f"{100.0*cont_reads[m]/clean[m]:.2f}\t{cont_rpm[m]:.0f}\n")
    out.write("#\npair\tratio_by_reads\tratio_by_length_weighted_depth\n")
    for x, y in [(T2, BI), (T1, BI), (T1, T2)]:
        out.write(f"{SHORT[x]}/{SHORT[y]}\t{cont_rpm[x]/cont_rpm[y]:.4f}\t{dep[x]/dep[y]:.4f}\n")
    out.write("#\n# the read-based Tib2/Biragzang ratio is used for the 1.5x window of shared vOTU\n")
R_CONT = cont_rpm[T2] / cont_rpm[BI]

# ---------------------------------------------------------------- общие vOTU
shared = [v for v in VOTUS if det[(T2, v)] == 1 and det[(BI, v)] == 1 and rpm[(BI, v)] > 0]
band_lo, band_hi = R_CONT / 1.5, R_CONT * 1.5
n_band = 0
with open(f"{REP}/votu_shared_ratios.tsv", "w") as out:
    out.write("# vOTU detected in BOTH MV Tib2 and MV Biragzang\n")
    out.write(f"# contaminant Tib2/Biragzang ratio (derived) = {R_CONT:.4f}\n")
    out.write(f"# 'within 1.5x' window = {band_lo:.4f} .. {band_hi:.4f}\n")
    out.write("votu\trpm_Tib2\trpm_Biragzang\tratio\twithin_1.5x_of_contaminant\n")
    for v in sorted(shared, key=lambda v: rpm[(T2, v)] / rpm[(BI, v)]):
        r = rpm[(T2, v)] / rpm[(BI, v)]
        ok = band_lo <= r <= band_hi
        n_band += ok
        out.write(f"{v}\t{rpm[(T2,v)]:.2f}\t{rpm[(BI,v)]:.2f}\t{r:.4f}\t{'yes' if ok else 'no'}\n")
    out.write(f"#\n# shared vOTU                    : {len(shared)}\n")
    out.write(f"# within 1.5x of the contaminant : {n_band}\n")

# ---------------------------------------------------------------- подтверждение клеточной фракцией
LEAK = 5e-4          # фон: доля MV-специфичной бактериальной ДНК, просачивающейся в клеточные библиотеки
CONFIRM = 40 * LEAK  # 0,02 — подтверждена
STRONG = 0.1         # подтверждена надёжно
conf = defaultdict(list)
with open(f"{REP}/votu_confirmed_cellular.tsv", "w") as out:
    out.write("# vOTU confirmed by the cellular fraction of their own site\n")
    out.write(f"# leak background {LEAK:g}; confirmed if cell/MV RPM ratio > {CONFIRM:g} with >=10 reads;\n")
    out.write(f"# reliable if ratio >= {STRONG:g}\n")
    out.write("site\tvotu\tcell_library\tcell_reads\tratio_cell_over_mv\tlevel"
              "\trpm_Tib1\trpm_Tib2\trpm_Biragzang\tmax_in\n")
    for m in MV:
        for v in VOTUS:
            if rpm[(m, v)] <= 0:
                continue
            for lib in CELLS[m]:
                if reads[(lib, v)] < 10:
                    continue
                ratio = rpm[(lib, v)] / rpm[(m, v)]
                if ratio <= CONFIRM:
                    continue
                lvl = "reliable" if ratio >= STRONG else "confirmed"
                conf[m].append((v, lvl))
                vals = [rpm[(x, v)] for x in MV]
                mx = SHORT[MV[vals.index(max(vals))]]
                out.write(f"{SHORT[m]}\t{v}\t{lib}\t{reads[(lib,v)]}\t{ratio:.4f}\t{lvl}\t"
                          f"{vals[0]:.1f}\t{vals[1]:.1f}\t{vals[2]:.1f}\t{mx}\n")
                break
    out.write("#\n")
    for m in MV:
        tot = len(conf[m])
        rel = sum(1 for _, l in conf[m] if l == "reliable")
        out.write(f"# {SHORT[m]}: confirmed {tot}, of which reliable {rel}\n")
    rel_t2 = [v for v, l in conf[T2] if l == "reliable"]
    peak_t1 = sum(1 for v in rel_t2
                  if max(MV, key=lambda m: rpm[(m, v)]) == T1)
    out.write(f"# of the {len(rel_t2)} reliable Tib2 vOTU, {peak_t1} peak in MV Tib1\n")

# ---------------------------------------------------------------- варианты матрицы детекции
def jaccard(a, b):
    A, B = set(a), set(b)
    return 100.0 * len(A & B) / len(A | B) if A | B else float("nan")


def bray(keep):
    out = {}
    for a, b in [(T1, T2), (T1, BI), (T2, BI)]:
        xa = [rpm[(a, v)] for v in keep]
        xb = [rpm[(b, v)] for v in keep]
        s = sum(xa) + sum(xb)
        out[f"{SHORT[a]}-{SHORT[b]}"] = (2 * sum(min(p, q) for p, q in zip(xa, xb)) / s) if s else float("nan")
    return out


follows = [v for v in shared if band_lo <= rpm[(T2, v)] / rpm[(BI, v)] <= band_hi]
versions = [
    ("(a) all 97 vOTU", VOTUS),
    ("(v1) minus dormant proviruses", [v for v in VOTUS if v not in dormant]),
    ("(v2) v1 minus vOTU following the contaminant",
     [v for v in VOTUS if v not in dormant and v not in follows]),
]
with open(f"{REP}/votu_matrix_versions.tsv", "w") as out:
    out.write("# detection matrix versions\n")
    out.write(f"# dormant proviruses excluded in v1/v2 ({len(dormant)}): {','.join(dormant)}\n")
    out.write(f"# induced prophages (kept, {len(induced)}): {','.join(induced)}\n")
    out.write(f"# follows the contaminant, excluded in v2 ({len(follows)}): {','.join(follows)}\n")
    out.write("version\tn_votu\tdet_Tib1\tdet_Tib2\tdet_Biragzang\tshared_all3"
              "\tJ_Tib1_Tib2\tJ_Tib1_Bi\tJ_Tib2_Bi\tBC_Tib1_Tib2\tBC_Tib1_Bi\tBC_Tib2_Bi\n")
    for name, keep in versions:
        d = {m: [v for v in keep if det[(m, v)] == 1] for m in MV}
        all3 = set(d[T1]) & set(d[T2]) & set(d[BI])
        bc = bray(keep)
        out.write(f"{name}\t{len(keep)}\t{len(d[T1])}\t{len(d[T2])}\t{len(d[BI])}\t{len(all3)}\t"
                  f"{jaccard(d[T1],d[T2]):.1f}\t{jaccard(d[T1],d[BI]):.1f}\t{jaccard(d[T2],d[BI]):.1f}\t"
                  f"{bc['Tib1-Tib2']:.2f}\t{bc['Tib1-Biragzang']:.2f}\t{bc['Tib2-Biragzang']:.2f}\n")

for f in ("votu_contaminant_ratio", "votu_shared_ratios",
          "votu_confirmed_cellular", "votu_matrix_versions"):
    print(f"=== {f}.tsv ===")
    print(open(f"{REP}/{f}.tsv").read().rstrip())
    print()
print("VOTU_STATS_DONE")
