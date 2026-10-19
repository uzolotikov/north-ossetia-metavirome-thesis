#!/usr/bin/env python3
"""Характеристики сборок: фагомы (metaSPAdes и SPAdes без --meta), клеточные фракции Тиб
и бактериальная сборка Бирагзанга.

Длины берутся по самим последовательностям, а не по именам контигов. Пути в столбце file
даны относительно рабочего каталога OUT; для сборки Бирагзанга — имя файла.

Запуск: 00_asm_stats.py > analysis/report/assembly_stats.tsv
"""
import os
from paths import OUT, ANALYSIS, BIRAGZANG

MV = ["MV-4122_S12", "MV-4125_S18", "MV-4128_S29"]
TIB_CELLULAR = ["4122-0-22_S11", "4122-0-44_S10", "4125-0-22_S17", "4125-0-44_S16"]


def lengths(path):
    out, n = [], None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if n is not None:
                    out.append(n)
                n = 0
            else:
                n += len(line.strip())
    if n is not None:
        out.append(n)
    return out


def stats(path):
    L = sorted(lengths(path), reverse=True)
    total, c, n50 = sum(L), 0, 0
    for x in L:
        c += x
        if c >= total / 2:
            n50 = x
            break
    return len(L), total, L[0], n50, sum(x >= 5000 for x in L), sum(x >= 10000 for x in L)


rows = []
for s in MV:
    rows.append((s, "metaSPAdes 4.3.0", f"{OUT}/assembly/no_NC/{s}/contigs.fasta"))
    rows.append((s, "SPAdes (no --meta)", f"{ANALYSIS}/assembly_spades/{s}/contigs.fasta"))
for s in TIB_CELLULAR:
    rows.append((s, "metaSPAdes 4.3.0", f"{OUT}/assembly/cellular/{s}/contigs.fasta"))
rows.append(("Biragzang.final", "external", BIRAGZANG))

print("assembly\tassembler\tcontigs\ttotal_bp\tlargest_bp\tN50\tge5kb\tge10kb\tfile")
for name, assembler, path in rows:
    shown = os.path.relpath(path, OUT) if path.startswith(OUT + "/") else os.path.basename(path)
    print("\t".join([name, assembler] + [str(x) for x in stats(path)] + [shown]))
