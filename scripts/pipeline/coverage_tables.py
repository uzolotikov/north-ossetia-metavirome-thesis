#!/usr/bin/env python3
"""Таблицы покрытия каталога vOTU для одной библиотеки по отфильтрованному BAM.

  <lib>.idxstats.tsv  samtools idxstats: vOTU, длина, картированные и некартированные прочтения
  <lib>.coverage.tsv  sample, votu, length, mapped_reads, mean_depth, breadth_1x_pct, breadth_5x_pct
                      глубина — samtools depth -a -s -q 20: все позиции, перекрывающиеся мейты пары
                      учитываются один раз, только основания с качеством >=20

Использование: coverage_tables.py <samtools> <filtered.bam> <sample> <out_dir>
"""
import subprocess, sys
from collections import defaultdict

st, bam, sample, out = sys.argv[1:5]
idx = subprocess.check_output([st, "idxstats", bam], universal_newlines=True)
with open("%s/%s.idxstats.tsv" % (out, sample), "w") as f:
    f.write(idx)
length, mapped = {}, {}
for line in idx.splitlines():
    c, L, m, _ = line.split("\t")
    if c != "*":
        length[c] = int(L); mapped[c] = int(m)
dsum, b1, b5 = defaultdict(int), defaultdict(int), defaultdict(int)
p = subprocess.Popen([st, "depth", "-a", "-s", "-q", "20", bam], stdout=subprocess.PIPE, universal_newlines=True)
for line in p.stdout:
    c, _, d = line.rstrip("\n").split("\t"); d = int(d)
    dsum[c] += d; b1[c] += d >= 1; b5[c] += d >= 5
p.wait()
with open("%s/%s.coverage.tsv" % (out, sample), "w") as f:
    f.write("sample\tvotu\tlength\tmapped_reads\tmean_depth\tbreadth_1x_pct\tbreadth_5x_pct\n")
    for c in length:
        L = length[c]
        f.write("%s\t%s\t%d\t%d\t%.4f\t%.2f\t%.2f\n" % (sample, c, L, mapped[c], dsum[c] / L, 100 * b1[c] / L, 100 * b5[c] / L))
