#!/bin/bash
# Шаг 8. Покрытие каталога vOTU во всех 11 библиотеках.
# 8а: по отфильтрованному BAM каждой библиотеки — idxstats и покрытие (samtools depth -a -s -q 20);
# 8б: сводная таблица — число прочтений, RPM, RPKM, средняя глубина, ширина, детекция.
set -euo pipefail
source "$(dirname "$0")/../common.sh"; cd "$OUT"
SAM=$ENVS/mapping-env/bin/samtools
for lib in $PILOT;      do mkdir -p mapping/no_NC/coverage;      $PY $PIPE/coverage_tables.py $SAM mapping/no_NC/filtered/$lib.filtered.bam $lib mapping/no_NC/coverage; done
for lib in $COMPANIONS; do mkdir -p mapping/companions/coverage; $PY $PIPE/coverage_tables.py $SAM mapping/companions/filtered/$lib.filtered.bam $lib mapping/companions/coverage; done
$PY $PIPE/rebuild_votu97.py "$OUT" qc/votu97_comparison
