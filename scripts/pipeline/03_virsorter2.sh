#!/bin/bash
# Шаг 3. Поиск вирусных контигов (VirSorter2 2.2.4). Конфигурация запуска — virsorter2_config.yaml.
set -euo pipefail
source "$(dirname "$0")/config.sh"; cd "$P"
VS=$ENVS/virsorter2-env/bin/virsorter
# база VirSorter2, устанавливается один раз (установлена 22.09.2026)
[ -e databases/virsorter2/Done_all_setup ] || $VS setup -d $P/databases/virsorter2 -j 4 --conda-frontend conda
for s in $MV; do
  $VS run -i assembly/no_NC/$s/contigs.fasta -w viral/no_NC/$s -d $P/databases/virsorter2 \
    --include-groups dsDNAphage,NCLDV,ssDNA,lavidaviridae --min-length 1500 --min-score 0.5 \
    --keep-original-seq -j 4 all --conda-frontend conda > logs/$s.virsorter2.log 2>&1
done
