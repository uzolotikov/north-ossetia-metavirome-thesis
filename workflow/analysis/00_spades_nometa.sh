#!/bin/bash
# Сборка трёх фагомов SPAdes без --meta для сравнения с metaSPAdes (таблица — 00_asm_stats.py).
source "$(dirname "$0")/../common.sh"
set +e   # сбой одной сборки не останавливает остальные
E=$ENVS/metaspades-env/bin
C=$OUT/clean/pilot
O=$ANALYSIS/assembly_spades
mkdir -p $O
for s in MV-4128_S29 MV-4125_S18 MV-4122_S12; do
  [ -s $O/$s/contigs.fasta ] && continue
  $E/spades.py -1 $C/${s}_R1.fastq.gz -2 $C/${s}_R2.fastq.gz -o $O/$s -t $THREADS -m 150 > $O/$s.log 2>&1 || echo "FAIL $s"
done
echo SPADES_DONE
