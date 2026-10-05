#!/bin/bash
set -uo pipefail
E=/NatureUsers/nsinyavskiy/miniconda3/envs/metaspades-env/bin
C=/NatureUsers/nsinyavskiy/metaviroms_kgt/clean/pilot
O=/NatureUsers/glitch/mv_review/assembly_spades
mkdir -p $O
for s in MV-4128_S29 MV-4125_S18 MV-4122_S12; do
  [ -s $O/$s/contigs.fasta ] && continue
  $E/spades.py -1 $C/${s}_R1.fastq.gz -2 $C/${s}_R2.fastq.gz -o $O/$s -t 32 -m 150 > $O/$s.log 2>&1 || echo "FAIL $s"
done
echo SPADES_DONE
