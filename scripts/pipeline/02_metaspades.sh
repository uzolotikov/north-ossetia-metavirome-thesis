#!/bin/bash
# Шаг 2. Сборка (metaSPAdes 4.3.0): три фагома и четыре клеточные фракции Тиб.
set -euo pipefail
source "$(dirname "$0")/config.sh"; cd "$P"
for s in $MV; do
  $ENVS/metaspades-env/bin/metaspades.py -1 clean/pilot/${s}_R1.fastq.gz -2 clean/pilot/${s}_R2.fastq.gz \
    -o assembly/no_NC/$s -t $THREADS -m 32 > logs/$s.metaspades.log 2>&1
done
for s in $TIB_CELLULAR; do
  $ENVS/metaspades-env/bin/metaspades.py -1 clean/companions/${s}_R1.fastq.gz -2 clean/companions/${s}_R2.fastq.gz \
    -o assembly/cellular/$s -t $THREADS -m 32 > logs/$s.cellular.metaspades.log 2>&1
done
