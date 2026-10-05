#!/bin/bash
# Шаг 1. Тримминг прочтений (fastp 0.23.2). Одинаковые параметры для всех 11 библиотек.
set -euo pipefail
source "$(dirname "$0")/config.sh"; cd "$P"
run(){ local lib=$1 grp=$2
  mkdir -p clean/$grp qc/fastp_$grp
  $CONDA/bin/fastp --in1 $RAW/${lib}_R1_001.fastq.gz --in2 $RAW/${lib}_R2_001.fastq.gz \
    --out1 clean/$grp/${lib}_R1.fastq.gz --out2 clean/$grp/${lib}_R2.fastq.gz \
    --detect_adapter_for_pe --trim_poly_g --poly_g_min_len 10 \
    --cut_tail --cut_tail_window_size 4 --cut_tail_mean_quality 20 \
    --qualified_quality_phred 20 --unqualified_percent_limit 40 --n_base_limit 5 \
    --length_required 75 --thread 2 \
    --html qc/fastp_$grp/$lib.html --json qc/fastp_$grp/$lib.json --dont_overwrite
}
for lib in $PILOT; do run $lib pilot; done
for lib in $COMPANIONS; do run $lib companions; done
