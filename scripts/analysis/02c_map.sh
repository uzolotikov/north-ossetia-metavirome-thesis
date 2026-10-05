#!/bin/bash
# usage: 02c_map.sh THREADS LIB [LIB...]
set -uo pipefail
T=$1; shift
K=/NatureUsers/nsinyavskiy/metaviroms_kgt; W=/NatureUsers/glitch/mv_review/step6
M=/NatureUsers/nsinyavskiy/miniconda3/envs/mapping-env/bin
mkdir -p $W/bam $W/cov
for lib in "$@"; do
  if [[ $lib == MV-* || $lib == NC1* ]]; then d=$K/clean/pilot; else d=$K/clean/companions; fi
  [ -s $W/cov/$lib.cov_mq0.tsv ] && continue
  $M/bowtie2 -p $T -x $W/ref_derep -1 $d/${lib}_R1.fastq.gz -2 $d/${lib}_R2.fastq.gz --no-unal 2> $W/bam/$lib.bowtie2.log \
   | $M/samtools view -u -F 0x904 -e "[NM]<=0.05*qlen" - | $M/samtools sort -@ 8 -m 2G -o $W/bam/$lib.bam - \
   && $M/samtools index $W/bam/$lib.bam \
   && $M/samtools idxstats $W/bam/$lib.bam > $W/cov/$lib.idxstats.tsv \
   && $M/samtools coverage $W/bam/$lib.bam > $W/cov/$lib.cov_mq0.tsv \
   && $M/samtools coverage -q 20 $W/bam/$lib.bam > $W/cov/$lib.cov_mq20.tsv \
   && echo "MAPPED $lib" || echo "MAPFAIL $lib"
done
echo MAP_BATCH_DONE
