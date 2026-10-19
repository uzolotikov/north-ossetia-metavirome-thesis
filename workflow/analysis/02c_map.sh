#!/bin/bash
# Картирование библиотек на общий референс (02a_build_ref.sh): bowtie2, первичные выравнивания,
# не более 5 % несовпадений; idxstats и samtools coverage без порога MAPQ и с MAPQ >=20.
# Использование: 02c_map.sh <библиотека> [<библиотека> ...]; без аргументов — все 11 библиотек.
source "$(dirname "$0")/../common.sh"
set +e   # сбой одной библиотеки не останавливает остальные
W=$ANALYSIS/step6
M=$ENVS/mapping-env/bin
mkdir -p $W/bam $W/cov
for lib in ${@:-$ALL_LIBS}; do
  if [[ $lib == MV-* || $lib == NC1* ]]; then d=$OUT/clean/pilot; else d=$OUT/clean/companions; fi
  [ -s $W/cov/$lib.cov_mq0.tsv ] && continue
  $M/bowtie2 -p $THREADS -x $W/ref_derep -1 $d/${lib}_R1.fastq.gz -2 $d/${lib}_R2.fastq.gz --no-unal 2> $W/bam/$lib.bowtie2.log \
   | $M/samtools view -u -F 0x904 -e "[NM]<=0.05*qlen" - | $M/samtools sort -@ 8 -m 2G -o $W/bam/$lib.bam - \
   && $M/samtools index $W/bam/$lib.bam \
   && $M/samtools idxstats $W/bam/$lib.bam > $W/cov/$lib.idxstats.tsv \
   && $M/samtools coverage $W/bam/$lib.bam > $W/cov/$lib.cov_mq0.tsv \
   && $M/samtools coverage -q 20 $W/bam/$lib.bam > $W/cov/$lib.cov_mq20.tsv \
   && echo "MAPPED $lib" || echo "MAPFAIL $lib"
done
echo MAP_BATCH_DONE
