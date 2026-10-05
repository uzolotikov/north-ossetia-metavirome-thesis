#!/bin/bash
# Шаг 7. Картирование всех 11 библиотек на каталог vOTU (bowtie2 2.5.5, samtools 1.24).
# Фильтр: MAPQ >= 20, только первичные картированные прочтения, не более 5 % несовпадений.
set -euo pipefail
source "$(dirname "$0")/config.sh"; cd "$P"
M=$ENVS/mapping-env/bin; IDX=mapping/no_NC/index/votu97
mkdir -p mapping/no_NC/index
$M/bowtie2-build votu/no_NC/catalog/votu_no_NC_95_85/representatives_ge5kb.fna $IDX > logs/votu97.bowtie2_build.log 2>&1
map(){ local lib=$1 grp=$2 out=$3
  mkdir -p mapping/$out/bam mapping/$out/filtered
  $M/bowtie2 --very-sensitive --end-to-end -p $THREADS -x $IDX \
      -1 clean/$grp/${lib}_R1.fastq.gz -2 clean/$grp/${lib}_R2.fastq.gz 2> logs/$lib.votu97.bowtie2.log \
    | $M/samtools sort -@ 2 -m 1G -o mapping/$out/bam/$lib.sorted.bam - 2> logs/$lib.votu97.sort.log
  $M/samtools view -@ 4 -b -q 20 -F 2820 -e 'mapq < 255 && qlen > 0 && [NM] <= 0.05*qlen' \
    -o mapping/$out/filtered/$lib.filtered.bam mapping/$out/bam/$lib.sorted.bam
  $M/samtools index mapping/$out/bam/$lib.sorted.bam
  $M/samtools index mapping/$out/filtered/$lib.filtered.bam
}
for lib in $PILOT; do map $lib pilot no_NC; done
for lib in $COMPANIONS; do map $lib companions companions; done
