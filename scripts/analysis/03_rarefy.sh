#!/bin/bash
W=/NatureUsers/glitch/mv_review; K=/NatureUsers/nsinyavskiy/metaviroms_kgt; M=/NatureUsers/nsinyavskiy/miniconda3/envs/mapping-env/bin
cd $W/matrix
declare -A FR=( [MV-4122_S12]=.6744 [MV-4125_S18]=1 [MV-4128_S29]=.6325 )
: > rarefied_coverage.tsv
for s in MV-4122_S12 MV-4125_S18 MV-4128_S29; do
  b=$K/mapping/no_NC/filtered/$s.filtered.bam
  $M/samtools coverage $b | awk -v s=$s -v r=full 'NR>1{print s"\t"r"\t"$1"\t"$3"\t"$6"\t"$7}' >> rarefied_coverage.tsv
  for seed in 1 2 3; do
    f=${FR[$s]}
    if [ "$f" = "1" ]; then t=$b; else t=tmp_${s}_$seed.bam; $M/samtools view -b -s ${seed}${f} -o $t $b; fi
    $M/samtools coverage $t | awk -v s=$s -v r=seed$seed 'NR>1{print s"\t"r"\t"$1"\t"$3"\t"$6"\t"$7}' >> rarefied_coverage.tsv
    [ "$t" != "$b" ] && rm -f $t
  done
done
echo RAREFY_DONE
