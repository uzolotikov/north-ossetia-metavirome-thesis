#!/bin/bash
# Рарефикация детекции vOTU: отфильтрованные BAM фагомов прореживаются до глубины самого бедного
# из них (MV-4125, 12,84 млн пар), три повтора; покрытие samtools coverage.
source "$(dirname "$0")/../common.sh"
set +e
M=$ENVS/mapping-env/bin
mkdir -p $ANALYSIS/matrix && cd $ANALYSIS/matrix
declare -A FR=( [MV-4122_S12]=.6744 [MV-4125_S18]=1 [MV-4128_S29]=.6325 )
: > rarefied_coverage.tsv
for s in $MV; do
  b=$OUT/mapping/no_NC/filtered/$s.filtered.bam
  $M/samtools coverage $b | awk -v s=$s -v r=full 'NR>1{print s"\t"r"\t"$1"\t"$3"\t"$6"\t"$7}' >> rarefied_coverage.tsv
  for seed in 1 2 3; do
    f=${FR[$s]}
    if [ "$f" = "1" ]; then t=$b; else t=tmp_${s}_$seed.bam; $M/samtools view -b -s ${seed}${f} -o $t $b; fi
    $M/samtools coverage $t | awk -v s=$s -v r=seed$seed 'NR>1{print s"\t"r"\t"$1"\t"$3"\t"$6"\t"$7}' >> rarefied_coverage.tsv
    [ "$t" != "$b" ] && rm -f $t
  done
done
echo RAREFY_DONE
