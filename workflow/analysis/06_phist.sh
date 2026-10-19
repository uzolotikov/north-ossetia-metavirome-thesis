#!/bin/bash
# Связи вирус — хозяин по общим k-мерам (PHIST 1.0.0): каждая vOTU против каждого контига >=5 кб
# клеточных сборок Тиб и сборки Бирагзанга.
source "$(dirname "$0")/../common.sh"
set +e
export PATH=$ENVS/mvtools/bin:$PATH
W=$ANALYSIS
D=$W/hosts/phist; rm -rf "${D:?}"; mkdir -p $D/virus $D/host $D/out
awk 'BEGIN{RS=">"} NR>1{n=index($0,"\n"); h=substr($0,1,n-1); split(h,a," "); f=D"/virus/"a[1]".fna"; printf ">%s\n%s",h,substr($0,n+1) > f; close(f)}' D=$D $OUT/votu/no_NC/catalog/votu_no_NC_95_85/representatives_ge5kb.fna
ext(){ awk -v P=$2 -v D=$D 'BEGIN{RS=">"} NR>1{n=index($0,"\n"); h=substr($0,1,n-1); split(h,a," "); s=substr($0,n+1); t=s; gsub(/\n/,"",t); if(length(t)>=5000){f=D"/host/"P"__"a[1]".fna"; printf ">%s__%s\n%s",P,a[1],s > f; close(f)}}' $1; }
export D
for a in $TIB_CELLULAR; do ext $OUT/assembly/cellular/$a/contigs.fasta cell-$a; done
ext $BIRAGZANG Biragzang
echo "viruses: $(ls $D/virus | wc -l) hosts: $(ls $D/host | wc -l)"
phist.py -t $THREADS $D/virus $D/host $D/out/common_kmers.csv $D/out/predictions.csv > $W/logs/phist_run.log 2>&1 && echo PHIST_DONE || echo PHIST_FAIL
