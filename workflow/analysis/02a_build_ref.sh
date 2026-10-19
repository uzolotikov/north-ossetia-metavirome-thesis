#!/bin/bash
# Общий референс для картирования: контиги >=5 кб трёх сборок фагомов, четырёх клеточных сборок Тиб
# и сборки Бирагзанга; megablast «все против всех» (>=99 %) и дерепликация (02b_derep.py),
# затем индекс bowtie2. Картирование — 02c_map.sh.
source "$(dirname "$0")/../common.sh"
B=$ENVS/blast-env/bin
M=$ENVS/mapping-env/bin
W=$ANALYSIS/step6
mkdir -p $W && cd $W
ext(){ awk -v M=5000 -v P=$2 'BEGIN{RS=">"} NR>1{n=index($0,"\n"); h=substr($0,1,n-1); split(h,hh," "); s=substr($0,n+1); gsub(/\n/,"",s); if(length(s)>=M) printf ">%s__%s\n%s\n",P,hh[1],s}' $1; }
: > ref_all.fna
for s in $MV; do ext $OUT/assembly/no_NC/$s/contigs.fasta $s >> ref_all.fna; done
for s in $TIB_CELLULAR; do ext $OUT/assembly/cellular/$s/contigs.fasta cell-$s >> ref_all.fna; done
ext $BIRAGZANG Biragzang >> ref_all.fna
grep -c ">" ref_all.fna
$B/makeblastdb -in ref_all.fna -dbtype nucl -out ref_all > /dev/null
$B/blastn -task megablast -query ref_all.fna -db ref_all -perc_identity 99 -evalue 1e-20 -num_threads $THREADS -max_target_seqs 50 \
  -outfmt "6 qseqid sseqid pident length qstart qend sstart send qlen slen" > allvsall_99.tsv
$PY $AN/02b_derep.py ref_all.fna allvsall_99.tsv ref_derep.fna derep_map.tsv
$M/bowtie2-build --threads $THREADS ref_derep.fna ref_derep > bt2build.log 2>&1
echo REF_DONE
