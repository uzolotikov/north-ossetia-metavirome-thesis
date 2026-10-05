#!/bin/bash
set -euo pipefail
K=/NatureUsers/nsinyavskiy/metaviroms_kgt; W=/NatureUsers/glitch/mv_review/step6
B=/NatureUsers/nsinyavskiy/miniconda3/envs/blast-env/bin
M=/NatureUsers/nsinyavskiy/miniconda3/envs/mapping-env/bin
cd $W
ext(){ awk -v M=5000 -v P=$2 'BEGIN{RS=">"} NR>1{n=index($0,"\n"); h=substr($0,1,n-1); split(h,hh," "); s=substr($0,n+1); gsub(/\n/,"",s); if(length(s)>=M) printf ">%s__%s\n%s\n",P,hh[1],s}' $1; }
: > ref_all.fna
for s in MV-4122_S12 MV-4125_S18 MV-4128_S29; do ext $K/assembly/no_NC/$s/contigs.fasta $s >> ref_all.fna; done
for s in 4122-0-22_S11 4122-0-44_S10 4125-0-22_S17 4125-0-44_S16; do ext $K/assembly/cellular/$s/contigs.fasta cell-$s >> ref_all.fna; done
ext /NatureUsers/oscypek/data/alania-meta/Biragzang.final.fna Biragzang >> ref_all.fna
grep -c ">" ref_all.fna
$B/makeblastdb -in ref_all.fna -dbtype nucl -out ref_all > /dev/null
$B/blastn -task megablast -query ref_all.fna -db ref_all -perc_identity 99 -evalue 1e-20 -num_threads 32 -max_target_seqs 50 \
  -outfmt "6 qseqid sseqid pident length qstart qend sstart send qlen slen" > allvsall_99.tsv
python3 /NatureUsers/glitch/mv_review/scripts/02b_derep.py ref_all.fna allvsall_99.tsv ref_derep.fna derep_map.tsv
$M/bowtie2-build --threads 32 ref_derep.fna ref_derep > bt2build.log 2>&1
echo REF_DONE
