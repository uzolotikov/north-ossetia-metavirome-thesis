#!/bin/bash
set -uo pipefail
B=/NatureUsers/nsinyavskiy/miniconda3/envs/blast-env/bin
A=/NatureUsers/nsinyavskiy/metaviroms_kgt/assembly/no_NC
R=/mss5/Reads/Metagenome/Ossetia-2020
O=/NatureUsers/glitch/mv_review/contam
ext(){ awk -v M=$2 'BEGIN{RS=">"} NR>1{n=index($0,"\n"); h=substr($0,1,n-1); s=substr($0,n+1); gsub(/\n/,"",s); if(length(s)>=M) printf ">%s\n%s\n",h,s}' $1; }
ext $A/MV-4128_S29/contigs.fasta 50000 | sed 's/^>/>MV4128big_/' > $O/ref_MV4128_ge50kb.fna
ext $A/MV-4122_S12/contigs.fasta 50000 | sed 's/^>/>MV4122big_/' > $O/ref_MV4122_ge50kb.fna
cat $O/ref_MV4128_ge50kb.fna $O/ref_MV4122_ge50kb.fna > $O/ref_runwide.fna
$B/makeblastdb -in $O/ref_runwide.fna -dbtype nucl -out $O/ref_runwide > /dev/null
N=20000
echo -e "library\treads_tested\thit_MV4128big\thit_MV4122big\tpct_MV4128big\tpct_MV4122big" > $O/runwide_presence.tsv
for f in $(ls $R/*_R1_001.fastq.gz $R/*_1.fastq.gz 2>/dev/null); do
  lib=$(basename $f | sed 's/_R1_001.fastq.gz//; s/_1.fastq.gz//')
  zcat $f | head -n $((N*4)) | awk 'NR%4==1{print ">"substr($1,2)} NR%4==2{print}' > $O/tmp_q.fa
  n=$(grep -c ">" $O/tmp_q.fa)
  $B/blastn -task megablast -query $O/tmp_q.fa -db $O/ref_runwide -outfmt "6 qseqid sseqid pident length" -evalue 1e-10 -max_target_seqs 1 -num_threads 8 2>/dev/null \
   | awk '$3>=97 && $4>=100' | sort -u -k1,1 > $O/tmp_hits.tsv
  a=$(grep -c "MV4128big_" $O/tmp_hits.tsv); b=$(grep -c "MV4122big_" $O/tmp_hits.tsv)
  awk -v l=$lib -v n=$n -v a=$a -v b=$b 'BEGIN{printf "%s\t%d\t%d\t%d\t%.3f\t%.3f\n",l,n,a,b,100*a/n,100*b/n}' >> $O/runwide_presence.tsv
done
rm -f $O/tmp_q.fa $O/tmp_hits.tsv
echo RUNWIDE_DONE
