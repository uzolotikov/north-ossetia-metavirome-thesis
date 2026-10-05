#!/bin/bash
# Шаг 10. Сравнение с бактериальными сборками тех же точек (BLAST+ 2.17.0, E <= 1e-10):
# каталог vOTU — с четырьмя клеточными сборками Тиб; кандидаты MV-4128 — со сборкой Biragzang.final.
set -euo pipefail
source "$(dirname "$0")/config.sh"; cd "$P"
B=$ENVS/blast-env/bin; Q=votu/no_NC/catalog/votu_no_NC_95_85/representatives_ge5kb.fna
mkdir -p databases/blast_cellular databases/blast_biragzang qc/cellular_comparison
for s in $TIB_CELLULAR; do
  $B/makeblastdb -in assembly/cellular/$s/contigs.fasta -dbtype nucl -out databases/blast_cellular/$s > logs/$s.makeblastdb.log 2>&1
  $B/blastn -task blastn -query $Q -db databases/blast_cellular/$s -num_threads $THREADS -evalue 1e-10 -max_target_seqs 10000 \
    -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore nident btop" \
    -out qc/cellular_comparison/votu97_vs_$s.tsv 2> logs/$s.votu97_blastn.log
done
awk '/^>/{keep=($0 ~ /^>MV-4128_S29__/)} keep' viral/filtered_no_NC/filtering_no_NC/strict_supported.fna \
  > qc/cellular_comparison/MV-4128.strict_supported.fna
$B/makeblastdb -in $BIRAGZANG -dbtype nucl -out databases/blast_biragzang/Biragzang > logs/Biragzang.makeblastdb.log 2>&1
$B/blastn -task blastn -query qc/cellular_comparison/MV-4128.strict_supported.fna -db databases/blast_biragzang/Biragzang \
  -num_threads $THREADS -evalue 1e-10 -max_target_seqs 10000 \
  -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore" \
  -out qc/cellular_comparison/MV-4128_vs_Biragzang.tsv 2> logs/MV-4128_vs_Biragzang.blastn.log
