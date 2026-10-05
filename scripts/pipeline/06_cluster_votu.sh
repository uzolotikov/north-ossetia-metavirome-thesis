#!/bin/bash
# Шаг 6. Кластеризация в vOTU: BLASTN «все против всех» (BLAST+ 2.17.0), 95 % идентичности, 85 % покрытия.
set -euo pipefail
source "$(dirname "$0")/config.sh"; cd "$P"
B=$ENVS/blast-env/bin; IN=viral/filtered_no_NC/filtering_no_NC/strict_ge5kb.fna; D=votu/no_NC/blastdb
OUT=votu/no_NC/catalog/votu_no_NC_95_85
mkdir -p $D
$B/makeblastdb -in $IN -dbtype nucl -out $D/candidates_ge5kb > logs/votu_ge5kb.makeblastdb.log 2>&1
$B/blastn -task blastn -query $IN -db $D/candidates_ge5kb -num_threads $THREADS -evalue 1e-10 -max_target_seqs 10000 \
  -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore nident btop" \
  -out $D/all_vs_all_ge5kb.tsv 2> logs/votu_ge5kb.blastn.log
$PY "$(dirname "$0")/cluster_votu.py" $IN $D/all_vs_all_ge5kb.tsv $OUT
$PY "$(dirname "$0")/catalog_subsets.py" $OUT
