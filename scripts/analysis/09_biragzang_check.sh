#!/bin/bash
# Matches of the vOTU catalogue and of the MV-4128 virus candidates against Biragzang.final,
# reported for a grid of identity / query-coverage thresholds.
set -euo pipefail
K=/NatureUsers/nsinyavskiy/metaviroms_kgt
W=/NatureUsers/glitch/mv_review
B=/NatureUsers/nsinyavskiy/miniconda3/envs/blast-env/bin
O=$W/biragzang
mkdir -p "$O"; cd "$O"

VOTU=$K/votu/no_NC/catalog/votu_no_NC_95_85/representatives_ge5kb.fna
BIRA=/NatureUsers/oscypek/data/alania-meta/Biragzang.final.fna

if [ ! -f bira.nsq ]; then
  "$B/makeblastdb" -in "$BIRA" -dbtype nucl -out bira > makeblastdb.log 2>&1
fi
"$B/blastn" -task megablast -query "$VOTU" -db bira \
  -evalue 1e-20 -num_threads 16 -max_target_seqs 20 \
  -outfmt "6 qseqid sseqid pident nident length qstart qend qlen slen bitscore" \
  > votu97_vs_Biragzang_full.tsv 2> blastn.log

python3 "$W/scripts/09_biragzang_eval.py" \
  votu97_vs_Biragzang_full.tsv \
  "$K/qc/cellular_comparison/MV-4128_vs_Biragzang.tsv" \
  > "$W/report/biragzang_match_check.tsv"
cat "$W/report/biragzang_match_check.tsv"
echo BIRAGZANG_DONE
