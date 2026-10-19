#!/bin/bash
# Совпадения каталога vOTU и вирусных кандидатов MV-4128 со сборкой Biragzang.final
# при разных порогах идентичности и покрытия запроса (оценка — 09_biragzang_eval.py).
source "$(dirname "$0")/../common.sh"
W=$ANALYSIS
B=$ENVS/blast-env/bin
O=$W/biragzang
mkdir -p "$O"; cd "$O"

VOTU=$OUT/votu/no_NC/catalog/votu_no_NC_95_85/representatives_ge5kb.fna
BIRA=$BIRAGZANG

if [ ! -f bira.nsq ]; then
  "$B/makeblastdb" -in "$BIRA" -dbtype nucl -out bira > makeblastdb.log 2>&1
fi
"$B/blastn" -task megablast -query "$VOTU" -db bira \
  -evalue 1e-20 -num_threads $THREADS -max_target_seqs 20 \
  -outfmt "6 qseqid sseqid pident nident length qstart qend qlen slen bitscore" \
  > votu97_vs_Biragzang_full.tsv 2> blastn.log

$PY "$AN/09_biragzang_eval.py" \
  votu97_vs_Biragzang_full.tsv \
  "$OUT/qc/cellular_comparison/MV-4128_vs_Biragzang.tsv" \
  > "$W/report/biragzang_match_check.tsv"
cat "$W/report/biragzang_match_check.tsv"
echo BIRAGZANG_DONE
