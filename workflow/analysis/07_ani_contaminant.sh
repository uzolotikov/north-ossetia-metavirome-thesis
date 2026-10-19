#!/bin/bash
# Идентичность MV-специфичного контаминанта между сборками MV-4125 и MV-4128.
# Таблицу step6/allvsall_99.tsv для этого использовать нельзя: она построена с -perc_identity 99, то есть
# содержит только HSP с идентичностью >=99 % и завышает оценку. Здесь megablast запускается без фильтра.
source "$(dirname "$0")/../common.sh"
W=$ANALYSIS
B=$ENVS/blast-env/bin
O=$W/ani
mkdir -p "$O"
cd "$O"

# Запрос: 46 контигов MV-4128 длиной >=50 кб, на которые приходится около 70 % прочтений MV-4128 (из 01_runwide.sh)
QUERY=$W/contam/ref_MV4128_ge50kb.fna
# База: вся сборка metaSPAdes MV-4125
SUBJ=$OUT/assembly/no_NC/MV-4125_S18/contigs.fasta

if [ ! -f mv4125.nsq ]; then
  "$B/makeblastdb" -in "$SUBJ" -dbtype nucl -out mv4125 > makeblastdb.log 2>&1
fi

# Без -perc_identity; nident запрашивается, чтобы идентичность считалась точно, без округления pident.
"$B/blastn" -task megablast -query "$QUERY" -db mv4125 \
  -evalue 1e-20 -num_threads $THREADS -max_target_seqs 20 \
  -outfmt "6 qseqid sseqid pident nident length mismatch gapopen qstart qend sstart send qlen slen bitscore" \
  > mv4128_ge50kb_vs_mv4125.tsv 2> blastn.log

$PY "$AN/07_ani_calc.py" mv4128_ge50kb_vs_mv4125.tsv "$QUERY" \
  > "$W/report/contaminant_ani.tsv"
cat "$W/report/contaminant_ani.tsv"
echo ANI_DONE
