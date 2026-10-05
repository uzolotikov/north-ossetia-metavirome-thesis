#!/bin/bash
# Reproduces §3.1: identity of the MV-specific contaminant contigs between MV-4125 and MV-4128.
# NOTE: step6/allvsall_99.tsv must NOT be used for this - it was built with -perc_identity 99,
# which keeps only >=99% HSPs and biases identity upward. This script runs an unfiltered megablast.
set -euo pipefail
K=/NatureUsers/nsinyavskiy/metaviroms_kgt
W=/NatureUsers/glitch/mv_review
B=/NatureUsers/nsinyavskiy/miniconda3/envs/blast-env/bin
O=$W/ani
mkdir -p "$O"
cd "$O"

# Query: the 46 MV-4128 contigs >=50 kb that carry ~70% of MV-4128 reads (built in 01_runwide.sh)
QUERY=$W/contam/ref_MV4128_ge50kb.fna
# Subject: the full MV-4125 metaSPAdes assembly
SUBJ=$K/assembly/no_NC/MV-4125_S18/contigs.fasta

if [ ! -f mv4125.nsq ]; then
  "$B/makeblastdb" -in "$SUBJ" -dbtype nucl -out mv4125 > makeblastdb.log 2>&1
fi

# No -perc_identity filter. nident is requested so identity is counted, not rounded.
"$B/blastn" -task megablast -query "$QUERY" -db mv4125 \
  -evalue 1e-20 -num_threads 16 -max_target_seqs 20 \
  -outfmt "6 qseqid sseqid pident nident length mismatch gapopen qstart qend sstart send qlen slen bitscore" \
  > mv4128_ge50kb_vs_mv4125.tsv 2> blastn.log

python3 "$W/scripts/07_ani_calc.py" mv4128_ge50kb_vs_mv4125.tsv "$QUERY" \
  > "$W/report/contaminant_ani.tsv"
cat "$W/report/contaminant_ani.tsv"
echo ANI_DONE
