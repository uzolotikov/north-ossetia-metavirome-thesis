#!/bin/bash
# Таксономия наиболее обильных генов 16S рРНК в каждой сборке фагома по базе NCBI 16S_ribosomal_RNA.
source "$(dirname "$0")/../common.sh"
W=$ANALYSIS
X=$ENVS/mvtools/bin
DB16=$DB/ncbi16S/16S_ribosomal_RNA
O=$W/report/MV_16S_blast.tsv
Q=$W/contam/q16s.fa

# В каждой сборке берутся три локуса 16S длиной >=800 б с наибольшим покрытием контига.
# Сами локусы — barrnap (шаг 9 конвейера) через 05_16S_extract.py, файл contam/MV_16S.fna.
$PY - "$W/contam/MV_16S.fna" > $Q <<'PY'
import re, sys
seqs, h = {}, None
for line in open(sys.argv[1]):
    line = line.rstrip()
    if line.startswith(">"):
        h = line[1:]; seqs[h] = ""
    else:
        seqs[h] += line
rows = []
for h, s in seqs.items():
    cov = float(re.search(r"cov([\d.]+)$", h).group(1))
    rows.append((h.split("|")[0], len(s), cov, h, s))
for samp in sorted({r[0] for r in rows}):
    cand = sorted([r for r in rows if r[0] == samp and r[1] >= 800], key=lambda r: -r[2])
    for r in cand[:3]:
        sys.stdout.write(">%s\n%s\n" % (r[3], r[4]))
PY

{
  echo "# most abundant (>=800 bp) 16S loci per MV assembly vs NCBI 16S_ribosomal_RNA"
  echo "# query ranked by assembly coverage of the host contig; top 3 subject hits each"
  printf "query\tlen\tcov\tpident\talen\tsubject\n"
  "$X/blastn" -task megablast -query $Q -db "$DB16" \
    -outfmt "6 qseqid pident length qlen stitle" -max_target_seqs 3 -num_threads 8 2>/dev/null \
  | awk -F'\t' 'BEGIN{OFS="\t"}{
      split($1,a,"|"); match($1,/len([0-9]+)\|/,L); match($1,/cov([0-9.]+)$/,C);
      print a[1], L[1], C[1], $2, $3, $5 }'
} > "$O"
rm -f $Q
cat "$O"
echo 16S_BLAST_DONE
