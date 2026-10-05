#!/bin/bash
# Taxonomy of the most abundant 16S rRNA loci in each MV assembly (NCBI 16S_ribosomal_RNA).
set -euo pipefail
W=/NatureUsers/glitch/mv_review
X=/NatureUsers/glitch/envs/mvtools/bin
DB=$W/db/ncbi16S/16S_ribosomal_RNA
O=$W/report/MV_16S_blast.tsv

# Pick, per MV assembly, the >=800 bp 16S loci ranked by contig coverage.
# 16S loci themselves come from barrnap via 05_16S_extract.py -> contam/MV_16S.fna
python3 - "$W/contam/MV_16S.fna" > /tmp/q16s_$$.fa <<'PY'
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
  "$X/blastn" -task megablast -query /tmp/q16s_$$.fa -db "$DB" \
    -outfmt "6 qseqid pident length qlen stitle" -max_target_seqs 3 -num_threads 8 2>/dev/null \
  | awk -F'\t' 'BEGIN{OFS="\t"}{
      split($1,a,"|"); match($1,/len([0-9]+)\|/,L); match($1,/cov([0-9.]+)$/,C);
      print a[1], L[1], C[1], $2, $3, $5 }'
} > "$O"
rm -f /tmp/q16s_$$.fa
cat "$O"
echo 16S_BLAST_DONE
