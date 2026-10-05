#!/bin/bash
# Шаг 9. Маркеры бактериальной ДНК в сборках фагомов: гены 16S рРНК (barrnap 1.10.5) и рибосомные белки
# (prodigal 2.6.3 + hmmsearch HMMER 3.4 по 156 профилям Pfam 38.2, список — profile_names.txt).
set -euo pipefail
source "$(dirname "$0")/config.sh"; cd "$P"
H=$ENVS/ribosomal-env/bin
# 156 рибосомных профилей из Pfam-A (список в databases/ribosomal/profile_names.txt)
if [ ! -s databases/ribosomal/ribosomal.hmm ]; then
  ln -sf $PFAM databases/ribosomal/Pfam-A.hmm
  $H/hmmfetch --index databases/ribosomal/Pfam-A.hmm
  $H/hmmfetch -f databases/ribosomal/Pfam-A.hmm databases/ribosomal/profile_names.txt > databases/ribosomal/ribosomal.hmm
fi
mkdir -p qc/barrnap_no_NC
for s in $MV; do
  for kingdom in bac arc; do
    $ENVS/barrnap-env/bin/barrnap --kingdom $kingdom --threads 4 assembly/no_NC/$s/contigs.fasta \
      > qc/barrnap_no_NC/$s.$kingdom.gff 2> logs/$s.barrnap.$kingdom.log
  done
  R=qc/ribosomal_no_NC/$s; mkdir -p $R
  $H/prodigal -i assembly/no_NC/$s/contigs.fasta -a $R/proteins.faa -o $R/genes.gff -f gff -p meta > logs/$s.prodigal.log 2>&1
  $H/hmmsearch -o $R/hmmsearch.txt --tblout $R/ribosomal.tblout --domtblout $R/ribosomal.domtblout \
    --noali --cut_ga --cpu 4 databases/ribosomal/ribosomal.hmm $R/proteins.faa
done
