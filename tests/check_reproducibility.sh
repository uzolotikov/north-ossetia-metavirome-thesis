#!/bin/bash
# Пересчёт результатов конвейера, не требующих сборки и картирования, во временной папке
# и побайтное сравнение с файлами проекта: отбор кандидатов, каталог vOTU, таблица покрытия,
# набор рибосомных профилей, сравнения с клеточными сборками и Biragzang.final.
# Использование: tests/check_reproducibility.sh [рабочая_папка]
set -uo pipefail
R=$(cd "$(dirname "$0")/.." && pwd)
source "$R/scripts/pipeline/config.sh"
SC=$R/scripts/pipeline
T=${1:-$(mktemp -d)}; mkdir -p $T/db
B=$ENVS/blast-env/bin; H=$ENVS/ribosomal-env/bin
pass=0; fail=0
same(){ if cmp -s "$1" "$2"; then echo "  ок        $3"; pass=$((pass+1)); else echo "  РАЗЛИЧИЕ  $3"; fail=$((fail+1)); fi; }
same_unordered(){ if diff -q <(sort "$1") <(sort "$2") >/dev/null; then echo "  ок        $3 (порядок строк не учитывается)"; pass=$((pass+1)); else echo "  РАЗЛИЧИЕ  $3"; fail=$((fail+1)); fi; }

echo "Шаг 5. Отбор кандидатов"
$PY $SC/filter_candidates.py $P/checkv/no_NC $P/viral/no_NC $T/filtering > /dev/null
F=$P/viral/filtered_no_NC/filtering_no_NC
for f in all_after_checkv.fna strict_supported.fna review.fna exclude_strict.fna strict_ge5kb.fna strict_ge10kb.fna short_HQ.fna \
         candidate_decisions.tsv sequence_manifest.tsv summary.tsv; do same $T/filtering/$f $F/$f $f; done

echo "Шаг 6. Кластеризация в vOTU"
$B/makeblastdb -in $T/filtering/strict_ge5kb.fna -dbtype nucl -out $T/db/candidates_ge5kb > /dev/null
$B/blastn -task blastn -query $T/filtering/strict_ge5kb.fna -db $T/db/candidates_ge5kb -num_threads $THREADS -evalue 1e-10 \
  -max_target_seqs 10000 -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore nident btop" \
  -out $T/all_vs_all_ge5kb.tsv
$PY $SC/cluster_votu.py $T/filtering/strict_ge5kb.fna $T/all_vs_all_ge5kb.tsv $T/catalog > /dev/null
$PY $SC/catalog_subsets.py $T/catalog
C=$P/votu/no_NC/catalog/votu_no_NC_95_85
for f in clusters.tsv membership.tsv representatives_ge5kb.fna representatives_ge10kb.fna \
         assembly_presence_ge5000.tsv assembly_presence_ge10000.tsv; do same $T/catalog/$f $C/$f $f; done
same_unordered $T/catalog/pairwise_metrics.tsv $C/pairwise_metrics.tsv pairwise_metrics.tsv

echo "Шаг 8. Покрытие каталога в 11 библиотеках"
SAM=$ENVS/mapping-env/bin/samtools
for x in $(for l in $PILOT; do echo $l:no_NC; done) $(for l in $COMPANIONS; do echo $l:companions; done); do
  l=${x%%:*}; g=${x##*:}; mkdir -p $T/proj/mapping/$g/coverage
  $PY $SC/coverage_tables.py $SAM $P/mapping/$g/filtered/$l.filtered.bam $l $T/proj/mapping/$g/coverage
  same $T/proj/mapping/$g/coverage/$l.coverage.tsv $P/mapping/$g/coverage/$l.coverage.tsv $l.coverage.tsv
  same $T/proj/mapping/$g/coverage/$l.idxstats.tsv $P/mapping/$g/coverage/$l.idxstats.tsv $l.idxstats.tsv
done
$PY $SC/rebuild_votu97.py $T/proj $T/proj/votu97 > /dev/null
same $T/proj/votu97/votu97_all_libraries.tsv $P/qc/votu97_comparison/votu97_all_libraries.tsv votu97_all_libraries.tsv
same $T/proj/votu97/library_summary.tsv $P/qc/votu97_comparison/library_summary.tsv library_summary.tsv

echo "Шаг 9. Набор рибосомных профилей"
ln -sf $PFAM $T/Pfam-A.hmm
$H/hmmfetch --index $T/Pfam-A.hmm > /dev/null
$H/hmmfetch -f $T/Pfam-A.hmm $P/databases/ribosomal/profile_names.txt > $T/ribosomal.hmm
same $T/ribosomal.hmm $P/databases/ribosomal/ribosomal.hmm ribosomal.hmm

echo "Шаг 10. Сравнение с бактериальными сборками"
for s in $TIB_CELLULAR; do
  $B/makeblastdb -in $P/assembly/cellular/$s/contigs.fasta -dbtype nucl -out $T/db/$s > /dev/null
  $B/blastn -task blastn -query $C/representatives_ge5kb.fna -db $T/db/$s -num_threads $THREADS -evalue 1e-10 -max_target_seqs 10000 \
    -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore nident btop" -out $T/votu97_vs_$s.tsv
  same_unordered $T/votu97_vs_$s.tsv $P/qc/cellular_comparison/votu97_vs_$s.tsv votu97_vs_$s.tsv
done
awk '/^>/{keep=($0 ~ /^>MV-4128_S29__/)} keep' $F/strict_supported.fna > $T/MV-4128.strict_supported.fna
$B/makeblastdb -in $BIRAGZANG -dbtype nucl -out $T/db/Biragzang > /dev/null
$B/blastn -task blastn -query $T/MV-4128.strict_supported.fna -db $T/db/Biragzang -num_threads $THREADS -evalue 1e-10 \
  -max_target_seqs 10000 -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore" \
  -out $T/MV-4128_vs_Biragzang.tsv
same_unordered $T/MV-4128_vs_Biragzang.tsv $P/qc/cellular_comparison/MV-4128_vs_Biragzang.tsv MV-4128_vs_Biragzang.tsv

echo; echo "Совпало: $pass, различий: $fail"; [ $fail -eq 0 ]
