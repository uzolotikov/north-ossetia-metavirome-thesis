#!/usr/bin/env bash
# Весь анализ по порядку: конвейер 01–10, затем анализ 00–14. Шаги можно запускать и по одному.
# Итоговые таблицы анализа появляются в $OUT/analysis/report; их копии лежат в results/tables.
source "$(dirname "$0")/common.sh"
R=$ANALYSIS/report

for s in 01_fastp 02_metaspades 03_virsorter2 04_checkv 05_filter_candidates 06_cluster_votu \
         07_map_catalog 08_votu_coverage 09_markers 10_blast_references; do
  log "конвейер: $s"
  bash $PIPE/$s.sh
done

log "анализ: сборки и контаминант"
bash $AN/00_spades_nometa.sh
$PY $AN/00_asm_stats.py > $R/assembly_stats.tsv
bash $AN/01_runwide.sh
log "анализ: общий референс"
bash $AN/02a_build_ref.sh
bash $AN/02c_map.sh
bash $AN/03_genomad.sh
$PY $AN/02d_summarize.py          # классы контигов требуют geNomad по общему референсу
bash $AN/03_rarefy.sh
log "анализ: профаги, маркеры, хозяева"
bash $AN/04_prophage.sh
$PY $AN/05_16S_extract.py $ANALYSIS/contam/MV_16S.fna
$PY $AN/05_tables.py
bash $AN/06_crispr.sh
bash $AN/06_phist.sh
bash $AN/07_ani_contaminant.sh
bash $AN/08_16S_blast.sh
bash $AN/09_biragzang_check.sh
log "анализ: итоговые таблицы"
$PY $AN/10_votu_stats.py
$PY $AN/11_catalog_flags.py
$PY $AN/12_coverage_tables.py
$PY $AN/13_virus_origin.py
$PY $AN/14_novelty_lysogeny.py
log "готово: $R"
