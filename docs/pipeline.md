# Конвейер обработки

Все команды находятся в скриптах; этот документ описывает порядок шагов, ключевые параметры и выходы. Конвейер (`scripts/pipeline/`) ведёт от сырых прочтений до каталога vOTU и таблиц покрытия; анализ (`scripts/analysis/`) использует эти результаты для оценки чистоты фракции, происхождения фагов и связей с хозяевами.

Пути и наборы библиотек задаются в `scripts/pipeline/config.sh`. Скрипты конвейера запускаются из любого места и работают в корне проекта на сервере labg (`/NatureUsers/nsinyavskiy/metaviroms_kgt`); скрипты анализа пишут в рабочую папку `/NatureUsers/glitch/mv_review`, итоговые таблицы которой собраны в `results/tables/`.

## Программы

| Программа | Версия | Окружение (`envs/`) |
|---|---|---|
| fastp | 0.23.2 | `fastp.yml` |
| SPAdes / metaSPAdes | 4.3.0 | `metaspades-env.yml` |
| VirSorter2 | 2.2.4 | `virsorter2-env.yml` |
| CheckV | 1.1.1 | `checkv-env.yml` |
| BLAST+ | 2.17.0 | `blast-env.yml`, `mvtools.yml` |
| bowtie2 | 2.5.5 | `mapping-env.yml` |
| samtools | 1.24 | `mapping-env.yml` |
| prodigal | 2.6.3 | `ribosomal-env.yml` |
| HMMER | 3.4 | `ribosomal-env.yml` |
| barrnap | 1.10.5 | `barrnap-env.yml` |
| geNomad | 1.12.0 | `mvtools.yml` |
| PHIST | 1.0.0 | `mvtools.yml` |
| MinCED | 0.4.2 | `mvtools.yml` |
| MEGAHIT | 1.2.9 | `mvtools.yml` |
| MMseqs2 | 18.8cc5c | `mvtools.yml` |
| Python | 3.11 | `mvtools.yml` |

Для каждого окружения приведены спецификация `*.yml` и, где доступна, точная спецификация пакетов `*.explicit.txt`.

## Базы и внешние данные

| Ресурс | Версия | Где используется |
|---|---|---|
| База VirSorter2 | установлена `virsorter setup` 22.09.2026 | шаг 3 |
| База CheckV | 1.5 | шаги 4, 5 |
| Pfam-A | 38.2 (30 134 семейства) | шаг 9 |
| База geNomad | 1.9 (Zenodo, запись 14886553) | анализ: `03_genomad.sh` |
| NCBI 16S_ribosomal_RNA | загружена 02.10.2026 | анализ: `08_16S_blast.sh` |
| NCBI nt / nr, удалённый BLAST | обращения 02.10.2026 | таксономия контигов-хозяев и рибосомных белков |
| `Biragzang.final.fna` | бактериальная метагеномная сборка точки Бирагзанг (oscypek, 2023) | шаг 10; анализ: `02a_*`, `06_*`, `09_*` |

## Шаги конвейера

| Шаг | Скрипт | Что делает | Ключевые параметры | Выход (в корне проекта) |
|---|---|---|---|---|
| 1 | `01_fastp.sh` | тримминг всех 11 библиотек | обрезка адаптеров и poly-G, скользящее окно 4 б с качеством 20, минимальная длина 75 б | `clean/{pilot,companions}/`, `qc/fastp_*/` |
| 2 | `02_metaspades.sh` | сборка трёх фагомов и четырёх клеточных фракций Тиб | metaSPAdes, `-m 32` | `assembly/no_NC/`, `assembly/cellular/` |
| 3 | `03_virsorter2.sh` | поиск вирусных контигов | группы dsDNAphage, NCLDV, ssDNA, lavidaviridae; длина ≥1500 б; score ≥0,5 | `viral/no_NC/` |
| 4 | `04_checkv.sh` | качество, полнота, обрезка провирусных краёв | `checkv end_to_end` | `checkv/no_NC/` |
| 5 | `05_filter_candidates.sh` | отбор кандидатов и повторный CheckV по строгому набору | правила — `filter_candidates_README.md` | `viral/filtered_no_NC/filtering_no_NC/`, `checkv/filtered_no_NC/strict_supported_pass2/` |
| 6 | `06_cluster_votu.sh` | кластеризация в vOTU | BLASTN `-task blastn`, E ≤10⁻¹⁰; 95 % идентичности, 85 % покрытия — `cluster_votu_README.md` | `votu/no_NC/catalog/votu_no_NC_95_85/` |
| 7 | `07_map_catalog.sh` | картирование 11 библиотек на каталог | bowtie2 `--very-sensitive --end-to-end`; MAPQ ≥20, первичные выравнивания, ≤5 % несовпадений | `mapping/{no_NC,companions}/{bam,filtered}/` |
| 8 | `08_votu_coverage.sh` | покрытие каталога в каждой библиотеке и сводная таблица | глубина — `samtools depth -a -s -q 20`; детекция — ширина ≥75 % и средняя глубина ≥1× | `mapping/*/coverage/`, `qc/votu97_comparison/` |
| 9 | `09_markers.sh` | гены 16S рРНК и рибосомные белки в сборках фагомов | barrnap для бактерий и архей; prodigal `-p meta`; hmmsearch `--cut_ga` по 156 профилям Pfam | `qc/barrnap_no_NC/`, `qc/ribosomal_no_NC/` |
| 10 | `10_blast_references.sh` | сравнение с бактериальными сборками тех же точек | BLASTN `-task blastn`, E ≤10⁻¹⁰ | `qc/cellular_comparison/` |

**Детали шагов.**

- **Шаг 8.** Число прочтений на vOTU — по `samtools idxstats` отфильтрованного BAM. Средняя глубина и ширина покрытия — по `samtools depth -a -s -q 20`: учитываются все позиции, перекрывающиеся мейты одной пары считаются один раз, основания с качеством ниже 20 не учитываются. RPM и RPKM нормированы на число прочтений после тримминга (оба мейта пары) по отчётам fastp.
- **Шаг 9.** Набор из 156 рибосомных профилей извлекается из Pfam-A по списку `databases/ribosomal/profile_names.txt` (`hmmfetch -f`). Из 162 кандидатных профилей исключены шесть, связанных с рибосомой, но не являющихся собственно рибосомными белками либо с неоднозначным составом семейства: AROS, Morc6_S5, PrmA, Ribosomal_L7Ae, Ribosomal_S30AE, Ribosom_S30AE_C.
- **Шаг 10.** Таблицы содержат все совпадения BLAST с E ≤10⁻¹⁰; пороги идентичности и покрытия для интерпретации задаются в анализе (`09_biragzang_eval.py`).

## Анализ

Скрипты запускаются в порядке номеров после шагов 1–10.

| Скрипт | Назначение | Таблицы в `results/` |
|---|---|---|
| `00_install_env.sh` | окружение mvtools (micromamba) | `envs/mvtools.yml` |
| `00_spades_one.sh`, `00_spades_nometa.sh`, `00_asm_stats.py` | сборка фагомов SPAdes без `--meta` и сравнение с metaSPAdes | `assembly_stats/assembly_stats.tsv` |
| `01_runwide.sh` | присутствие MV-специфичного материала во всех 56 библиотеках прогона | `contamination/runwide_presence.tsv` |
| `02a_build_ref.sh` → `02b_derep.py` → `02c_map.sh` → `02d_summarize.py` | общий референс (8 754 контига ≥5 кб из всех сборок, дерепликация 99 % / 90 %), картирование 11 библиотек, классы контигов | `step6/`, `tables/step6_*.tsv` |
| `03_genomad.sh` | geNomad по каталогу и по общему референсу | `genomad/` |
| `03_rarefy.sh` | рарефикация детекции vOTU до 12,84 млн пар | `matrix/rarefied_coverage.tsv` |
| `04_prophage.sh` | vOTU, встроенные в бактериальные контиги; отношение покрытия «внутри / фланги» | `prophage/` |
| `05_16S_extract.py`, `05_tables.py` | локусы 16S рРНК; таблицы качества и чистоты фракции | `contamination/MV_16S.fna`, `tables/qc_table.tsv`, `tables/purity_table.tsv` |
| `06_crispr.sh`, `06_phist.sh` | спейсеры CRISPR; связи вирус — хозяин по общим k-мерам | `hosts/` |
| `07_ani_contaminant.sh` → `07_ani_calc.py` | идентичность MV-специфичного контаминанта между MV-4125 и MV-4128 | `tables/contaminant_ani.tsv` |
| `08_16S_blast.sh` | таксономия наиболее обильных генов 16S рРНК | `tables/MV_16S_blast.tsv` |
| `09_biragzang_check.sh` → `09_biragzang_eval.py` | совпадения каталога и кандидатов MV-4128 с `Biragzang.final` при разных порогах | `tables/biragzang_match_check.tsv` |
| `10_votu_stats.py` | отношение обилий контаминанта, общие vOTU, подтверждение клеточной фракцией, версии матрицы детекции | `tables/votu_*.tsv` |
| `11_catalog_flags.py` | помеченный каталог с рекомендацией по каждой vOTU | `tables/votu97_catalog_flags.tsv` |
| `12_coverage_tables.py` | сводка покрытия по библиотекам и матрица 97 × 11 | `tables/votu_coverage_*.tsv` |
| `13_virus_origin.py` | категория происхождения каждой vOTU | `tables/votu_origin.tsv` |
| `14_novelty_lysogeny.py` | новизна (AAI к базе CheckV) и гены лизогении | `tables/votu_novelty_lysogeny.tsv` |

Таксономия контигов-хозяев и рибосомных белков получена удалённым BLAST в NCBI; запросы и ответы сохранены: `results/contamination/top_rp.faa`, `top_rp_remote.tsv`, `results/hosts/phist_strong_hosts.fna`, `phist_strong_hosts_nt.tsv`, `results/hosts/remote_blast/`.

## Проверка воспроизводимости

`tests/check_reproducibility.sh` пересчитывает во временной папке шаги, не требующие сборки и картирования, — отбор кандидатов, кластеризацию, таблицы покрытия по отфильтрованным BAM, набор рибосомных профилей и сравнения с бактериальными сборками — и побайтно сравнивает результат с файлами проекта (47 файлов). Последний запуск 04.10.2026: все 47 файлов совпадают.
