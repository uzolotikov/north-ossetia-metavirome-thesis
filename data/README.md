# Входные данные

Каталог пуст. Сюда кладутся прочтения и сборка Бирагзанга; в `config.example.sh` уже указаны пути
`data/reads` и `data/Biragzang.final.fna`.

| Что | Где взять |
|---|---|
| Сырые прочтения 11 библиотек, `reads/<библиотека>_R{1,2}_001.fastq.gz`, 20 ГБ | NCBI SRA, _[ЗАПОЛНИТЬ: номер BioProject и номера прочтений]_ |
| `Biragzang.final.fna` — бактериальная метагеномная сборка точки Бирагзанг, 67 Мб | _[ЗАПОЛНИТЬ: номер в NCBI или DOI на Zenodo]_ |

Список библиотек и их файлов — `metadata/libraries.tsv`. Контрольные суммы MD5 прочтений — в
`metadata/raw_md5.txt`:

```bash
cd data/reads && md5sum -c ../../metadata/raw_md5.txt
```

Скрипт `workflow/analysis/01_runwide.sh` проверяет все 56 библиотек прогона HY5G3DRXX, а не только 11 библиотек
работы. Если в `data/reads` лежат только они, таблица `runwide_presence.tsv` получится на 11 строк вместо 56.

Сборка Бирагзанга получена вне этой работы, и её построение в конвейер не входит. Её характеристики — в
`results/assembly_stats/assembly_stats.tsv`.
