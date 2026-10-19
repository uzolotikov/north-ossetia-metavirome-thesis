# Настройки конвейера. Скопируйте файл в config.sh и укажите свои пути:
#   cp config.example.sh config.sh
# Относительные пути считаются от корня репозитория.

THREADS=16

# Исходные прочтения прогона HY5G3DRXX: <библиотека>_R1_001.fastq.gz и <библиотека>_R2_001.fastq.gz.
# Какие файлы нужны, описано в data/README.md.
RAW=data/reads

# Бактериальная метагеномная сборка точки Бирагзанг (см. data/README.md).
BIRAGZANG=data/Biragzang.final.fna

# Каталог баз: VirSorter2, CheckV 1.5, geNomad 1.9, NCBI 16S_ribosomal_RNA. Состав — в README.
DB=db

# Рабочий каталог: очищенные прочтения, сборки, BAM и промежуточные таблицы, около 80 ГБ.
# Скрипты анализа пишут в его подкаталог analysis/.
OUT=out

# Каталог с окружениями conda из envs/; программы вызываются как $ENVS/<окружение>/bin/<программа>.
ENVS="$(conda info --base)/envs"
