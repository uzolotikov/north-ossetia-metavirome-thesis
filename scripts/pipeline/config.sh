# Общие пути и наборы библиотек конвейера. Подключается всеми скриптами: source config.sh
P=/NatureUsers/nsinyavskiy/metaviroms_kgt            # корень проекта; все пути ниже — относительно него
RAW=/mss5/Reads/Metagenome/Ossetia-2020               # исходные прочтения прогона HY5G3DRXX
BIRAGZANG=/NatureUsers/oscypek/data/alania-meta/Biragzang.final.fna
PFAM=/NatureUsers/nsinyavskiy/db/Pfam/Pfam-A.hmm      # Pfam 38.2
CONDA=/NatureUsers/nsinyavskiy/miniconda3
ENVS=$CONDA/envs
THREADS=8
PY=/NatureUsers/glitch/envs/mvtools/bin/python   # Python 3.11 (окружение mvtools)

MV="MV-4122_S12 MV-4125_S18 MV-4128_S29"                                   # фагомы
PILOT="$MV NC1_S64"                                                        # фагомы и отрицательный контроль
COMPANIONS="4122-0-22_S11 4122-0-44_S10 4125-0-22_S17 4125-0-44_S16 4128-0-22_S28 4128-0-44_S27 4128_S30"
TIB_CELLULAR="4122-0-22_S11 4122-0-44_S10 4125-0-22_S17 4125-0-44_S16"   # клеточные фракции Тиб (собираются)
