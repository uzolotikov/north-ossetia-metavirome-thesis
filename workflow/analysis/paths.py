"""Рабочие каталоги для скриптов анализа на Python.

Значения берутся из переменных окружения, которые задаёт workflow/common.sh:
  OUT       рабочий каталог конвейера (OUT в config.sh)
  ANALYSIS  рабочий каталог анализа, $OUT/analysis
  BIRAGZANG бактериальная сборка точки Бирагзанг (BIRAGZANG в config.sh)
"""
import os
import sys

try:
    OUT = os.environ["OUT"]
    ANALYSIS = os.environ["ANALYSIS"]
    BIRAGZANG = os.environ["BIRAGZANG"]
except KeyError:
    sys.exit("Не заданы OUT, ANALYSIS и BIRAGZANG. Запускайте скрипт через workflow/run_all.sh "
             "или из оболочки после source workflow/common.sh.")
REPORT = os.path.join(ANALYSIS, "report")
