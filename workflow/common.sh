# Общие настройки и наборы библиотек: читается всеми скриптами конвейера и анализа.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
if [ ! -f config.sh ]; then
  echo "Нет config.sh. Скопируйте config.example.sh в config.sh и укажите пути." >&2
  exit 1
fi
# shellcheck source=/dev/null
source config.sh

# Относительные пути из config.sh считаются от корня репозитория.
abspath() { case $1 in /*) echo "$1" ;; *) echo "$ROOT/$1" ;; esac; }
RAW=$(abspath "$RAW"); BIRAGZANG=$(abspath "$BIRAGZANG"); DB=$(abspath "$DB"); OUT=$(abspath "$OUT")
ANALYSIS=$OUT/analysis                       # рабочий каталог скриптов анализа
PIPE=$ROOT/workflow/pipeline; AN=$ROOT/workflow/analysis
PY=$ENVS/mvtools/bin/python                  # Python 3.11 из окружения mvtools
export OUT ANALYSIS BIRAGZANG                # их читает workflow/analysis/paths.py
mkdir -p "$DB" "$OUT/logs" "$ANALYSIS/logs" "$ANALYSIS/report"

MV="MV-4122_S12 MV-4125_S18 MV-4128_S29"                                   # фагомы
PILOT="$MV NC1_S64"                                                        # фагомы и отрицательный контроль
COMPANIONS="4122-0-22_S11 4122-0-44_S10 4125-0-22_S17 4125-0-44_S16 4128-0-22_S28 4128-0-44_S27 4128_S30"
TIB_CELLULAR="4122-0-22_S11 4122-0-44_S10 4125-0-22_S17 4125-0-44_S16"   # клеточные фракции Тиб (собираются)
ALL_LIBS="$MV 4122-0-22_S11 4122-0-44_S10 4125-0-22_S17 4125-0-44_S16 4128-0-22_S28 4128-0-44_S27 4128_S30 NC1_S64"

log() { echo "[$(date '+%F %T')] $*" >&2; }
