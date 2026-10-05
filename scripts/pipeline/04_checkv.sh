#!/bin/bash
# Шаг 4. Качество и обрезка провирусных краёв (CheckV 1.1.1, база 1.5).
set -euo pipefail
source "$(dirname "$0")/config.sh"; cd "$P"
for s in $MV; do
  $ENVS/checkv-env/bin/checkv end_to_end viral/no_NC/$s/final-viral-combined.fa checkv/no_NC/$s \
    -d $P/databases/checkv-db-v1.5 -t 4 > logs/$s.checkv.log 2>&1
done
