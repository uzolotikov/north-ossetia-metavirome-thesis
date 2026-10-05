#!/bin/bash
# Шаг 5. Отбор вирусных кандидатов (правила — filter_candidates_README.md) и повторный CheckV по строгому набору.
set -euo pipefail
source "$(dirname "$0")/config.sh"; cd "$P"
F=viral/filtered_no_NC/filtering_no_NC
$PY "$(dirname "$0")/filter_candidates.py" checkv/no_NC viral/no_NC $F
$ENVS/checkv-env/bin/checkv end_to_end $F/strict_supported.fna checkv/filtered_no_NC/strict_supported_pass2 \
  -d $P/databases/checkv-db-v1.5 -t 4 > logs/strict_supported.checkv_pass2.log 2>&1
