#!/usr/bin/env bash
set -euo pipefail
ASMSV_NEXT_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
export PYTHONPATH="$ASMSV_NEXT_ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
exec "${ASMSV_CALL_PY:-python3}" -B -m asmsv_next.cli "$@"
