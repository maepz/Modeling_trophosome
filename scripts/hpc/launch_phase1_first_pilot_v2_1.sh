#!/usr/bin/env bash
set -euo pipefail

# Thin HPC wrapper for the fixed-pool model-2.1 first pilot.  Override the
# defaults with TROPHOSOME_MAMBA_ENV and TROPHOSOME_PILOT_JOBS if needed.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPOSITORY="$(cd "$SCRIPT_DIR/../.." && pwd)"
MAMBA_ENVIRONMENT="${TROPHOSOME_MAMBA_ENV:-trophosome}"
PILOT_JOBS="${TROPHOSOME_PILOT_JOBS:-8}"

if ! command -v mamba >/dev/null 2>&1; then
  echo "mamba was not found in PATH." >&2
  exit 2
fi

exec mamba run --no-capture-output -n "$MAMBA_ENVIRONMENT" \
  python "$REPOSITORY/scripts/run_phase1_first_pilot_v2_1.py" \
  --repository "$REPOSITORY" \
  --jobs "$PILOT_JOBS" \
  "$@"
