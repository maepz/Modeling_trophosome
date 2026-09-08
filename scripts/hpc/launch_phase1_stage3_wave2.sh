#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPOSITORY="$(cd "$SCRIPT_DIR/../.." && pwd)"
MAMBA_ENVIRONMENT="${TROPHOSOME_MAMBA_ENV:-trophosome}"
WAVE_JOBS="${TROPHOSOME_STAGE3_WAVE2_JOBS:-8}"

source "$SCRIPT_DIR/_activate_environment.sh"
source "$SCRIPT_DIR/_completion_email.sh"
if ! trophosome_select_python "$MAMBA_ENVIRONMENT"; then
  exit 2
fi

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
COMMAND=(
  "$PYTHON_EXECUTABLE" "$REPOSITORY/scripts/run_phase1_stage3_wave2.py"
  --repository "$REPOSITORY" --jobs "$WAVE_JOBS" "$@"
)
if trophosome_invocation_is_job "$@"; then
  trophosome_run_with_completion_notification \
    "Phase 1 Stage 3 Wave 2" "$REPOSITORY" "${COMMAND[@]}"
else
  exec "${COMMAND[@]}"
fi
