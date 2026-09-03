#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPOSITORY="$(cd "$SCRIPT_DIR/../.." && pwd)"
MAMBA_ENVIRONMENT="${TROPHOSOME_MAMBA_ENV:-trophosome}"
WAVE_JOBS="${TROPHOSOME_STAGE3_JOBS:-8}"

source "$SCRIPT_DIR/_activate_environment.sh"
if ! trophosome_select_python "$MAMBA_ENVIRONMENT"; then
  exit 2
fi

# Limit BLAS to one thread inside each worker to avoid CPU oversubscription.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
exec "$PYTHON_EXECUTABLE" "$REPOSITORY/scripts/run_phase1_stage3_wave1.py" \
  --repository "$REPOSITORY" --jobs "$WAVE_JOBS" "$@"
