#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPOSITORY="$(cd "$SCRIPT_DIR/../.." && pwd)"
MAMBA_ENVIRONMENT="${TROPHOSOME_MAMBA_ENV:-trophosome}"
WAVE_JOBS="${TROPHOSOME_STAGE3_JOBS:-8}"

if ! MAMBA_EXECUTABLE="$(command -v mamba)"; then
  echo "mamba was not found in PATH." >&2
  exit 2
fi
if ! MAMBA_SHELL_HOOK="$("$MAMBA_EXECUTABLE" shell hook -s bash)"; then
  echo "Could not initialise mamba for this shell." >&2
  exit 2
fi
eval "$MAMBA_SHELL_HOOK"
# Conda compiler activation hooks read optional unset backup variables.
set +u
MAMBA_ACTIVATION_STATUS=0
mamba activate "$MAMBA_ENVIRONMENT" || MAMBA_ACTIVATION_STATUS=$?
set -u
if ((MAMBA_ACTIVATION_STATUS != 0)); then
  echo "Could not activate mamba environment: $MAMBA_ENVIRONMENT" >&2
  exit 2
fi
if ! PYTHON_EXECUTABLE="$(command -v python)"; then
  echo "Python was not found after activation." >&2
  exit 2
fi

# Limit BLAS to one thread inside each worker to avoid CPU oversubscription.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
exec "$PYTHON_EXECUTABLE" "$REPOSITORY/scripts/run_phase1_stage3_wave1.py" \
  --repository "$REPOSITORY" --jobs "$WAVE_JOBS" "$@"
