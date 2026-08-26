#!/usr/bin/env bash
set -euo pipefail

# HPC wrapper for the model-2.1 equilibrium-and-precision second pilot.
# Override concurrency with TROPHOSOME_SECOND_PILOT_JOBS if needed.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPOSITORY="$(cd "$SCRIPT_DIR/../.." && pwd)"
MAMBA_ENVIRONMENT="${TROPHOSOME_MAMBA_ENV:-trophosome}"
PILOT_JOBS="${TROPHOSOME_SECOND_PILOT_JOBS:-8}"

if ! MAMBA_EXECUTABLE="$(command -v mamba)"; then
  echo "mamba was not found in PATH." >&2
  exit 2
fi

# Avoid `mamba run`, which generates an invalid `exec --` wrapper on the HPC.
if ! MAMBA_SHELL_HOOK="$("$MAMBA_EXECUTABLE" shell hook -s bash)"; then
  echo "Could not initialise mamba for a non-interactive Bash shell." >&2
  exit 2
fi
eval "$MAMBA_SHELL_HOOK"

# Conda activation scripts can read optional backup variables while nounset is
# active. Suspend it only for activation and restore strict checking afterward.
set +u
MAMBA_ACTIVATION_STATUS=0
mamba activate "$MAMBA_ENVIRONMENT" || MAMBA_ACTIVATION_STATUS=$?
set -u
if ((MAMBA_ACTIVATION_STATUS != 0)); then
  echo "Could not activate mamba environment: $MAMBA_ENVIRONMENT" >&2
  exit 2
fi

if ! PYTHON_EXECUTABLE="$(command -v python)"; then
  echo "Python was not found after activating: $MAMBA_ENVIRONMENT" >&2
  exit 2
fi

exec "$PYTHON_EXECUTABLE" \
  "$REPOSITORY/scripts/run_phase1_second_pilot.py" \
  --repository "$REPOSITORY" \
  --jobs "$PILOT_JOBS" \
  "$@"
