#!/usr/bin/env bash
set -euo pipefail

# Thin HPC wrapper for the fixed-pool model-2.1 first pilot.  Override the
# defaults with TROPHOSOME_MAMBA_ENV and TROPHOSOME_PILOT_JOBS if needed.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPOSITORY="$(cd "$SCRIPT_DIR/../.." && pwd)"
MAMBA_ENVIRONMENT="${TROPHOSOME_MAMBA_ENV:-trophosome}"
PILOT_JOBS="${TROPHOSOME_PILOT_JOBS:-8}"

if ! MAMBA_EXECUTABLE="$(command -v mamba)"; then
  echo "mamba was not found in PATH." >&2
  exit 2
fi

# On some HPC installations, including the tested mamba 2.1 and 2.2 setups,
# `mamba run` generates a temporary script containing `exec --`, which Bash
# rejects.  Activate the named environment in this non-interactive shell
# instead, then invoke its Python directly.
if ! MAMBA_SHELL_HOOK="$("$MAMBA_EXECUTABLE" shell hook -s bash)"; then
  echo "Could not initialise mamba for a non-interactive Bash shell." >&2
  exit 2
fi
eval "$MAMBA_SHELL_HOOK"

if ! mamba activate "$MAMBA_ENVIRONMENT"; then
  echo "Could not activate mamba environment: $MAMBA_ENVIRONMENT" >&2
  exit 2
fi

if ! PYTHON_EXECUTABLE="$(command -v python)"; then
  echo "Python was not found after activating: $MAMBA_ENVIRONMENT" >&2
  exit 2
fi

exec "$PYTHON_EXECUTABLE" \
  "$REPOSITORY/scripts/run_phase1_first_pilot_v2_1.py" \
  --repository "$REPOSITORY" \
  --jobs "$PILOT_JOBS" \
  "$@"
