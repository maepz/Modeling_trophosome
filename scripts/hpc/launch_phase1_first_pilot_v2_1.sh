#!/usr/bin/env bash
set -euo pipefail

# Thin HPC wrapper for the fixed-pool model-2.1 first pilot.  Override the
# defaults with TROPHOSOME_MAMBA_ENV and TROPHOSOME_PILOT_JOBS if needed.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPOSITORY="$(cd "$SCRIPT_DIR/../.." && pwd)"
MAMBA_ENVIRONMENT="${TROPHOSOME_MAMBA_ENV:-trophosome}"
PILOT_JOBS="${TROPHOSOME_PILOT_JOBS:-8}"

source "$SCRIPT_DIR/_activate_environment.sh"
source "$SCRIPT_DIR/_completion_email.sh"
if ! trophosome_select_python "$MAMBA_ENVIRONMENT"; then
  exit 2
fi

COMMAND=(
  "$PYTHON_EXECUTABLE"
  "$REPOSITORY/scripts/run_phase1_first_pilot_v2_1.py"
  --repository "$REPOSITORY"
  --jobs "$PILOT_JOBS"
  "$@"
)
if trophosome_invocation_is_job "$@"; then
  trophosome_run_with_completion_notification \
    "Phase 1 first pilot" "$REPOSITORY" "${COMMAND[@]}"
else
  exec "${COMMAND[@]}"
fi
