#!/usr/bin/env bash

# Shared by the HPC launchers. This file is sourced, so it deliberately does
# not change the caller's strict-shell settings permanently.
trophosome_select_python() {
  local requested_environment="$1"
  local active_environment="${CONDA_DEFAULT_ENV:-${MAMBA_DEFAULT_ENV:-}}"
  local active_prefix_name=""
  local mamba_executable=""
  local mamba_shell_hook=""
  local activation_status=0

  if [[ -n "${CONDA_PREFIX:-}" ]]; then
    active_prefix_name="${CONDA_PREFIX##*/}"
  fi

  # An activated environment already supplies the correct Python. Some HPC
  # installations do not keep the base-environment `mamba` executable on PATH
  # after activation, so do not require mamba in this case.
  if [[ "$active_environment" == "$requested_environment" \
        || "$active_prefix_name" == "$requested_environment" ]]; then
    if [[ -n "${CONDA_PREFIX:-}" && -x "$CONDA_PREFIX/bin/python" ]]; then
      PYTHON_EXECUTABLE="$CONDA_PREFIX/bin/python"
    elif PYTHON_EXECUTABLE="$(command -v python)"; then
      :
    else
      echo "Python was not found in the active environment: $requested_environment" >&2
      return 2
    fi
    return 0
  fi

  if mamba_executable="$(command -v mamba)"; then
    :
  elif [[ -n "${MAMBA_EXE:-}" && -x "$MAMBA_EXE" ]]; then
    mamba_executable="$MAMBA_EXE"
  else
    echo "The '$requested_environment' environment is not active and mamba was not found." >&2
    echo "Activate the environment first, or make mamba available in PATH." >&2
    return 2
  fi

  # Avoid `mamba run`, which generates an invalid `exec --` wrapper on the
  # tested HPC installation. Initialise this shell and activate directly.
  if ! mamba_shell_hook="$("$mamba_executable" shell hook -s bash)"; then
    echo "Could not initialise mamba for a non-interactive Bash shell." >&2
    return 2
  fi
  eval "$mamba_shell_hook"

  # Third-party compiler hooks can inspect optional unset variables. Suspend
  # nounset only while those hooks run, then restore it immediately.
  set +u
  mamba activate "$requested_environment" || activation_status=$?
  set -u
  if ((activation_status != 0)); then
    echo "Could not activate mamba environment: $requested_environment" >&2
    return 2
  fi

  if ! PYTHON_EXECUTABLE="$(command -v python)"; then
    echo "Python was not found after activating: $requested_environment" >&2
    return 2
  fi
}
