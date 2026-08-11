#!/usr/bin/env bash
set -uo pipefail

# Read-only environment probe for the trophosome HPC workflow.
#
# The only file this script creates is the requested report. It does not make
# network requests, install packages, activate environments, or change shell
# configuration. Run it from the shell environment normally used to launch
# model jobs so that the report describes the relevant mamba/Python setup.

SCRIPT_VERSION="1.0.0"
OUTPUT_PATH=""
PROJECT_ROOT=""
QUIET=false

usage() {
    cat <<'EOF'
Usage:
  bash scripts/hpc/probe_hpc_environment.sh [options]

Options:
  --output PATH        Report path. The default is a timestamped file in the
                       current directory.
  --project-root PATH  HPC project root to inspect. The default is
                       $CRF_PROJECT_ROOT or $HOME/data/CRF_project.
  --quiet              Write only to the report instead of also printing it.
  -h, --help           Show this help message.

Example:
  bash scripts/hpc/probe_hpc_environment.sh \
    --project-root "$HOME/data/CRF_project"

Run the script after entering the shell or mamba environment normally used for
jobs. It is safe to run on a login node and does not start a simulation.
EOF
}

fail() {
    printf 'error: %s\n' "$*" >&2
    exit 2
}

while (($# > 0)); do
    case "$1" in
        --output)
            (($# >= 2)) || fail "--output requires a path"
            OUTPUT_PATH=$2
            shift 2
            ;;
        --project-root)
            (($# >= 2)) || fail "--project-root requires a path"
            PROJECT_ROOT=$2
            shift 2
            ;;
        --quiet)
            QUIET=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            fail "unknown option: $1"
            ;;
    esac
done

if [[ -z "$PROJECT_ROOT" ]]; then
    if [[ -n "${CRF_PROJECT_ROOT:-}" ]]; then
        PROJECT_ROOT=$CRF_PROJECT_ROOT
    elif [[ -n "${HOME:-}" ]]; then
        PROJECT_ROOT="${HOME}/data/CRF_project"
    else
        PROJECT_ROOT="$(pwd -P)"
    fi
fi

if [[ -z "$OUTPUT_PATH" ]]; then
    REPORT_TIMESTAMP=$(date -u '+%Y%m%dT%H%M%SZ' 2>/dev/null || date '+%Y%m%dT%H%M%S')
    OUTPUT_PATH="hpc_environment_probe_${REPORT_TIMESTAMP}.txt"
fi

OUTPUT_PARENT=$(dirname -- "$OUTPUT_PATH")
[[ -d "$OUTPUT_PARENT" ]] || fail "output directory does not exist: $OUTPUT_PARENT"
[[ ! -e "$OUTPUT_PATH" ]] || fail "refusing to overwrite existing report: $OUTPUT_PATH"

if [[ "$QUIET" == true ]]; then
    exec >"$OUTPUT_PATH" 2>&1
else
    exec > >(tee "$OUTPUT_PATH") 2>&1
fi

section() {
    printf '\n================================================================================\n'
    printf '%s\n' "$1"
    printf '================================================================================\n'
}

print_command() {
    printf '$'
    printf ' %q' "$@"
    printf '\n'
}

run_command() {
    local label=$1
    shift
    printf '\n--- %s ---\n' "$label"
    print_command "$@"
    "$@" 2>&1
    local command_status=$?
    printf '[exit_status=%d]\n' "$command_status"
    return 0
}

run_shell() {
    local label=$1
    local shell_code=$2
    printf '\n--- %s ---\n' "$label"
    printf '$ %s\n' "$shell_code"
    bash -o pipefail -c "$shell_code" 2>&1
    local command_status=$?
    printf '[exit_status=%d]\n' "$command_status"
    return 0
}

sanitize_stream() {
    # Remove URL user information and common query-string credentials while
    # preserving channel hostnames and paths needed for environment diagnosis.
    sed -E \
        -e 's#(https?://)[^/@[:space:]]+@#\1[REDACTED]@#g' \
        -e 's#([?&](token|access_token|auth|key)=)[^&[:space:]]+#\1[REDACTED]#g'
}

run_sanitized() {
    local label=$1
    shift
    printf '\n--- %s ---\n' "$label"
    print_command "$@"
    "$@" 2>&1 | sanitize_stream
    local command_status=${PIPESTATUS[0]}
    printf '[exit_status=%d]\n' "$command_status"
    return 0
}

command_path() {
    if command -v "$1" >/dev/null 2>&1; then
        command -v "$1"
    else
        printf 'NOT_FOUND'
    fi
}

read_if_available() {
    local label=$1
    local file_path=$2
    printf '\n--- %s ---\n' "$label"
    printf 'file: %s\n' "$file_path"
    if [[ -r "$file_path" ]]; then
        sed -n '1,240p' "$file_path"
        printf '[read_status=0]\n'
    else
        printf 'NOT_READABLE_OR_NOT_PRESENT\n'
        printf '[read_status=1]\n'
    fi
}

section "Probe metadata"
printf 'probe_version: %s\n' "$SCRIPT_VERSION"
printf 'started_utc: %s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ' 2>/dev/null || date)"
printf 'hostname: %s\n' "$(hostname 2>/dev/null || printf 'UNKNOWN')"
printf 'working_directory: %s\n' "$(pwd -P)"
printf 'project_root_requested: %s\n' "$PROJECT_ROOT"
printf 'report_path: %s\n' "$OUTPUT_PATH"
printf 'network_operations_performed: false\n'
printf 'software_installation_performed: false\n'

section "Operating system and shell"
run_command "Kernel and architecture" uname -a
if [[ -r /etc/os-release ]]; then
    read_if_available "Operating-system release" /etc/os-release
elif command -v sw_vers >/dev/null 2>&1; then
    run_command "Operating-system release" sw_vers
else
    printf '\nOperating-system release information was not found.\n'
fi
run_command "Bash version" bash --version
if command -v locale >/dev/null 2>&1; then
    run_command "Locale" locale
fi
run_shell "Current resource limits" "ulimit -a"
if command -v getconf >/dev/null 2>&1; then
    run_command "Configured online processors" getconf _NPROCESSORS_ONLN
    run_command "Memory page size" getconf PAGE_SIZE
    run_command "Clock ticks per second" getconf CLK_TCK
fi

section "CPU information"
if command -v lscpu >/dev/null 2>&1; then
    run_command "Detailed CPU information" lscpu
fi
if command -v nproc >/dev/null 2>&1; then
    run_command "Available processor count" nproc --all
fi
if command -v sysctl >/dev/null 2>&1; then
    run_command "macOS logical CPU count, when supported" sysctl -n hw.logicalcpu
    run_command "macOS physical CPU count, when supported" sysctl -n hw.physicalcpu
fi
if command -v taskset >/dev/null 2>&1; then
    run_command "Probe-shell CPU affinity" taskset -pc "$$"
fi

section "Memory and process limits"
if command -v free >/dev/null 2>&1; then
    run_command "Memory summary" free -h
fi
read_if_available "Linux memory information" /proc/meminfo
if command -v vm_stat >/dev/null 2>&1; then
    run_command "macOS virtual-memory summary" vm_stat
fi
if command -v ps >/dev/null 2>&1; then
    run_command \
        "Probe process fields" \
        ps -o pid=,ppid=,rss=,vsz=,etime=,command= -p "$$"
fi

section "Control groups and memory-accounting interfaces"
read_if_available "Process control-group membership" /proc/self/cgroup
for cgroup_file in \
    /sys/fs/cgroup/memory.max \
    /sys/fs/cgroup/memory.current \
    /sys/fs/cgroup/memory.peak \
    /sys/fs/cgroup/memory.swap.max \
    /sys/fs/cgroup/memory/memory.limit_in_bytes \
    /sys/fs/cgroup/memory/memory.usage_in_bytes \
    /sys/fs/cgroup/memory/memory.max_usage_in_bytes
do
    read_if_available "Cgroup value" "$cgroup_file"
done
if command -v mount >/dev/null 2>&1; then
    run_shell "Mounted cgroup filesystems" "mount | grep -E 'cgroup|cgroup2' || true"
fi

section "Project filesystem"
if [[ -e "$PROJECT_ROOT" ]]; then
    printf 'project_root_status: EXISTS\n'
    run_command "Project-root disk space" df -h "$PROJECT_ROOT"
    run_command "Project-root inode space" df -i "$PROJECT_ROOT"
    if df -T "$PROJECT_ROOT" >/dev/null 2>&1; then
        run_command "Project-root filesystem type" df -T "$PROJECT_ROOT"
    fi
else
    printf 'project_root_status: NOT_FOUND\n'
fi
for expected_directory in work scratch data; do
    expected_path="${PROJECT_ROOT}/${expected_directory}"
    if [[ -d "$expected_path" ]]; then
        printf '%s: PRESENT\n' "$expected_path"
    else
        printf '%s: ABSENT\n' "$expected_path"
    fi
done

section "Background-job and monitoring commands"
for utility in \
    nohup setsid tmux screen flock ps pgrep pkill nice renice ionice timeout \
    time rsync git mamba micromamba conda python python3 pip pip3
do
    printf '%-16s %s\n' "$utility" "$(command_path "$utility")"
done
if [[ -x /usr/bin/time ]]; then
    run_command "System time implementation" /usr/bin/time --version
    run_command "Verbose resource-accounting sample" /usr/bin/time -v bash -c ':'
else
    printf '\n/usr/bin/time: NOT_EXECUTABLE_OR_NOT_PRESENT\n'
fi
if command -v git >/dev/null 2>&1; then
    run_command "Git version" git --version
fi

section "Mamba and Conda"
ENVIRONMENT_MANAGER=""
for candidate in mamba micromamba conda; do
    if command -v "$candidate" >/dev/null 2>&1; then
        ENVIRONMENT_MANAGER=$candidate
        break
    fi
done

if [[ -n "$ENVIRONMENT_MANAGER" ]]; then
    printf 'environment_manager_selected: %s\n' "$ENVIRONMENT_MANAGER"
    run_command \
        "Environment-manager version" \
        "$ENVIRONMENT_MANAGER" --version
    run_sanitized \
        "Environment-manager information" \
        "$ENVIRONMENT_MANAGER" info
    run_sanitized \
        "Available environments" \
        "$ENVIRONMENT_MANAGER" env list
    run_sanitized \
        "Installed packages in the active environment" \
        "$ENVIRONMENT_MANAGER" list
else
    printf 'environment_manager_selected: NOT_FOUND\n'
fi
printf 'CONDA_DEFAULT_ENV: %s\n' "${CONDA_DEFAULT_ENV:-UNSET}"
printf 'CONDA_PREFIX: %s\n' "${CONDA_PREFIX:-UNSET}"

section "Python runtime"
PYTHON_COMMAND=""
for candidate in python python3; do
    if command -v "$candidate" >/dev/null 2>&1; then
        PYTHON_COMMAND=$candidate
        break
    fi
done

if [[ -n "$PYTHON_COMMAND" ]]; then
    run_command "Python version" "$PYTHON_COMMAND" --version
    run_command "Pip version" "$PYTHON_COMMAND" -m pip --version
    printf '\n--- Python platform and relevant package versions ---\n'
    print_command "$PYTHON_COMMAND" -
    "$PYTHON_COMMAND" - <<'PY'
from __future__ import annotations

import importlib.metadata
import multiprocessing
import os
import platform
import sys

print(f"executable: {sys.executable}")
print(f"version: {sys.version.replace(chr(10), ' ')}")
print(f"prefix: {sys.prefix}")
print(f"base_prefix: {sys.base_prefix}")
print(f"platform: {platform.platform()}")
print(f"machine: {platform.machine()}")
print(f"processor: {platform.processor()}")
print(f"python_implementation: {platform.python_implementation()}")
print(f"os_cpu_count: {os.cpu_count()}")
print(f"multiprocessing_start_method: {multiprocessing.get_start_method(allow_none=True)}")

for variable in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
):
    print(f"{variable}: {os.environ.get(variable, 'UNSET')}")

for distribution in (
    "numpy",
    "networkx",
    "psutil",
    "trophosome-model",
    "pytest",
    "ruff",
):
    try:
        version = importlib.metadata.version(distribution)
    except importlib.metadata.PackageNotFoundError:
        version = "NOT_INSTALLED"
    print(f"package.{distribution}: {version}")

try:
    import numpy as np
except ImportError:
    print("numpy_configuration: NOT_AVAILABLE")
else:
    print("numpy_configuration:")
    np.show_config()
PY
    printf '[exit_status=%d]\n' "$?"
else
    printf 'python_runtime: NOT_FOUND\n'
fi

section "Probe completion"
printf 'completed_utc: %s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ' 2>/dev/null || date)"
printf 'report_path: %s\n' "$OUTPUT_PATH"
printf 'result: COMPLETE\n'
printf '\nTransfer this report back to the Mac for review. Inspect it before sharing\n'
printf 'outside the research project because it contains hostnames and local paths.\n'
