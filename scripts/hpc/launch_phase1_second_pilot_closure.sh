#!/usr/bin/env bash
set -euo pipefail

# Launch only the eight matched seed blocks added for Stage 2 closure.
# The shared launcher still requires all 120 populations before rebuilding the
# combined Stage 2 report, so the 72 original populations must remain in scratch.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

for argument in "$@"; do
  case "$argument" in
    --seed-block|--seed-block=*)
      echo "The closure launcher fixes seed blocks to sb0013--sb0020." >&2
      exit 2
      ;;
  esac
done

exec "$SCRIPT_DIR/launch_phase1_second_pilot.sh" \
  --seed-block sb0013 \
  --seed-block sb0014 \
  --seed-block sb0015 \
  --seed-block sb0016 \
  --seed-block sb0017 \
  --seed-block sb0018 \
  --seed-block sb0019 \
  --seed-block sb0020 \
  "$@"
