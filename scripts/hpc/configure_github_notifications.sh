#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPOSITORY="$(cd "$SCRIPT_DIR/../.." && pwd)"
CONFIG_ROOT="${XDG_CONFIG_HOME:-$HOME/.config}"
TOKEN_FILE="${TROPHOSOME_GITHUB_TOKEN_FILE:-$CONFIG_ROOT/trophosome/github-token}"

if [[ ! -t 0 ]]; then
  echo "Run this setup interactively so the token can be entered without echoing it." >&2
  exit 2
fi

if command -v python >/dev/null 2>&1; then
  PYTHON_COMMAND="$(command -v python)"
elif command -v python3 >/dev/null 2>&1; then
  PYTHON_COMMAND="$(command -v python3)"
else
  echo "Python is required to configure GitHub notifications." >&2
  exit 2
fi

echo "Paste the fine-grained GitHub token created for maepz/Modeling_trophosome."
printf "The token will be stored privately in %s\n" "$TOKEN_FILE"
printf "Token: "
IFS= read -r -s TOKEN
printf "\n"
if [[ -z "$TOKEN" ]]; then
  echo "No token was entered; nothing was changed." >&2
  exit 2
fi

umask 077
mkdir -p "$(dirname "$TOKEN_FILE")"
printf '%s\n' "$TOKEN" > "$TOKEN_FILE"
unset TOKEN
chmod 600 "$TOKEN_FILE"

echo "Saved the token with owner-only permissions. Sending a test notice..."
"$PYTHON_COMMAND" "$SCRIPT_DIR/send_github_completion.py" \
  --repository "$REPOSITORY" --test
