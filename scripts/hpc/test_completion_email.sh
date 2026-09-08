#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPOSITORY="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$SCRIPT_DIR/_completion_email.sh"

RECIPIENT="$(trophosome_notification_recipient "$REPOSITORY" || true)"
if [[ -z "$RECIPIENT" ]]; then
  echo "No notification address is configured." >&2
  echo "Set TROPHOSOME_NOTIFY_EMAIL, then run this test again." >&2
  exit 2
fi
if ! TRANSPORT="$(trophosome_email_transport)"; then
  echo "No mail client was found. Ask the HPC administrator about mailx or sendmail." >&2
  exit 2
fi

HOST="$(hostname 2>/dev/null || printf 'unknown')"
SUBJECT="trophosome email test from $HOST"
printf -v BODY \
  'This is a test of trophosome HPC completion notifications.\n\nIf you received this message, future long-running launchers can notify you using %s.' \
  "$TRANSPORT"
trophosome_send_completion_email "$RECIPIENT" "$SUBJECT" "$BODY"
echo "Test email submitted to $RECIPIENT using $TRANSPORT."
