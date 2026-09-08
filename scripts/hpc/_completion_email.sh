#!/usr/bin/env bash

# Shared completion-email support for long-running HPC launchers. This file is
# sourced by launch scripts; notification failures never change a model job's
# exit status.

trophosome_invocation_is_job() {
  local argument
  for argument in "$@"; do
    case "$argument" in
      --prepare-only|--dry-run|--check-smoke|--assess-only|--summarize-only|--report-only)
        return 1
        ;;
    esac
  done
  return 0
}

trophosome_notification_recipient() {
  local repository="$1"
  local configured="${TROPHOSOME_NOTIFY_EMAIL-}"

  case "$configured" in
    [Oo][Ff][Ff]|[Nn][Oo][Nn][Ee]|[Dd][Ii][Ss][Aa][Bb][Ll][Ee][Dd])
      return 1
      ;;
  esac
  if [[ -n "$configured" ]]; then
    printf '%s\n' "$configured"
    return 0
  fi
  git -C "$repository" config --get user.email 2>/dev/null
}

trophosome_email_transport() {
  local candidate
  for candidate in mailx mail sendmail; do
    if command -v "$candidate" >/dev/null 2>&1; then
      printf '%s\n' "$candidate"
      return 0
    fi
  done
  return 1
}

trophosome_send_completion_email() {
  local recipient="$1"
  local subject="$2"
  local body="$3"
  local transport=""

  if [[ "$recipient" == *$'\n'* || "$recipient" == *$'\r'* ]]; then
    echo "Completion email was not sent: invalid recipient." >&2
    return 2
  fi
  if ! transport="$(trophosome_email_transport)"; then
    echo "Completion email was not sent: mailx, mail and sendmail are unavailable." >&2
    return 2
  fi
  case "$transport" in
    mailx|mail)
      printf '%s\n' "$body" | "$transport" -s "$subject" "$recipient"
      ;;
    sendmail)
      {
        printf 'To: %s\n' "$recipient"
        printf 'Subject: %s\n' "$subject"
        printf 'Content-Type: text/plain; charset=UTF-8\n'
        printf '\n%s\n' "$body"
      } | sendmail "$recipient"
      ;;
  esac
}

trophosome_run_with_completion_email() {
  local label="$1"
  local repository="$2"
  shift 2
  local recipient=""
  local transport=""
  local started_utc=""
  local finished_utc=""
  local started_epoch=0
  local finished_epoch=0
  local elapsed=0
  local status=0
  local outcome="SUCCESS"
  local command_text=""
  local body=""
  local subject=""
  local host=""
  local revision=""

  recipient="$(trophosome_notification_recipient "$repository" || true)"
  if [[ -n "$recipient" ]]; then
    if transport="$(trophosome_email_transport)"; then
      echo "A completion email will be sent to $recipient using $transport." >&2
    else
      echo "Warning: completion email requested for $recipient, but no mail client is available." >&2
      echo "Run scripts/hpc/test_completion_email.sh before a long job." >&2
      recipient=""
    fi
  fi

  started_utc="$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
  started_epoch="$(date '+%s')"
  if "$@"; then
    status=0
  else
    status=$?
  fi
  finished_utc="$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
  finished_epoch="$(date '+%s')"
  elapsed=$((finished_epoch - started_epoch))
  if ((status != 0)); then
    outcome="FAILED"
  fi

  if [[ -n "$recipient" ]]; then
    printf -v command_text '%q ' "$@"
    host="$(hostname 2>/dev/null || printf 'unknown')"
    revision="$(git -C "$repository" rev-parse --short HEAD 2>/dev/null || printf 'unknown')"
    subject="trophosome $outcome: $label"
    printf -v body '%s\n\nStatus: %s (exit %d)\nHost: %s\nStarted UTC: %s\nFinished UTC: %s\nElapsed seconds: %d\nGit revision: %s\nCommand: %s\n' \
      "$label finished." "$outcome" "$status" "$host" "$started_utc" \
      "$finished_utc" "$elapsed" "$revision" "$command_text"
    if ! trophosome_send_completion_email "$recipient" "$subject" "$body"; then
      echo "Warning: the job finished, but its completion email could not be sent." >&2
    fi
  fi
  return "$status"
}
