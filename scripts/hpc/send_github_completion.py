#!/usr/bin/env python3
"""Request the GitHub workflow that reports an HPC job completion."""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

WORKFLOW_FILE = "hpc-completion-notification.yml"
DEFAULT_REF = "main"
API_VERSION = "2022-11-28"


def default_token_path() -> Path:
    configured = os.environ.get("TROPHOSOME_GITHUB_TOKEN_FILE")
    if configured:
        return Path(configured).expanduser()
    config_home = os.environ.get("XDG_CONFIG_HOME")
    root = Path(config_home).expanduser() if config_home else Path.home() / ".config"
    return root / "trophosome" / "github-token"


def read_token() -> str:
    token = os.environ.get("TROPHOSOME_GITHUB_TOKEN", "").strip()
    if token:
        return token
    path = default_token_path()
    try:
        if path.stat().st_mode & 0o077:
            raise RuntimeError(
                f"GitHub notification token is readable by other users: {path}; "
                f"run chmod 600 {path}"
            )
        token = path.read_text(encoding="utf-8").strip()
    except RuntimeError:
        raise
    except OSError as exc:
        raise RuntimeError(
            f"GitHub notification token is unavailable at {path}"
        ) from exc
    if not token:
        raise RuntimeError(f"GitHub notification token is empty: {path}")
    return token


def repository_slug(repository: Path) -> str:
    configured = os.environ.get("TROPHOSOME_GITHUB_REPOSITORY", "").strip()
    if configured:
        candidate = configured
    else:
        try:
            candidate = subprocess.run(
                ["git", "-C", str(repository), "remote", "get-url", "origin"],
                check=True,
                capture_output=True,
                text=True,
            ).stdout.strip()
        except (OSError, subprocess.CalledProcessError) as exc:
            raise RuntimeError("cannot determine the GitHub origin repository") from exc

    patterns = (
        r"https://github\.com/(?P<slug>[^/\s]+/[^/\s]+?)(?:\.git)?$",
        r"git@github\.com:(?P<slug>[^/\s]+/[^/\s]+?)(?:\.git)?$",
        r"ssh://git@github\.com/(?P<slug>[^/\s]+/[^/\s]+?)(?:\.git)?$",
        r"(?P<slug>[^/\s]+/[^/\s]+)$",
    )
    for pattern in patterns:
        match = re.fullmatch(pattern, candidate)
        if match:
            return match.group("slug")
    raise RuntimeError(f"origin is not a supported GitHub repository: {candidate}")


def dispatch_url(slug: str) -> str:
    owner, name = slug.split("/", maxsplit=1)
    workflow = urllib.parse.quote(WORKFLOW_FILE, safe="")
    return (
        "https://api.github.com/repos/"
        f"{urllib.parse.quote(owner, safe='')}/{urllib.parse.quote(name, safe='')}"
        f"/actions/workflows/{workflow}/dispatches"
    )


def dispatch(repository: Path, inputs: dict[str, str]) -> str:
    slug = repository_slug(repository)
    payload = json.dumps(
        {
            "ref": os.environ.get("TROPHOSOME_GITHUB_REF", DEFAULT_REF),
            "inputs": inputs,
        }
    ).encode("utf-8")
    request = urllib.request.Request(
        dispatch_url(slug),
        data=payload,
        method="POST",
        headers={
            "Accept": "application/vnd.github+json",
            "Authorization": f"Bearer {read_token()}",
            "Content-Type": "application/json",
            "User-Agent": "trophosome-hpc-notifier",
            "X-GitHub-Api-Version": API_VERSION,
        },
    )
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            status = response.status
    except urllib.error.HTTPError as exc:
        try:
            detail = exc.read().decode("utf-8", errors="replace")
        except OSError:
            detail = ""
        suffix = f": {detail[:500]}" if detail else ""
        raise RuntimeError(f"GitHub returned HTTP {exc.code}{suffix}") from exc
    except urllib.error.URLError as exc:
        raise RuntimeError(f"cannot reach GitHub: {exc.reason}") from exc
    if status != 204:
        raise RuntimeError(f"GitHub returned unexpected HTTP status {status}")
    return slug


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Ask GitHub to create a success/failure notification for an HPC job. "
            "The token is read without printing it."
        )
    )
    parser.add_argument("--repository", type=Path, required=True)
    parser.add_argument(
        "--check", action="store_true", help="validate local setup only"
    )
    parser.add_argument(
        "--test", action="store_true", help="send a harmless test notice"
    )
    parser.add_argument("--label")
    parser.add_argument("--outcome", choices=("SUCCESS", "FAILED", "TEST"))
    parser.add_argument("--exit-status", type=int)
    parser.add_argument("--host")
    parser.add_argument("--started-utc")
    parser.add_argument("--finished-utc")
    parser.add_argument("--elapsed-seconds", type=int)
    parser.add_argument("--revision")
    parser.add_argument("--command")
    return parser


def main() -> int:
    parser = _parser()
    args = parser.parse_args()
    repository = args.repository.resolve()
    if args.check and args.test:
        parser.error("--check and --test cannot be combined")
    try:
        slug = repository_slug(repository)
        read_token()
        workflow = repository / ".github" / "workflows" / WORKFLOW_FILE
        if not workflow.is_file():
            raise RuntimeError(f"GitHub notification workflow is missing: {workflow}")
        if args.check:
            print(f"GitHub notification configuration is present for {slug}.")
            return 0

        if args.test:
            inputs = {
                "label": "HPC notification test",
                "outcome": "TEST",
                "exit_status": "0",
                "host": "setup-test",
                "started_utc": "not applicable",
                "finished_utc": "not applicable",
                "elapsed_seconds": "0",
                "revision": "not applicable",
                "command": "notification setup test",
            }
        else:
            required = {
                "label": args.label,
                "outcome": args.outcome,
                "exit_status": args.exit_status,
                "host": args.host,
                "started_utc": args.started_utc,
                "finished_utc": args.finished_utc,
                "elapsed_seconds": args.elapsed_seconds,
                "revision": args.revision,
                "command": args.command,
            }
            missing = [name for name, value in required.items() if value is None]
            if missing:
                parser.error("missing job fields: " + ", ".join(missing))
            inputs = {name: str(value) for name, value in required.items()}
        accepted = dispatch(repository, inputs)
    except RuntimeError as exc:
        print(f"GitHub notification was not sent: {exc}", file=sys.stderr)
        return 2
    print(
        f"GitHub accepted the notification for {accepted}. "
        "The workflow will create the notice shortly."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
