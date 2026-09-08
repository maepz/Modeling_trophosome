from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from scripts.hpc import send_github_completion as github_notification

REPOSITORY = Path(__file__).resolve().parents[1]


class _Response:
    status = 204

    def __enter__(self) -> _Response:
        return self

    def __exit__(self, *args: object) -> None:
        return None


class GitHubCompletionNotificationTests(unittest.TestCase):
    def test_every_phase_launcher_uses_the_shared_notification_wrapper(self) -> None:
        launchers = sorted((REPOSITORY / "scripts/hpc").glob("launch_phase*.sh"))
        self.assertGreaterEqual(len(launchers), 5)
        for launcher in launchers:
            with self.subTest(launcher=launcher.name):
                content = launcher.read_text(encoding="utf-8")
                direct = (
                    'source "$SCRIPT_DIR/_completion_email.sh"' in content
                    and "trophosome_run_with_completion_notification" in content
                )
                delegated = 'exec "$SCRIPT_DIR/launch_phase' in content
                self.assertTrue(
                    direct or delegated,
                    "a phase launcher must notify directly or delegate to a "
                    "notifying launcher",
                )

    def test_repository_slug_accepts_https_and_ssh_origins(self) -> None:
        cases = (
            "https://github.com/maepz/Modeling_trophosome.git",
            "git@github.com:maepz/Modeling_trophosome.git",
            "ssh://git@github.com/maepz/Modeling_trophosome.git",
        )
        for origin in cases:
            with (
                self.subTest(origin=origin),
                patch.dict(
                    os.environ,
                    {"TROPHOSOME_GITHUB_REPOSITORY": origin},
                    clear=False,
                ),
            ):
                self.assertEqual(
                    github_notification.repository_slug(REPOSITORY),
                    "maepz/Modeling_trophosome",
                )

    def test_dispatch_uses_workflow_api_without_exposing_token_in_payload(
        self,
    ) -> None:
        inputs = {
            "label": "Synthetic run",
            "outcome": "SUCCESS",
            "exit_status": "0",
            "host": "node01",
            "started_utc": "start",
            "finished_utc": "finish",
            "elapsed_seconds": "3",
            "revision": "abc1234",
            "command": "synthetic command",
        }
        with (
            patch.dict(
                os.environ,
                {
                    "TROPHOSOME_GITHUB_REPOSITORY": "maepz/Modeling_trophosome",
                    "TROPHOSOME_GITHUB_TOKEN": "secret-test-token",
                },
                clear=False,
            ),
            patch.object(
                github_notification.urllib.request,
                "urlopen",
                return_value=_Response(),
            ) as urlopen,
        ):
            slug = github_notification.dispatch(REPOSITORY, inputs)

        self.assertEqual(slug, "maepz/Modeling_trophosome")
        request = urlopen.call_args.args[0]
        self.assertTrue(
            request.full_url.endswith("/hpc-completion-notification.yml/dispatches")
        )
        payload = json.loads(request.data.decode("utf-8"))
        self.assertEqual(payload, {"ref": "main", "inputs": inputs})
        self.assertNotIn("secret-test-token", request.data.decode("utf-8"))

    def test_check_validates_local_configuration_without_network(self) -> None:
        with (
            patch.dict(
                os.environ,
                {
                    "TROPHOSOME_GITHUB_REPOSITORY": "maepz/Modeling_trophosome",
                    "TROPHOSOME_GITHUB_TOKEN": "secret-test-token",
                },
                clear=False,
            ),
            patch.object(
                sys,
                "argv",
                [
                    "send_github_completion.py",
                    "--repository",
                    str(REPOSITORY),
                    "--check",
                ],
            ),
            patch.object(github_notification, "dispatch") as dispatch,
        ):
            self.assertEqual(github_notification.main(), 0)
        dispatch.assert_not_called()

    def test_shell_wrapper_preserves_job_failure_when_github_is_used(self) -> None:
        helper = REPOSITORY / "scripts/hpc/_completion_email.sh"
        with tempfile.TemporaryDirectory() as temporary_directory:
            temporary = Path(temporary_directory)
            capture = temporary / "github-arguments.txt"
            fake_python = temporary / "python"
            fake_python.write_text(
                "#!/bin/bash\n"
                'if [[ " $* " == *" --check "* ]]; then exit 0; fi\n'
                'printf \'%s\\n\' "$@" > "$FAKE_GITHUB_ARGUMENTS"\n',
                encoding="utf-8",
            )
            fake_python.chmod(0o755)
            environment = os.environ.copy()
            environment.update(
                {
                    "FAKE_GITHUB_ARGUMENTS": str(capture),
                    "PYTHON_EXECUTABLE": str(fake_python),
                    "TROPHOSOME_NOTIFY_EMAIL": "off",
                }
            )
            result = subprocess.run(
                [
                    "/bin/bash",
                    "-c",
                    (
                        f"source {helper!s}; "
                        "trophosome_run_with_completion_notification "
                        f"'Synthetic failure' {REPOSITORY!s} /bin/bash -c 'exit 7'"
                    ),
                ],
                cwd=REPOSITORY,
                env=environment,
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(result.returncode, 7)
            arguments = capture.read_text(encoding="utf-8").splitlines()
            self.assertIn("--outcome", arguments)
            self.assertIn("FAILED", arguments)
            self.assertIn("--exit-status", arguments)
            self.assertIn("7", arguments)


if __name__ == "__main__":
    unittest.main()
