"""Command-line interface for validation and reproducible runs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from trophosome import __version__
from trophosome.config import load_config
from trophosome.simulation import run_simulation


class _HelpFormatter(
    argparse.RawDescriptionHelpFormatter,
    argparse.ArgumentDefaultsHelpFormatter,
):
    """Preserve examples while also showing useful default values."""


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="trophosome",
        formatter_class=_HelpFormatter,
        description="""\
Simulate the evolution of environmentally acquired microbial symbionts.

Trophosome follows labelled symbiont strains through a repeated biological cycle:

  environmental reservoir -> host infection -> within-host reproduction
  -> symbiont release -> environmental reservoir

The software was inspired by tubeworm symbioses, but the same cycle can be
parameterised for other host--microbial associations with environmental
acquisition.

Typical workflow:
  1. Copy a TOML configuration from the configs/ folder.
  2. Edit its biological parameters for the experiment.
  3. Check it with `trophosome validate`.
  4. Run it with `trophosome run` and give the results a new output folder.
""",
        epilog="""\
Examples:
  trophosome validate configs/smoke.toml
  trophosome run configs/smoke.toml --output results/my_first_run

Use `trophosome COMMAND -h` for help with a particular command.
The README explains installation, parameters and result files in plain language.
""",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
        help="show the installed trophosome version and exit",
    )
    commands = parser.add_subparsers(
        dest="command",
        required=True,
        title="commands",
        metavar="COMMAND",
    )

    validate = commands.add_parser(
        "validate",
        help="check an experiment configuration without running the model",
        formatter_class=_HelpFormatter,
        description="""\
Check an experiment configuration before committing time and computing resources.

This command reads the TOML file, checks that parameter values are valid and
mutually consistent, and prints the complete interpreted configuration. It also
reports feasibility warnings, such as an unusually large expected mutation
workload. It does not simulate hosts and does not create a results folder.
""",
        epilog="""\
Example:
  trophosome validate configs/my_experiment.toml

If no error is reported, use the same file with `trophosome run`.
""",
    )
    validate.add_argument(
        "config",
        type=Path,
        metavar="CONFIG.toml",
        help="TOML file containing the biological and computational parameters",
    )

    run = commands.add_parser(
        "run",
        help="run an experiment and write its results",
        formatter_class=_HelpFormatter,
        description="""\
Run the maintained strain-count model described by a plain-text TOML
configuration file.

The model independently infects each host from the environmental reservoir,
simulates within-host reproduction, mutation and optional selection, pools the
released symbionts, and updates the environmental strain composition. It repeats
this cycle for the requested host generations and stochastic replicates.

The output directory contains CSV tables that can be opened in R, Python or a
spreadsheet program, together with the resolved configuration, provenance and
recovery information. Use a different output directory for every experiment.
If an interrupted run has a valid recovery checkpoint, repeat the original
command with --resume and the same output directory.
""",
        epilog="""\
Example:
  trophosome run configs/my_experiment.toml --output results/my_experiment

Resume an interrupted run:
  trophosome run configs/my_experiment.toml --output results/my_experiment --resume

Recommended first run:
  trophosome validate configs/smoke.toml
  trophosome run configs/smoke.toml --output results/smoke

Configuration sections:
  [environment]  starting environmental strains and effective reservoir size
  [host]         host number, infection, within-host growth and release
  [evolution]    mutation, fitness effects and selection switches
  [output]       how frequently and for which hosts detailed results are saved
  [execution]    number of CPU workers and host-processing batch size

Start interpretation with:
  host_generation_summary.csv  one row per replicate and host generation
  environment_counts.csv       abundance of every environmental strain over time
  resolved_config.json         exact parameter values used for the run

See the README for a guide to all output tables and configuration parameters.
""",
    )
    run.add_argument(
        "config",
        type=Path,
        metavar="CONFIG.toml",
        help="validated TOML file describing the experiment",
    )
    run.add_argument(
        "--output",
        type=Path,
        required=True,
        metavar="DIRECTORY",
        help="new directory in which to save all results",
    )
    run.add_argument(
        "--repository",
        type=Path,
        default=Path.cwd(),
        metavar="PROJECT_DIRECTORY",
        help=(
            "software project folder used only to record which code version "
            "produced the results; normally omit this option"
        ),
    )
    run.add_argument(
        "--resume",
        action="store_true",
        help=(
            "continue from the newest valid recovery checkpoint, safely "
            "discarding any incomplete rows written after it"
        ),
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    config = load_config(args.config)
    if args.command == "validate":
        print(json.dumps(config.to_dict(), indent=2, sort_keys=True))
        for warning in config.feasibility_warnings():
            print(f"WARNING: {warning}")
        return 0
    summaries = run_simulation(
        config,
        args.output,
        args.repository,
        resume=args.resume,
    )
    print(
        json.dumps(
            {
                "output": str(args.output.resolve()),
                "host_generations": len(summaries),
                "warnings": config.feasibility_warnings(),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
