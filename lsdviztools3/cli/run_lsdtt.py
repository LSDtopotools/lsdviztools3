"""CLI: lsdviz3-lsdtt — run an lsdtt3 tool with a parameter file.

Also the place to point lsdviztools3 at your lsdtt3 executables:
    lsdviz3-lsdtt --set-lsdtt-dir /path/to/lsdtt3/bin
"""

from __future__ import annotations

import argparse

from ..driver import run_lsdtt, write_user_config


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for ``lsdviz3-lsdtt``."""
    parser = argparse.ArgumentParser(
        prog="lsdviz3-lsdtt",
        description="Run an lsdtt3 command-line tool, or configure its location.",
    )
    parser.add_argument(
        "tool",
        nargs="?",
        help="lsdtt3 tool name, e.g. lsdtt-channel-extraction.",
    )
    parser.add_argument(
        "--params",
        help="Path to an existing lsdtt parameter (.driver) file.",
    )
    parser.add_argument(
        "--data-path",
        help="The <path> argument lsdtt3 expects (defaults to the params file dir).",
    )
    parser.add_argument(
        "--lsdtt-dir",
        help="Directory containing the lsdtt-* executables (overrides env/config).",
    )
    parser.add_argument(
        "--set-lsdtt-dir",
        metavar="DIR",
        help="Persist the lsdtt3 bin directory to the user config file and exit.",
    )
    parser.add_argument("--cwd", help="Working directory to run the tool in.")
    return parser


def main(argv=None) -> int:
    """Entry point: configure the lsdtt3 location or run a tool via the driver."""
    args = build_parser().parse_args(argv)

    if args.set_lsdtt_dir:
        path = write_user_config(args.set_lsdtt_dir)
        print(f"Wrote lsdtt3 directory to {path}")
        return 0

    if not args.tool or not args.params:
        raise SystemExit("error: a tool name and --params are required to run lsdtt3")

    result = run_lsdtt(
        args.tool,
        args.params,
        data_path=args.data_path,
        lsdtt_dir=args.lsdtt_dir,
        cwd=args.cwd,
        check=False,
    )
    if result.stdout:
        print(result.stdout, end="")
    if result.stderr:
        print(result.stderr, end="")
    return result.returncode


if __name__ == "__main__":
    raise SystemExit(main())
