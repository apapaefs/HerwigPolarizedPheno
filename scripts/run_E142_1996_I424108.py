#!/usr/bin/env python3
"""E142 campaign wrapper; --output-root can place results outside this repository."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

import run_experimental_campaign as campaign

MEASUREMENT = "E142_1996_I424108"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--output-root", type=Path)
    options, remaining = parser.parse_known_args(sys.argv[1:] if argv is None else argv)
    if not remaining or remaining == ["--help"]:
        print(f"Usage: {Path(__file__).name} [--output-root DIRECTORY] COMMAND [campaign options]")
        print(f"Measurement is fixed to {MEASUREMENT}.")
        print("Results are stored at DIRECTORY/E142_1996_I424108/TAG.")
        campaign.make_parser().print_help()
        return 0
    if any(arg == "--measurement" or arg.startswith("--measurement=") for arg in remaining):
        parser.error("This wrapper fixes --measurement to E142_1996_I424108")
    if options.output_root is not None:
        campaign.CAMPAIGN_ROOT = options.output_root.expanduser().resolve()
    command, *arguments = remaining
    if command == "list":
        return campaign.main([command, *arguments])
    return campaign.main([command, "--measurement", MEASUREMENT, *arguments])


if __name__ == "__main__":
    raise SystemExit(main())
