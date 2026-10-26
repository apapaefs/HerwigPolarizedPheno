#!/usr/bin/env python3
"""PHENIX2023 campaign wrapper; --output-root can place results outside this repository."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

import run_phenomenology_campaign as campaign

MEASUREMENT = "PHENIX_2023_I2033856"


def run_measurement(measurement: str, argv: list[str] | None = None,
                    program: str | None = None) -> int:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--output-root", type=Path)
    options, remaining = parser.parse_known_args(sys.argv[1:] if argv is None else argv)
    if not remaining or remaining == ["--help"]:
        print(f"Usage: {program or Path(__file__).name} [--output-root DIRECTORY] COMMAND [campaign options]")
        print(f"Measurement is fixed to {measurement}.")
        print(f"Results are stored at DIRECTORY/{measurement}/TAG.")
        print("Optional --hard-pt-edges 4,8,12,16,20,inf generates and adds Born-pT slices (GeV).")
        print("--lo-events and --shards apply per helicity, component, and slice.")
        print("--hard-pt-events CHANNEL=N1,N2,... overrides events per helicity in each slice.")
        print("Repeat for photon/fragmentation/gammajet_loop; other components retain --lo-events.")
        if measurement in {MEASUREMENT, MEASUREMENT + "_WithLoop"}:
            print("--photon-enhancement C enhances QCD-dijet photon showering with fixed detuning 2.")
            print("Use C=1 for the control; WithLoop defaults to C=12, other variants preserve base settings.")
        if measurement == MEASUREMENT + "_WithLoop":
            print("Default: five slices including the open tail; 125M events total at central scale.")
            print("--smoke uses 100 events per helicity/component/slice (6000 total).")
            print("Explicit --lo-events or --hard-pt-edges replaces the default event allocation.")
        print("See docs/PHENIX_2023_I2033856.md for budgets, prerequisites, and examples.")
        campaign.make_parser().print_help()
        return 0
    if any(arg == "--measurement" or arg.startswith("--measurement=") for arg in remaining):
        parser.error(f"This wrapper fixes --measurement to {measurement}")
    if options.output_root is not None:
        campaign.CAMPAIGN_ROOT = options.output_root.expanduser().resolve()
    command, *arguments = remaining
    if command == "list":
        return campaign.main([command, *arguments])
    return campaign.main([command, "--measurement", measurement, *arguments])


def main(argv: list[str] | None = None) -> int:
    return run_measurement(MEASUREMENT, argv)


if __name__ == "__main__":
    raise SystemExit(main())
