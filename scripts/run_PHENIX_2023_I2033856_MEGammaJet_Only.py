#!/usr/bin/env python3
"""Reproduce the original photon+jet definition in a separate campaign."""
from pathlib import Path
from run_PHENIX_2023_I2033856 import run_measurement

MEASUREMENT = "PHENIX_2023_I2033856_MEGammaJet-Only"


def main(argv=None):
    return run_measurement(MEASUREMENT, argv, Path(__file__).name)


if __name__ == "__main__":
    raise SystemExit(main())
