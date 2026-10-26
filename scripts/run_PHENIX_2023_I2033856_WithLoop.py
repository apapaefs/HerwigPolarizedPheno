#!/usr/bin/env python3
"""Three-component PHENIX campaign, using the existing cross-section sum.

--hard-pt-events can set photon, fragmentation, and gammajet_loop budgets.
--photon-enhancement applies only to fragmentation, with detuning fixed at 2.
"""
from pathlib import Path
from run_PHENIX_2023_I2033856 import run_measurement

MEASUREMENT = "PHENIX_2023_I2033856_WithLoop"


def main(argv=None):
    return run_measurement(MEASUREMENT, argv, Path(__file__).name)


if __name__ == "__main__":
    raise SystemExit(main())
