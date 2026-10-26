#!/usr/bin/env python3
"""Run only the polarized gg -> gamma g loop component with expanded PHENIX cuts.

The channel key 'photon' here refers to MEPP2GammaJetLoop, not MEGammaJet.
Use --hard-pt-events photon=... for per-slice budgets. The usual runner
normalizes each sample to its own cross section; no external loop fraction.
"""
from pathlib import Path
from run_PHENIX_2023_I2033856 import run_measurement

MEASUREMENT = "PHENIX_2023_I2033856_MEGammaJetLoop-Only"


def main(argv=None):
    return run_measurement(MEASUREMENT, argv, Path(__file__).name)


if __name__ == "__main__":
    raise SystemExit(main())
