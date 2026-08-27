#!/usr/bin/env python3
"""Generate the finite COMPASS SIDIS target/helicity/NLO base-card matrix."""

from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
POLARIZATIONS = {
    "PP": (1, 1), "PM": (1, -1), "MP": (-1, 1), "MM": (-1, -1),
    "00": (0, 0),
}
CONTRIBUTIONS = {"POSNLO": "PositiveNLO", "NEGNLO": "NegativeNLO"}


def _card(measurement: str, target: str, helicity: str, contribution: str) -> str:
    first, second = POLARIZATIONS[helicity]
    lines = [
        "# -*- ThePEG-repository -*-",
        f"read {measurement}-Common.in",
    ]
    if target == "N":
        lines.extend([
            "set /Herwig/EventHandlers/EventHandler:BeamB /Herwig/Particles/n0",
            "set /Herwig/EventHandlers/FixedTargetLuminosity:TargetParticle /Herwig/Particles/n0",
        ])
    lines.extend([
        f"set /Herwig/Partons/EPPolarizedExtractor:FirstLongitudinalPolarization {first}",
        f"set /Herwig/Partons/EPPolarizedExtractor:SecondLongitudinalPolarization {second}",
        f"set /Herwig/MatrixElements/PowhegMEDISNCPol:Contribution {CONTRIBUTIONS[contribution]}",
        "cd /Herwig/Generators",
        f"saverun {measurement}_{target}_{helicity}-{contribution} EventGenerator",
        "",
    ])
    return "\n".join(lines)


def main() -> None:
    matrices = {
        "COMPASS_2009_I820721": ("PP", "PM", "MP", "MM"),
        "COMPASS_2017_I1444985": ("00",),
        "COMPASS_2017_I1483098": ("00",),
    }
    for measurement, helicities in matrices.items():
        directory = ROOT / "cards" / "phenomenology" / measurement
        directory.mkdir(parents=True, exist_ok=True)
        for target in ("P", "N"):
            for helicity in helicities:
                for contribution in CONTRIBUTIONS:
                    destination = directory / f"{measurement}_{target}_{helicity}-{contribution}.in"
                    destination.write_text(
                        _card(measurement, target, helicity, contribution),
                        encoding="utf-8",
                    )


if __name__ == "__main__":
    main()
