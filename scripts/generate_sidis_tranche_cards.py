#!/usr/bin/env python3
"""Generate the finite target/helicity/signed-NLO cards for SIDIS analyses."""

from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE_MULTIPLICITY = ROOT / "cards/phenomenology/COMPASS_2017_I1444985/COMPASS_2017_I1444985-Common.in"
BASE_ASYMMETRY = ROOT / "cards/phenomenology/COMPASS_2009_I820721/COMPASS_2009_I820721-Common.in"
POLARIZATIONS = {
    "PP": (1, 1), "PM": (1, -1), "MP": (-1, 1), "MM": (-1, -1),
    "00": (0, 0),
}
CONTRIBUTIONS = {"POSNLO": "PositiveNLO", "NEGNLO": "NegativeNLO"}
SPECIFICATIONS = {
    "COMPASS_2026_I3096394": {"base": BASE_MULTIPLICITY, "targets": ("P","N"),
                               "helicities": ("00",), "maxy": ".70", "minw2": "25.0"},
    "COMPASS_2025_I2840545": {"base": BASE_MULTIPLICITY, "targets": ("P",),
                               "helicities": ("00",), "maxy": ".70", "minw2": "25.0"},
    "COMPASS_2010_I862410": {"base": BASE_ASYMMETRY, "targets": ("P",),
                              "helicities": ("PP","PM","MP","MM"),
                              "maxy": ".90", "minw2": "0.0"},
    "HERMES_2013_I1208547": {"base": BASE_MULTIPLICITY, "targets": ("P","N"),
                              "helicities": ("00",), "maxy": ".85", "minw2": "10.0",
                              "hermes": True},
    "COMPASS_2018_I1624692": {"base": BASE_MULTIPLICITY, "targets": ("P","N"),
                               "helicities": ("00",), "maxy": ".90", "minw2": "25.0"},
    "COMPASS_2020_I1788430": {"base": BASE_MULTIPLICITY, "targets": ("P","N"),
                               "helicities": ("00",), "maxy": "1.00", "minw2": "25.0"},
}


def common(measurement: str, specification: dict[str, object]) -> str:
    base_path = Path(specification["base"])
    base_id = "COMPASS_2009_I820721" if base_path == BASE_ASYMMETRY else "COMPASS_2017_I1444985"
    text = base_path.read_text(encoding="utf-8").replace(base_id, measurement)
    text = text.replace("set /Herwig/Cuts/NeutralCurrentCut:Maxy 0.70",
                        f"set /Herwig/Cuts/NeutralCurrentCut:Maxy {specification['maxy']}")
    text = text.replace("set /Herwig/Cuts/NeutralCurrentCut:Maxy 0.90",
                        f"set /Herwig/Cuts/NeutralCurrentCut:Maxy {specification['maxy']}")
    text = text.replace("set /Herwig/Cuts/NeutralCurrentCut:MinW2 25.0*GeV2",
                        f"set /Herwig/Cuts/NeutralCurrentCut:MinW2 {specification['minw2']}*GeV2")
    text = text.replace("set /Herwig/Cuts/NeutralCurrentCut:MinW2 0.0*GeV2",
                        f"set /Herwig/Cuts/NeutralCurrentCut:MinW2 {specification['minw2']}*GeV2")
    if specification.get("hermes"):
        text = text.replace("COMPASS 2017 pion/unidentified-hadron multiplicities at 160 GeV",
                            "HERMES identified-hadron multiplicities at 27.6 GeV")
        text = text.replace("/Herwig/Particles/mu-", "/Herwig/Particles/e-")
        text = text.replace("/Herwig/Particles/mu+", "/Herwig/Particles/e+")
        text = text.replace("160.0*GeV", "27.6*GeV")
    return text


def card(measurement: str, target: str, helicity: str, contribution: str) -> str:
    first, second = POLARIZATIONS[helicity]
    lines = ["# -*- ThePEG-repository -*-", f"read {measurement}-Common.in"]
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
    for measurement, specification in SPECIFICATIONS.items():
        directory = ROOT / "cards/phenomenology" / measurement
        directory.mkdir(parents=True, exist_ok=True)
        (directory / f"{measurement}-Common.in").write_text(
            common(measurement, specification), encoding="utf-8"
        )
        for target in specification["targets"]:
            for helicity in specification["helicities"]:
                for contribution in CONTRIBUTIONS:
                    (directory / f"{measurement}_{target}_{helicity}-{contribution}.in").write_text(
                        card(measurement, target, helicity, contribution), encoding="utf-8"
                    )


if __name__ == "__main__":
    main()
