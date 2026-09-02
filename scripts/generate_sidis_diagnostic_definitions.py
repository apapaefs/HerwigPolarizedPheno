#!/usr/bin/env python3
"""Generate the exact-bin internal definition for HERMES azimuthal moments."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
MEASUREMENT = "HERMES_2013_I1111237"
SPECIES = ("hplus", "hminus", "piplus", "piminus", "kplus", "kminus")
TARGETS = {"H": {"P": 1.0}, "D": {"P": .5, "N": .5}}


def definition() -> dict[str, Any]:
    edges = [float(index) for index in range(901)]
    datasets: dict[str, Any] = {}
    observables: dict[str, Any] = {}
    for target, weights in TARGETS.items():
        for species in SPECIES:
            for harmonic in (1, 2):
                identifier = f"{target}_{species}_cos{harmonic}_xyzt"
                raw_identifier = f"{species}_cos{harmonic}"
                numerator = f"MomentNumerator_{raw_identifier}"
                raw_objects = {
                    "numerator": numerator,
                    "denominator": f"MomentDenominator_{raw_identifier}",
                    "covariance_positive": f"CovariancePositive_{raw_identifier}",
                    "covariance_negative": f"CovarianceNegative_{raw_identifier}",
                }
                datasets[identifier] = {
                    "id": identifier,
                    "data_available": False,
                    "published_target": target,
                    "target_weights": weights,
                    "species": species,
                    "harmonic": harmonic,
                    "edges": edges,
                    "rivet_path": f"/{MEASUREMENT}/DIAGNOSTICS/{identifier}",
                    "raw_objects": raw_objects,
                }
                observables[identifier] = {
                    "raw_object": numerator,
                    "edges": edges,
                    "helicity_resolved": False,
                    "description": (
                        f"{target} {species} direct <cos({harmonic} phi_h)> "
                        "in the exact flattened 5x5x6x6 HERMES cell order"
                    ),
                }
    return {
        "schema_version": 1,
        "kind": "internal_observable_definition",
        "measurement": MEASUREMENT,
        "observable": "<cos(phi_h)>_UU and <cos(2phi_h)>_UU",
        "data_status": {
            "available": False,
            "four_dimensional_cells_per_sample": 900,
            "sample_groups": 12,
            "expected_moment_values": 21600,
            "expected_covariance": True,
            "official_endpoint": "https://www-hermes.desy.de/cosnphi/",
            "policy": "no numerical values or covariance entries are fabricated",
        },
        "beam": {"pid": "e+", "energy_gev": 27.6},
        "selection": {
            "q2_min_gev2": 1.0,
            "w2_min_gev2": 10.0,
            "x": [.023, .6],
            "y": [.2, .85],
            "z": [.2, 1.0],
            "xF_min": .2,
            "unidentified_and_pion_momentum_gev": [1.0, 15.0],
            "kaon_momentum_gev": [2.0, 15.0],
        },
        "binning": {
            "x": [.023, .042, .078, .145, .27, .6],
            "y": [.2, .3, .45, .6, .7, .85],
            "z": [.2, .3, .4, .5, .6, .75, 1.0],
            "pt_gev": [.05, .2, .35, .5, .7, 1.0, 1.3],
            "cells_per_target_species": 900,
            "flattening": "((x_bin * 5) + y_bin) * 36 + z_bin * 6 + pt_bin",
        },
        "estimator": {
            "numerator": "sum cos(n phi_h)",
            "denominator": "identified-hadron count",
            "same_event_covariance": True,
        },
        "target_outputs": TARGETS,
        "datasets": datasets,
        "observables": observables,
        "provenance": {
            "source_manifest": (
                "data/phenomenology/HERMES_2013_I1111237/source-manifest.json"
            ),
            "paper_source_sha256": (
                "a09eca719c7fbde602270425c25c86107697d1c1f5ed472f575f8ed03ae5f56a"
            ),
            "paper": "arXiv:1204.4161",
        },
        "interpretation": (
            "stable-hadron generator diagnostic without detector unfolding, "
            "radiative corrections, diffractive subtraction, or a claim of "
            "controlled TMD/twist-3 accuracy"
        ),
    }


def main() -> None:
    destination = (
        ROOT / "data/phenomenology/HERMES_2013_I1111237/analysis-definition.json"
    )
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(definition(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
