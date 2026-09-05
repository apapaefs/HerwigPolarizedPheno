#!/usr/bin/env python3
"""Generate schema-6 registry descriptors for the SIDIS tranches."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
COMMON_SCALES = [
    "/Herwig/MatrixElements/MEDISNC",
    "/Herwig/MatrixElements/MEDISNCPol",
    "/Herwig/MatrixElements/PowhegMEDISNC",
    "/Herwig/MatrixElements/PowhegMEDISNCPol",
]
SPECIFICATIONS: dict[str, dict[str, Any]] = {
    "COMPASS_2026_I3096394": {
        "title": "COMPASS corrected isoscalar pion, kaon and charged-hadron SIDIS multiplicities",
        "kind": "unpolarized_sidis", "postprocessor": "compass_sidis_multiplicity",
        "targets": {"P": "proton", "N": "neutron"},
        "outputs": {"D": {"P": .5, "N": .5}}, "helicities": {"00": [0,0]},
        "doi": "10.17182/hepdata.169859.v1", "source": "https://www.hepdata.net/record/ins3096394?version=1",
        "label": "h+/-, pi+/-, K+/- corrected multiplicities", "seed": 3096394,
    },
    "COMPASS_2025_I2840545": {
        "title": "COMPASS hydrogen pion, kaon and charged-hadron SIDIS multiplicities",
        "kind": "unpolarized_sidis", "postprocessor": "compass_sidis_multiplicity",
        "targets": {"P": "proton"}, "outputs": {"H": {"P": 1.}},
        "helicities": {"00": [0,0]}, "doi": "10.17182/hepdata.159544.v1",
        "source": "https://www.hepdata.net/record/ins2840545?version=1",
        "label": "h+/-, pi+/-, K+/- proton multiplicities", "seed": 2840545,
    },
    "COMPASS_2010_I862410": {
        "title": "COMPASS proton identified-hadron longitudinal SIDIS asymmetries",
        "kind": "polarized_sidis", "postprocessor": "compass_sidis_a1",
        "targets": {"P": "proton"}, "outputs": {"H": {"P": 1.}},
        "helicities": {"PP": [1,1], "PM": [1,-1], "MP": [-1,1], "MM": [-1,-1]},
        "source": "https://arxiv.org/abs/1007.4061", "label": "A1p pi+/-, K+/-",
        "seed": 862410, "longitudinal_target_scale": 1.,
    },
    "HERMES_2013_I1208547": {
        "title": "HERMES hydrogen/deuterium pion and kaon SIDIS multiplicities",
        "kind": "unpolarized_sidis", "postprocessor": "hermes_sidis_multiplicity",
        "targets": {"P": "proton", "N": "neutron"},
        "outputs": {"H": {"P": 1.}, "D": {"P": .5, "N": .5}},
        "helicities": {"00": [0,0]}, "doi": "10.17182/hepdata.62097.v1",
        "source": "https://www.hepdata.net/record/ins1208547?version=1",
        "label": "H/D pi+/-, K+/- five-binning multiplicities", "seed": 1208547,
        "projections": [
            "published_flattened_cells",
            "separately_integrated_numerator_denominator_projections",
        ],
    },
    "COMPASS_2018_I1624692": {
        "title": "COMPASS transverse-momentum-dependent charged-hadron SIDIS multiplicities",
        "kind": "unpolarized_sidis", "postprocessor": "compass_sidis_multiplicity",
        "targets": {"P": "proton", "N": "neutron"},
        "outputs": {"D": {"P": .5, "N": .5}}, "helicities": {"00": [0,0]},
        "doi": "10.17182/hepdata.83542.v1",
        "source": "https://www.hepdata.net/record/ins1624692?version=1",
        "label": "h+/- multiplicities in z and PhT2", "seed": 1624692,
    },
    "COMPASS_2020_I1788430": {
        "title": "COMPASS high-z antiproton/proton and negative/positive kaon SIDIS ratios",
        "kind": "unpolarized_sidis", "postprocessor": "compass_sidis_charge_ratio",
        "targets": {"P": "proton", "N": "neutron"},
        "outputs": {"D": {"P": .5, "N": .5}}, "helicities": {"00": [0,0]},
        "source": "https://arxiv.org/abs/2003.11791",
        "raw_source": "arXiv-2003.11791-source.tar.gz",
        "label": "high-z antiproton/proton and K-/K+ ratios", "seed": 1788430,
        "observable": "charge_ratio",
        "observable_level": "stable identified hadron",
    },
}


def descriptor(measurement: str, specification: dict[str, Any]) -> dict[str, Any]:
    support = [
        "analyses/rivet/dis/COMPASSInclusiveDIS.hh",
        "analyses/rivet/dis/COMPASSSIDIS.hh",
    ]
    if measurement != "COMPASS_2010_I862410":
        support.append("analyses/rivet/dis/SIDISTrancheBinning.hh")
    if measurement in {"COMPASS_2025_I2840545", "COMPASS_2026_I3096394"}:
        support.append("analyses/rivet/dis/COMPASSModernMultiplicity.hh")
        support.append("analyses/rivet/dis/COMPASSMultiplicityFiducial.hh")
    reference: dict[str, Any] = {
        "snapshot": f"data/phenomenology/{measurement}/reference.json",
        "source_manifest": f"data/phenomenology/{measurement}/source-manifest.json",
        "source_url": specification["source"],
    }
    if "doi" in specification:
        reference.update({
            "raw_snapshot": f"data/phenomenology/{measurement}/raw/record.json",
            "record_doi": specification["doi"],
        })
    else:
        reference["raw_snapshot"] = (
            f"data/phenomenology/{measurement}/raw/"
            f"{specification.get('raw_source', 'arXiv-1007.4061-source.tar.gz')}"
        )
    polarized = specification["kind"] == "polarized_sidis"
    config: dict[str, Any] = {"target_outputs": specification["outputs"]}
    if polarized:
        config["longitudinal_target_scale"] = specification["longitudinal_target_scale"]
    return {
        "schema_version": 6, "id": measurement, "title": specification["title"],
        "process_kind": specification["kind"],
        "analysis": {
            "name": measurement, "source": f"analyses/rivet/dis/{measurement}.cc",
            "info": f"analyses/rivet/dis/{measurement}.info",
            "plot": f"analyses/rivet/dis/{measurement}.plot", "support_files": support,
            "plugin": f"Rivet{measurement}.so",
            "reference_yoda": f"analyses/rivet/dis/{measurement}.yoda.gz",
        },
        "reference": reference,
        "cards": {
            "directory": f"cards/phenomenology/{measurement}",
            "common": f"{measurement}-Common.in",
            "stem_pattern": f"{measurement}_{{target}}_{{helicity}}-{{contribution}}",
            "polarized_pdf_object": "/Herwig/Partons/COMPASSSIDISDiffPDF",
            "unpolarized_pdf_object": "/Herwig/Partons/COMPASSSIDISPDF",
            "helicities": specification["helicities"],
            "contributions": {"POSNLO": "PositiveNLO", "NEGNLO": "NegativeNLO"},
            "target_components": specification["targets"], "scale_objects": COMMON_SCALES,
        },
        "channels": {"sidis": {
            "label": specification["label"], "targets": list(specification["outputs"]),
            "projections": specification.get(
                "projections", ["published_flattened_cells", "readable_slices"]
            ),
            "observable_level": specification.get(
                "observable_level", "stable charged hadron"
            ),
        }},
        "families": {"nominal": {
            "label": "POWHEG NLO+PS" + ("" if polarized else " (explicit 00 mode)"),
            "helicities": list(specification["helicities"]), "mpi": "off",
            "observable_level": "stable_particle",
        }},
        "campaign": {
            "default_jobs": 4, "default_shards": 1,
            "default_seed_base": int(specification["seed"]),
            "default_events": {"POSNLO": 300000, "NEGNLO": 30000},
            "smoke_events": 100,
        },
        "postprocessor": specification["postprocessor"],
        "postprocess_config": config,
        "scales": {"central": 1., "points": [.5,1.,2.],
                   "mapping": "DIS ScaleFactor = mu/mu0"},
        "pdf_ensembles": {
            "polarized": {"set": "NNPDFpol20_nlo_as_01180", "central_member": 0,
                          "replica_members": [1,100], "active": polarized},
            "unpolarized": {"set": "NNPDF40_nlo_pch_as_01180", "central_member": 0,
                            "replica_members": [1,100], "active": True},
            "combination": ("independent polarized/unpolarized replicas plus scale envelope"
                            if polarized else "active unpolarized replicas plus scale envelope"),
        },
        "goodness_of_fit": {
            "primary_observable": specification.get(
                "observable", "A1" if polarized else "multiplicity"
            ),
            "source_contract": "see checksum-pinned normalized reference",
        },
        "physics": {
            "beam": "fixed-target lepton DIS", "exchange": "photon",
            "target": specification["outputs"],
            "hard_process_accuracy": "POWHEG NLO+PS",
            "event_simulation": "QCD shower, hadronization and remnants; QED shower and MPI off",
        },
    }


def main() -> None:
    directory = ROOT / "config/phenomenology"
    for measurement, specification in SPECIFICATIONS.items():
        path = directory / f"{measurement}.json"
        path.write_text(json.dumps(descriptor(measurement, specification),
                                   indent=2, sort_keys=False) + "\n",
                        encoding="utf-8")


if __name__ == "__main__":
    main()
