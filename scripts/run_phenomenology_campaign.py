#!/usr/bin/env python3
"""Unified polarized-data campaign runner for fixed-target DIS and RHIC.

The fixed-target descriptors are executed through the established experimental
engine.  Registry-native proton--proton, polarized-jet, and SIDIS descriptors
share immutable manifests, sharding, tracker, safe-resume, and fresh-seed
recovery while retaining explicit process, target, observable-level, family,
helicity, perturbative-contribution, PDF, hard-scale, and MPI axes.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import copy
import csv
import hashlib
import html
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import run_experimental_campaign as experimental
import compass_sidis_postprocess as compass_sidis
import polarized_sidis_postprocess as sidis
import runtime_provenance as provenance
import star_comparison_policy as star_policy
from phenomenology_reference_data import (
    REFERENCE_YODA_PATHS,
    ReferenceDataError,
    fetch_and_validate,
    validate_vendored,
    write_reference_yoda,
)


SCRIPT_PATH = Path(__file__).resolve()
DISPOL_ROOT = SCRIPT_PATH.parents[1]
REGISTRY_DIR = DISPOL_ROOT / "config" / "phenomenology"
CAMPAIGN_ROOT = DISPOL_ROOT / "campaigns" / "phenomenology"
TAG_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


MC_POLJETSHAPES_RATIO_PREFIXES = (
    "SigmaPP_", "SigmaPM_", "SigmaMP_", "SigmaMM_", "SigmaUU_",
    "ShapePP_", "ShapePM_", "ShapeMP_", "ShapeMM_", "ShapeUU_",
    "R32_", "R43_", "ThirdJetVeto_",
)


# This curated hierarchy is independent of the pilot sensitivity
# ranking.  It prevents a noisy pilot or a look-elsewhere fluctuation from
# silently deciding which production plots are promoted to the focused page.
MC_POLJETSHAPES_FOCUS_SECTIONS = (
    {
        "title": "1. Resolved-radiation plane moments",
        "why": (
            "These directly ask whether a resolved third or fourth jet retains "
            "a coherent cos(2 psi) modulation. The x dependence separates "
            "soft extra radiation from genuinely hard additional jets."
        ),
        "plots": (
            (
                "C2UU_resolved_dphi31_vs_pt31",
                "C2UU of the jet-3 plane versus pT3/pT1",
                "difference",
            ),
            (
                "C2LL_resolved_dphi31_vs_pt31",
                "C2LL of the jet-3 plane versus pT3/pT1",
                "difference",
            ),
            (
                "C2UU_resolved_dpsi34_vs_pt41",
                "C2UU of the jet-3/jet-4 planes versus pT4/pT1",
                "difference",
            ),
            (
                "C2LL_resolved_dpsi34_vs_pt41",
                "C2LL of the jet-3/jet-4 planes versus pT4/pT1",
                "difference",
            ),
            (
                "S2UU_resolved_dphi31_vs_pt31",
                "Sine null for the jet-3 plane",
                "difference",
            ),
            (
                "S2UU_resolved_dpsi34_vs_pt41",
                "Sine null for the jet-3/jet-4 planes",
                "difference",
            ),
        ),
    },
    {
        "title": "2. Hard energy-sharing declusterings",
        "why": (
            "The 0.25<z1<0.40 and z2>0.35 requirement suppresses strongly "
            "ordered soft branchings. The 1 and 2 GeV kT selections show "
            "whether any modulation strengthens in the perturbative core."
        ),
        "plots": (
            (
                "A2UU_dpsi12_j1_hardshare_kt1",
                "A2UU: jet 1, hard sharing, kT > 1 GeV",
                "difference",
            ),
            (
                "A2LL_dpsi12_j1_hardshare_kt1",
                "A2LL: jet 1, hard sharing, kT > 1 GeV",
                "difference",
            ),
            (
                "A2UU_dpsi12_j2_hardshare_kt1",
                "A2UU: jet 2, hard sharing, kT > 1 GeV",
                "difference",
            ),
            (
                "A2LL_dpsi12_j2_hardshare_kt1",
                "A2LL: jet 2, hard sharing, kT > 1 GeV",
                "difference",
            ),
            (
                "A2UU_dpsi12_j1_hardshare_kt2",
                "A2UU: jet 1, hard sharing, kT > 2 GeV",
                "difference",
            ),
            (
                "A2UU_dpsi12_j2_hardshare_kt2",
                "A2UU: jet 2, hard sharing, kT > 2 GeV",
                "difference",
            ),
        ),
    },
    {
        "title": "3. Hard-sharing squeezed EEEC",
        "why": (
            "The same hard-splitting gate is applied before constructing the "
            "all-particle squeezed three-point correlator. This checks the "
            "effect with an energy-weighted analyser rather than one selected "
            "declustering pair."
        ),
        "plots": (
            (
                "ShapeUU_eeec_squeezed_j1_hardshare_kt1",
                "Squeezed EEEC: jet 1, hard sharing, kT > 1 GeV",
                "ratio",
            ),
            (
                "ShapeUU_eeec_squeezed_j2_hardshare_kt1",
                "Squeezed EEEC: jet 2, hard sharing, kT > 1 GeV",
                "ratio",
            ),
            (
                "ShapeUU_eeec_squeezed_j1_hardshare_kt2",
                "Squeezed EEEC: jet 1, hard sharing, kT > 2 GeV",
                "ratio",
            ),
            (
                "ShapeUU_eeec_squeezed_j2_hardshare_kt2",
                "Squeezed EEEC: jet 2, hard sharing, kT > 2 GeV",
                "ratio",
            ),
        ),
    },
    {
        "title": "4. Particle-level gluon-enriched proxies",
        "why": (
            "High constituent multiplicity and low pTD are experimentally "
            "constructible gluon-enriched categories. Their complementary "
            "bins are shown beside them; no truth-flavour label is used."
        ),
        "plots": (
            (
                "A2UU_dpsi12_j1_hardshare_kt1_nconst_high",
                "A2UU: jet 1, Nconst >= 8",
                "difference",
            ),
            (
                "A2UU_dpsi12_j1_hardshare_kt1_nconst_low",
                "A2UU: jet 1, Nconst < 8",
                "difference",
            ),
            (
                "A2UU_dpsi12_j2_hardshare_kt1_nconst_high",
                "A2UU: jet 2, Nconst >= 8",
                "difference",
            ),
            (
                "A2UU_dpsi12_j2_hardshare_kt1_nconst_low",
                "A2UU: jet 2, Nconst < 8",
                "difference",
            ),
            (
                "A2UU_dpsi12_j1_hardshare_kt1_ptd_low",
                "A2UU: jet 1, pTD < 0.45",
                "difference",
            ),
            (
                "A2UU_dpsi12_j1_hardshare_kt1_ptd_high",
                "A2UU: jet 1, pTD >= 0.45",
                "difference",
            ),
        ),
    },
    {
        "title": "5. Inclusive splitting-plane moments",
        "why": (
            "Primary tests of a coherent cos(2 psi) shower-spin modulation. "
            "The companion panel is spin on minus spin off; moment ratios are "
            "intentionally avoided because moments may cross zero."
        ),
        "plots": (
            ("A2UU_dpsi12_j1_loose", "A2UU: jet 1, loose", "difference"),
            ("A2LL_dpsi12_j1_loose", "A2LL: jet 1, loose", "difference"),
            ("A2UU_dpsi12_j2_loose", "A2UU: jet 2, loose", "difference"),
            ("A2LL_dpsi12_j2_loose", "A2LL: jet 2, loose", "difference"),
            (
                "A2UU_dpsi12_j1_symmetric_secondary",
                "A2UU: jet 1, symmetric secondary",
                "difference",
            ),
            (
                "A2LL_dpsi12_j1_symmetric_secondary",
                "A2LL: jet 1, symmetric secondary",
                "difference",
            ),
            (
                "A2UU_dpsi12_j2_symmetric_secondary",
                "A2UU: jet 2, symmetric secondary",
                "difference",
            ),
            (
                "A2LL_dpsi12_j2_symmetric_secondary",
                "A2LL: jet 2, symmetric secondary",
                "difference",
            ),
        ),
    },
    {
        "title": "6. Inclusive splitting-plane shapes",
        "why": (
            "The normalized UU shapes expose the expected smooth even angular "
            "pattern while removing the inclusive normalization.  Delta-sigma "
            "LL is shown as a difference rather than a ratio because it is signed."
        ),
        "plots": (
            ("ShapeUU_dpsi12_j1_loose", "UU shape: jet 1, loose", "ratio"),
            ("ShapeUU_dpsi12_j2_loose", "UU shape: jet 2, loose", "ratio"),
            (
                "ShapeUU_dpsi12_j1_symmetric_secondary",
                "UU shape: jet 1, symmetric secondary",
                "ratio",
            ),
            (
                "ShapeUU_dpsi12_j2_symmetric_secondary",
                "UU shape: jet 2, symmetric secondary",
                "ratio",
            ),
            (
                "DeltaSigmaLL_dpsi12_j1_loose",
                "Delta-sigma LL: jet 1, loose",
                "difference",
            ),
            (
                "DeltaSigmaLL_dpsi12_j2_loose",
                "Delta-sigma LL: jet 2, loose",
                "difference",
            ),
            (
                "DeltaSigmaLL_dpsi12_j1_symmetric_secondary",
                "Delta-sigma LL: jet 1, symmetric secondary",
                "difference",
            ),
            (
                "DeltaSigmaLL_dpsi12_j2_symmetric_secondary",
                "Delta-sigma LL: jet 2, symmetric secondary",
                "difference",
            ),
        ),
    },
    {
        "title": "7. Independent angular confirmation",
        "why": (
            "These use different particle-level angular analysers.  A credible "
            "effect should be coherent across more than one construction and "
            "should not be driven by one sparse bin."
        ),
        "plots": (
            (
                "ShapeUU_interjet_dpsi11_kt05",
                "Inter-jet primary-plane angle, kT > 0.5 GeV",
                "ratio",
            ),
            (
                "ShapeUU_hardplane_primary_j1_kt05",
                "Hard-plane angle: jet 1, kT > 0.5 GeV",
                "ratio",
            ),
            (
                "ShapeUU_hardplane_primary_j2_kt05",
                "Hard-plane angle: jet 2, kT > 0.5 GeV",
                "ratio",
            ),
            ("ShapeUU_eeec_squeezed_j1", "Squeezed EEEC: jet 1", "ratio"),
            ("ShapeUU_eeec_squeezed_j2", "Squeezed EEEC: jet 2", "ratio"),
            ("ShapeUU_bz_angle", "Four-jet Bengtsson-Zerwas angle", "ratio"),
        ),
    },
    {
        "title": "8. Sine and symmetry null tests",
        "why": (
            "The B2 sine moments should be compatible with zero.  The S2 "
            "quadrupole spectra are supporting symmetry controls; unexpected "
            "structure here should be understood before a cosine signal is claimed."
        ),
        "plots": (
            ("B2UU_dpsi12_j1_loose", "B2UU: jet 1, loose", "difference"),
            ("B2LL_dpsi12_j1_loose", "B2LL: jet 1, loose", "difference"),
            ("B2UU_dpsi12_j2_loose", "B2UU: jet 2, loose", "difference"),
            ("B2LL_dpsi12_j2_loose", "B2LL: jet 2, loose", "difference"),
            ("SigmaUU_s2_beta1_j1", "UU S2,beta=1: jet 1", "ratio"),
            ("SigmaUU_s2_beta2_j1", "UU S2,beta=2: jet 1", "ratio"),
            ("SigmaUU_s2_beta1_j2", "UU S2,beta=1: jet 2", "ratio"),
            ("SigmaUU_s2_beta2_j2", "UU S2,beta=2: jet 2", "ratio"),
        ),
    },
    {
        "title": "9. Resolved third- and fourth-jet spectra",
        "why": (
            "These are the direct particle-level analogues of the resolved-"
            "radiation hardness tests that were useful in polarized DIS.  "
            "The UU spectrum tests the amount of additional radiation, while "
            "Delta-sigma LL and A_LL distinguish that from a longitudinal-"
            "spin asymmetry."
        ),
        "plots": (
            ("SigmaUU_jet3_pt", "Third-jet pT: UU cross section", "ratio"),
            (
                "DeltaSigmaLL_jet3_pt",
                "Third-jet pT: Delta-sigma LL",
                "difference",
            ),
            ("ALL_jet3_pt", "Third-jet pT: A_LL", "difference"),
            ("SigmaUU_jet4_pt", "Fourth-jet pT: UU cross section", "ratio"),
            (
                "DeltaSigmaLL_jet4_pt",
                "Fourth-jet pT: Delta-sigma LL",
                "difference",
            ),
            ("ALL_jet4_pt", "Fourth-jet pT: A_LL", "difference"),
        ),
    },
    {
        "title": "10. Radiation-rate scans",
        "why": (
            "These test indirect changes in resolved radiation.  They are "
            "secondary to the angular observables because shower spin "
            "correlations are expected to affect azimuthal structure more "
            "directly than inclusive emission probabilities."
        ),
        "plots": (
            ("R32_UU", "Inclusive R32", "ratio"),
            ("ThirdJetVeto_UU", "Third-jet veto efficiency", "ratio"),
            (
                "SigmaUU_pt31_cumulative_tail",
                "Cumulative pT3 / pT1 tail",
                "ratio",
            ),
            ("R43_UU", "Inclusive R43", "ratio"),
            (
                "SigmaUU_pt41_cumulative_tail",
                "Cumulative pT4 / pT1 tail",
                "ratio",
            ),
        ),
    },
    {
        "title": "11. Inclusive controls",
        "why": (
            "Large changes in these broad spectra would be surprising and "
            "should first trigger a configuration, normalization, or "
            "statistics check rather than a shower-spin interpretation."
        ),
        "plots": (
            ("SigmaUU_jet1_pt", "Leading-jet pT", "ratio"),
            ("SigmaUU_jet2_pt", "Subleading-jet pT", "ratio"),
            ("SigmaUU_jet1_eta", "Leading-jet eta", "ratio"),
            ("SigmaUU_jet2_eta", "Subleading-jet eta", "ratio"),
            ("SigmaUU_delta_r", "Leading-dijet Delta R", "ratio"),
        ),
    },
)


class CampaignError(RuntimeError):
    """A user-facing phenomenology campaign error."""


def _load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise CampaignError(f"Could not read {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise CampaignError(f"Expected a JSON object in {path}")
    return value


def _measurement_snapshot(measurement: Mapping[str, Any]) -> dict[str, Any]:
    """Load either external reference data or an internal MC definition.

    Internal measurements deliberately have no fabricated data points or
    reference YODA.  Their checked snapshot instead fixes the selection,
    raw-object names, and binnings that define the observable contract.
    """

    reference = measurement["reference"]
    if reference.get("kind") != "internal_observable_definition":
        return validate_vendored(str(measurement["id"]))

    snapshot_path = DISPOL_ROOT / str(reference["snapshot"])
    snapshot = _load_json(snapshot_path)
    if snapshot.get("kind") != "internal_observable_definition":
        raise CampaignError(
            f"{snapshot_path} is not an internal observable definition"
        )
    if snapshot.get("measurement") != measurement["id"]:
        raise CampaignError(
            f"{snapshot_path} belongs to {snapshot.get('measurement')!r}, "
            f"not {measurement['id']!r}"
        )
    definitions = snapshot.get("observables")
    if not isinstance(definitions, dict) or not definitions:
        raise CampaignError(f"{snapshot_path} defines no observables")

    configured: dict[str, tuple[str, str]] = {}
    helicity_resolved: set[str] = set()
    for channel_id, channel in measurement["channels"].items():
        raw_objects = channel.get("raw_objects")
        if not isinstance(raw_objects, dict) or not raw_objects:
            raise CampaignError(
                f"{measurement['id']}/{channel_id} defines no raw objects"
            )
        for observable, raw_object in raw_objects.items():
            if observable in configured:
                raise CampaignError(
                    f"{measurement['id']} repeats internal observable "
                    f"{observable!r} across channels"
                )
            source_analysis, object_name = _raw_object_source(
                measurement, raw_object
            )
            configured[str(observable)] = (source_analysis, object_name)
        resolved = channel.get("helicity_resolved_observables", [])
        if (
            not isinstance(resolved, list)
            or any(not isinstance(item, str) for item in resolved)
            or len(resolved) != len(set(resolved))
        ):
            raise CampaignError(
                f"{measurement['id']}/{channel_id} has invalid "
                "helicity_resolved_observables"
            )
        unknown_resolved = sorted(set(resolved) - set(raw_objects))
        if unknown_resolved:
            raise CampaignError(
                f"{measurement['id']}/{channel_id} resolves unknown "
                f"observables: {', '.join(unknown_resolved)}"
            )
        repeated_resolved = sorted(helicity_resolved.intersection(resolved))
        if repeated_resolved:
            raise CampaignError(
                f"{measurement['id']} repeats helicity-resolved observables: "
                f"{', '.join(repeated_resolved)}"
            )
        helicity_resolved.update(resolved)
    if set(definitions) != set(configured):
        raise CampaignError(
            f"{snapshot_path} observable keys differ from the campaign descriptor"
        )
    for observable, (source_analysis, raw_object) in configured.items():
        definition = definitions[observable]
        if not isinstance(definition, dict):
            raise CampaignError(
                f"{snapshot_path} observable {observable!r} must be an object"
            )
        if definition.get("raw_object") != raw_object:
            raise CampaignError(
                f"{snapshot_path} raw object for {observable!r} differs from "
                "the campaign descriptor"
            )
        defined_analysis = str(
            definition.get("source_analysis", measurement["analysis"]["name"])
        )
        if defined_analysis != source_analysis:
            raise CampaignError(
                f"{snapshot_path} source analysis for {observable!r} differs "
                "from the campaign descriptor"
            )
        edges = definition.get("edges")
        if (
            not isinstance(edges, list)
            or len(edges) < 2
            or any(not isinstance(value, (int, float)) for value in edges)
            or any(not math.isfinite(float(value)) for value in edges)
            or any(float(high) <= float(low) for low, high in zip(edges, edges[1:]))
        ):
            raise CampaignError(
                f"{snapshot_path} has invalid bin edges for {observable!r}"
            )
        resolved_flag = definition.get("helicity_resolved", False)
        if not isinstance(resolved_flag, bool):
            raise CampaignError(
                f"{snapshot_path} has a non-boolean helicity_resolved flag "
                f"for {observable!r}"
            )
        if resolved_flag != (observable in helicity_resolved):
            raise CampaignError(
                f"{snapshot_path} helicity-resolved contract differs from "
                f"the campaign descriptor for {observable!r}"
            )
    return snapshot


def _comparison_pair(measurement: Mapping[str, Any]) -> tuple[str, str]:
    """Legacy descriptors retain the original azimuthal-spin comparison."""
    pair = measurement.get("comparison_pair", ["nominal", "shower_spin_off"])
    if (not isinstance(pair, (list, tuple)) or len(pair) != 2
            or pair[0] != "nominal" or not isinstance(pair[1], str)
            or pair[1] == "nominal"):
        raise CampaignError("comparison_pair must be [nominal, distinct control]")
    if "comparison_pair" in measurement and any(
        family not in measurement["families"] for family in pair
    ):
        raise CampaignError("comparison_pair names an unconfigured family")
    return tuple(pair)


def _comparison_label(measurement: Mapping[str, Any], kind: str) -> str:
    defaults = {"ratio": "Shower spin on / off",
                "difference": "Spin-on minus spin-off"}
    return str(measurement.get("comparison_labels", {}).get(kind, defaults[kind]))


def _validate_descriptor(measurement: Mapping[str, Any], path: Path) -> None:
    required = {"id", "schema_version", "process_kind", "analysis", "reference",
                "cards", "channels", "families", "campaign", "postprocessor",
                "scales", "pdf_ensembles"}
    missing = sorted(required-set(measurement))
    if missing:
        raise CampaignError(f"{path} is missing: {', '.join(missing)}")
    _analysis_spec_map(measurement)
    process_kind = str(measurement["process_kind"])
    schema_version = int(measurement["schema_version"])
    if schema_version not in {3, 4, 5, 6} or process_kind not in {
        "polarized_pp", "polarized_pp_jets", "polarized_sidis",
        "unpolarized_sidis",
    }:
        raise CampaignError(f"Unsupported schema/process kind in {path}")
    if path.stem != measurement["id"]:
        raise CampaignError(f"Registry filename must match id in {path}")
    helicities = set(measurement["cards"]["helicities"])
    if process_kind == "unpolarized_sidis":
        if helicities != {"00"}:
            raise CampaignError(
                f"{measurement['id']} unpolarized SIDIS must define only 00"
            )
    elif helicities != {"PP", "PM", "MP", "MM"}:
        raise CampaignError(f"{measurement['id']} must define four physical helicities")
    contributions = set(measurement["cards"]["contributions"])
    if process_kind in {"polarized_pp", "polarized_pp_jets"} and contributions != {"LO"}:
        raise CampaignError(f"{measurement['id']} RHIC hard processes must be labelled LO")
    if process_kind in {"polarized_sidis", "unpolarized_sidis"}:
        if contributions != {"POSNLO", "NEGNLO"}:
            raise CampaignError(f"{measurement['id']} SIDIS requires POSNLO and NEGNLO")
        components = set(measurement["cards"].get("target_components", {}))
        if schema_version < 6:
            if components != {"P", "N"}:
                raise CampaignError(
                    f"{measurement['id']} schema-5 SIDIS requires proton and neutron components"
                )
        elif not components or not components <= {"P", "N"}:
            raise CampaignError(
                f"{measurement['id']} schema-6 SIDIS target components must be a nonempty subset of P/N"
            )
        if schema_version == 6:
            config = measurement.get("postprocess_config")
            if not isinstance(config, dict):
                raise CampaignError(
                    f"{measurement['id']} schema-6 SIDIS requires postprocess_config"
                )
            outputs = config.get("target_outputs")
            if not isinstance(outputs, dict) or not outputs:
                raise CampaignError(
                    f"{measurement['id']} target_outputs must be a nonempty object"
                )
            for output, weights in outputs.items():
                if not str(output) or not isinstance(weights, dict) or not weights:
                    raise CampaignError(
                        f"{measurement['id']} has an invalid target output {output!r}"
                    )
                if not set(weights) <= components:
                    raise CampaignError(
                        f"{measurement['id']} target output {output!r} uses an unconfigured component"
                    )
                if any(not isinstance(weight, (int, float)) or
                       not math.isfinite(float(weight)) for weight in weights.values()):
                    raise CampaignError(
                        f"{measurement['id']} target output {output!r} has non-finite weights"
                    )
            if process_kind == "polarized_sidis":
                scale = config.get("longitudinal_target_scale")
                if not isinstance(scale, (int, float)) or not math.isfinite(float(scale)):
                    raise CampaignError(
                        f"{measurement['id']} requires a finite longitudinal_target_scale"
                    )
    if "nominal" not in measurement["families"]:
        raise CampaignError(f"{measurement['id']} has no nominal family")
    _comparison_pair(measurement)
    for family_id, family in measurement["families"].items():
        if family.get("hard_process_spin") not in {None, "on", "off"}:
            raise CampaignError("hard_process_spin must be on or off")
        shower_spin = family.get("shower_spin_correlations")
        if shower_spin not in {None, "on", "off"}:
            raise CampaignError(
                f"{measurement['id']} family {family_id!r} has invalid "
                "shower_spin_correlations; expected 'on' or 'off'"
            )
        plot_options = family.get("plot_options", {})
        if not isinstance(plot_options, dict):
            raise CampaignError(
                f"{measurement['id']} family {family_id!r} plot_options "
                "must be an object"
            )
        for key, value in plot_options.items():
            if not re.fullmatch(r"[A-Za-z][A-Za-z0-9]*", str(key)):
                raise CampaignError(
                    f"{measurement['id']} family {family_id!r} has invalid "
                    f"plot option name {key!r}"
                )
            if ":" in str(value):
                raise CampaignError(
                    f"{measurement['id']} family {family_id!r} plot option "
                    f"{key!r} cannot contain ':'"
                )
    if process_kind == "polarized_pp_jets":
        generator_cuts = measurement.get("generator_cuts", {})
        allowed = generator_cuts.get("jet_kt_min_scan_gev")
        nominal = generator_cuts.get("jet_kt_min_gev")
        if nominal not in {3.0, 4.0, 5.0} or allowed != [3.0, 4.0, 5.0]:
            raise CampaignError(
                f"{measurement['id']} must pin the 3, 4, and 5 GeV jet-cut scan"
            )
    for axis in ("polarized", "unpolarized"):
        ensemble = measurement["pdf_ensembles"][axis]
        if not isinstance(ensemble.get("active", True), bool):
            raise CampaignError(
                f"{measurement['id']} {axis} active flag must be boolean"
            )
        if int(ensemble["central_member"]) != 0 or list(ensemble["replica_members"]) != [1, 100]:
            raise CampaignError(f"{measurement['id']} {axis} ensemble is not the pinned 101-member set")


def discover_pp_registry() -> dict[str, dict[str, Any]]:
    registry: dict[str, dict[str, Any]] = {}
    for path in sorted(REGISTRY_DIR.glob("*.json")):
        value = _load_json(path)
        _validate_descriptor(value, path)
        value["_registry_path"] = str(path)
        registry[str(value["id"])] = value
    return registry


def discover_all() -> dict[str, dict[str, Any]]:
    registry = experimental.discover_registry()
    for value in registry.values():
        value["process_kind"] = "fixed_target_dis"
    for identifier, value in discover_pp_registry().items():
        if identifier in registry:
            raise CampaignError(f"Duplicate measurement id {identifier}")
        registry[identifier] = value
    return dict(sorted(registry.items()))


def _campaign_dir(measurement: str, tag: str) -> Path:
    if not TAG_PATTERN.fullmatch(tag):
        raise CampaignError(f"Invalid campaign tag {tag!r}")
    return CAMPAIGN_ROOT / measurement / tag


def _member_selector(value: str | None, profile: str) -> list[int]:
    token = value or ("all" if profile == "paper" else "central")
    if token == "central":
        return [0]
    if token == "all":
        return list(range(101))
    try:
        members = sorted({int(item) for item in token.split(",")})
    except ValueError as exc:
        raise CampaignError(f"Invalid PDF member selector {token!r}") from exc
    if not members or members[0] < 0 or members[-1] > 100:
        raise CampaignError("PDF members must be in the inclusive range 0--100")
    return members


def _scale_selector(value: str | None, profile: str) -> list[float]:
    token = value or ("all" if profile == "paper" else "central")
    if token == "central":
        return [1.0]
    if token == "all":
        return [0.5, 1.0, 2.0]
    try:
        scales = sorted({float(item) for item in token.split(",")})
    except ValueError as exc:
        raise CampaignError(f"Invalid hard-scale selector {token!r}") from exc
    if not scales or any(value not in {0.5, 1.0, 2.0} for value in scales):
        raise CampaignError("Hard scales must be selected from 0.5, 1, and 2")
    return scales


def _variation_points(
    args: argparse.Namespace,
    measurement: Mapping[str, Any] | None = None,
) -> list[tuple[int, int, float]]:
    ensembles = measurement.get("pdf_ensembles", {}) if measurement else {}
    polarized_active = bool(ensembles.get("polarized", {}).get("active", True))
    unpolarized_active = bool(ensembles.get("unpolarized", {}).get("active", True))
    if polarized_active:
        polarized = _member_selector(args.polarized_pdf_members, args.profile)
    else:
        requested = _member_selector(args.polarized_pdf_members, "central")
        if requested != [0]:
            raise CampaignError(
                "Unpolarized SIDIS accepts only the central polarized-PDF member"
            )
        polarized = [0]
    if unpolarized_active:
        unpolarized = _member_selector(args.unpolarized_pdf_members, args.profile)
    else:
        requested = _member_selector(args.unpolarized_pdf_members, "central")
        if requested != [0]:
            raise CampaignError("This descriptor accepts only the central unpolarized-PDF member")
        unpolarized = [0]
    scales = _scale_selector(args.scales, args.profile)
    # The two PDF ensembles are varied independently; a Cartesian product
    # would not represent the prescribed uncertainty construction.
    points = {(0, 0, 1.0)}
    if polarized_active:
        points.update((member, 0, 1.0) for member in polarized)
    if unpolarized_active:
        points.update((0, member, 1.0) for member in unpolarized)
    points.update((0, 0, scale) for scale in scales)
    return sorted(points, key=lambda item: (item[2] != 1.0, item[0] != 0,
                                             item[1] != 0, item))


def _scale_token(value: float) -> str:
    return {0.5: "0p5", 1.0: "1", 2.0: "2"}[float(value)]


def _analysis_spec_map(
    measurement: Mapping[str, Any],
) -> dict[str, Mapping[str, Any]]:
    try:
        specs = experimental.analysis_specs(measurement)
    except experimental.CampaignError as exc:
        raise CampaignError(str(exc)) from exc
    return {str(spec["name"]): spec for spec in specs}


def _analysis_options(
    measurement: Mapping[str, Any], family_id: str,
    spec: Mapping[str, Any],
) -> Mapping[str, Any]:
    if spec is measurement["analysis"]:
        options = measurement["families"][family_id].get(
            "analysis_options", {}
        )
    else:
        options = spec.get("analysis_options", {})
    if not isinstance(options, Mapping):
        raise CampaignError(
            f"{measurement['id']} analysis options for {spec['name']} "
            "must be an object"
        )
    return options


def _analysis_instance_for_spec(
    measurement: Mapping[str, Any], family_id: str,
    spec: Mapping[str, Any],
) -> str:
    """Return one exact Rivet analysis identifier written into YODA paths."""

    name = str(spec["name"])
    options = _analysis_options(measurement, family_id, spec)
    if not options:
        return name
    return name + ":" + ":".join(
        f"{key}={value}" for key, value in sorted(options.items())
    )


def _analysis_instances(
    measurement: Mapping[str, Any], family_id: str
) -> dict[str, str]:
    return {
        str(spec["name"]): _analysis_instance_for_spec(
            measurement, family_id, spec
        )
        for spec in experimental.analysis_specs(measurement)
    }


def _analysis_instance(
    measurement: Mapping[str, Any], family_id: str
) -> str:
    """Return the primary Rivet analysis identifier (legacy helper)."""

    return _analysis_instance_for_spec(
        measurement, family_id, measurement["analysis"]
    )


def _raw_object_source(
    measurement: Mapping[str, Any], raw_object: Any
) -> tuple[str, str]:
    """Resolve a backward-compatible raw-object descriptor."""

    primary = str(measurement["analysis"]["name"])
    if isinstance(raw_object, str):
        source_analysis, object_name = primary, raw_object
    elif isinstance(raw_object, Mapping):
        source_analysis = str(raw_object.get("analysis", primary))
        object_name = str(raw_object.get("object", ""))
    else:
        raise CampaignError(
            f"{measurement['id']} raw objects must be strings or objects"
        )
    if source_analysis not in _analysis_spec_map(measurement):
        raise CampaignError(
            f"{measurement['id']} raw object names unknown analysis "
            f"{source_analysis!r}"
        )
    if not object_name or object_name.startswith("/"):
        raise CampaignError(
            f"{measurement['id']} raw object has invalid name {object_name!r}"
        )
    return source_analysis, object_name


def _family_selector(
    value: str | None, measurement: Mapping[str, Any]
) -> list[str]:
    families = measurement["families"]
    token = value or "nominal"
    if token == "nominal":
        selected = ["nominal"]
    elif token == "all":
        selected = list(families)
    else:
        selected = [item.strip() for item in token.split(",") if item.strip()]
    unknown = sorted(set(selected)-set(families))
    if unknown:
        raise CampaignError(
            f"Unknown families for {measurement['id']}: {', '.join(unknown)}"
        )
    if not selected:
        raise CampaignError("At least one physics family must be selected")
    return selected


def _resolved_options(args: argparse.Namespace, measurement: Mapping[str, Any]) -> dict[str, Any]:
    defaults = measurement["campaign"]
    smoke = bool(args.smoke)
    process_kind = str(measurement["process_kind"])
    if process_kind in {"polarized_sidis", "unpolarized_sidis"}:
        posnlo = int(defaults["smoke_events"] if smoke else
                     (getattr(args, "posnlo_events", None) or
                      defaults["default_events"]["POSNLO"]))
        negnlo = int(defaults["smoke_events"] if smoke else
                     (getattr(args, "negnlo_events", None) or
                      defaults["default_events"]["NEGNLO"]))
        events_by_contribution = {"POSNLO": posnlo, "NEGNLO": negnlo}
        lo_events = 0
    else:
        lo_events = int(defaults["smoke_events"] if smoke else
                        (getattr(args, "lo_events", None) or
                         defaults["default_events"]["LO"]))
        events_by_contribution = {"LO": lo_events}
    requested_jet_cut = getattr(args, "jet_kt_min_gev", None)
    if process_kind == "polarized_pp_jets":
        generator_cuts = measurement["generator_cuts"]
        jet_kt_min_gev = float(
            generator_cuts["jet_kt_min_gev"]
            if requested_jet_cut is None else requested_jet_cut
        )
        allowed_cuts = {
            float(value) for value in generator_cuts["jet_kt_min_scan_gev"]
        }
        if jet_kt_min_gev not in allowed_cuts:
            raise CampaignError(
                "Polarized-pp jet generator cuts must be selected from "
                "3, 4, and 5 GeV"
            )
    else:
        if requested_jet_cut is not None:
            raise CampaignError(
                "--jet-kt-min-gev is only valid for polarized pp jet measurements"
            )
        jet_kt_min_gev = None

    selected_families = _family_selector(
        getattr(args, "families", None), measurement
    )
    if (
        getattr(args, "nominal_prediction", None) is not None
        and "nominal" in selected_families
    ):
        raise CampaignError(
            "--nominal-prediction requires a comparison-family-only campaign"
        )
    if (
        getattr(args, "nominal_prediction", None) is not None
        and not bool(getattr(args, "plot_comparisons", False))
    ):
        raise CampaignError(
            "--nominal-prediction requires --plot-comparisons"
        )
    if (
        getattr(args, "nominal_prediction", None) is not None
        and getattr(args, "seed_base", None) is None
    ):
        raise CampaignError(
            "--nominal-prediction requires an explicit --seed-base disjoint "
            "from the nominal campaign"
        )
    options = {
        "profile": args.profile,
        "families": selected_families,
        "jobs": int(args.jobs or defaults["default_jobs"]),
        "shards": int(args.shards or defaults["default_shards"]),
        "seed_base": int(args.seed_base or defaults["default_seed_base"]),
        "lo_events": lo_events,
        "events_by_contribution": events_by_contribution,
        "smoke": smoke,
        # Keep the manifest representation JSON-native so a resumed campaign
        # compares equal after the on-disk JSON has been reloaded.
        "variation_points": [
            list(point) for point in _variation_points(args, measurement)
        ],
    }
    if jet_kt_min_gev is not None:
        options["jet_kt_min_gev"] = jet_kt_min_gev
    if (options["jobs"] <= 0 or options["shards"] <= 0 or
            any(int(value) <= 0 for value in events_by_contribution.values())):
        raise CampaignError("Jobs, shards, and event counts must be positive")
    return options


def build_job_matrix(measurement: Mapping[str, Any], options: Mapping[str, Any]) -> list[dict[str, Any]]:
    families = measurement["families"]
    has_companions = bool(measurement["analysis"].get("companions"))
    selected_families = list(options.get("families", ["nominal"]))
    jobs: list[dict[str, Any]] = []
    seed_slot = 0
    for family_id in selected_families:
        family = families[family_id]
        variation_points = options["variation_points"] if family_id == "nominal" else [(0, 0, 1.0)]
        for channel in measurement["channels"]:
            for polarized_member, unpolarized_member, scale in variation_points:
                targets = (measurement["cards"].get("target_components", {"none": "none"})
                           if measurement["process_kind"] in {
                               "polarized_sidis", "unpolarized_sidis"
                           }
                           else {"none": "none"})
                contributions = measurement["cards"]["contributions"]
                for target_component in targets:
                    for helicity in family["helicities"]:
                        for contribution in contributions:
                            event_splits = experimental.split_events(
                                int(options["events_by_contribution"][contribution]),
                                int(options["shards"]),
                            )
                            level = str(family.get("observable_level", "default"))
                            analysis_instance = _analysis_instance(
                                measurement, family_id
                            )
                            analysis_instances = _analysis_instances(
                                measurement, family_id
                            )
                            logical = (
                                f"{channel}-target{target_component}-level{level}-"
                                f"{family_id}-{helicity}-{contribution}-"
                                f"p{polarized_member:03d}-u{unpolarized_member:03d}-"
                                f"mu{_scale_token(scale)}-mpi{family['mpi']}"
                            )
                            if measurement["process_kind"] == "polarized_pp_jets":
                                logical += (
                                    f"-kt{float(options['jet_kt_min_gev']):.0f}gev"
                                )
                            run_stem = f"{measurement['id']}_{logical}"
                            for shard, events in enumerate(event_splits, start=1):
                                job_id = (
                                    f"{logical}-s{shard:03d}-of-{len(event_splits):03d}"
                                )
                                job = {
                                    "id": job_id, "measurement": measurement["id"],
                                    "channel": channel, "target_component": target_component,
                                    "observable_level": level, "family": family_id,
                                    "family_label": family["label"], "helicity": helicity,
                                    "analysis_instance": analysis_instance,
                                    "contribution": contribution, "order": contribution,
                                    "polarized_pdf_member": polarized_member,
                                    "unpolarized_pdf_member": unpolarized_member,
                                    "scale": scale, "mpi": family["mpi"],
                                    "shower_spin_correlations": family.get(
                                        "shower_spin_correlations", "on"
                                    ),
                                    "jet_kt_min_gev": options.get("jet_kt_min_gev"),
                                    "shard": shard, "shards": len(event_splits),
                                    "events": events,
                                    "seed": int(options["seed_base"])+seed_slot,
                                    "initial_seed": int(options["seed_base"])+seed_slot,
                                    "attempt": 0, "status": "planned", "stem": run_stem,
                                    "card_input": f"{run_stem}.in",
                                    "run_file": f"runs/{run_stem}.run",
                                    "output_yoda": f"yoda/{job_id}.yoda",
                                }
                                if has_companions:
                                    job["analysis_instances"] = analysis_instances
                                if "hard_process_spin" in family:
                                    job["hard_process_spin"] = family["hard_process_spin"]
                                jobs.append(job)
                                seed_slot += 1
    ids = [job["id"] for job in jobs]
    seeds = [job["seed"] for job in jobs]
    if len(ids) != len(set(ids)) or len(seeds) != len(set(seeds)):
        raise CampaignError("Job identities or initial seeds are not unique")
    return jobs


def _plan(measurement: Mapping[str, Any], args: argparse.Namespace) -> dict[str, Any]:
    options = _resolved_options(args, measurement)
    jobs = build_job_matrix(measurement, options)
    logical = {(job["channel"], job["target_component"], job["observable_level"],
                job["family"], job["helicity"], job["contribution"],
                job["polarized_pdf_member"], job["unpolarized_pdf_member"],
                job["scale"], job["mpi"]) for job in jobs}
    return {"measurement": measurement["id"], "tag": args.tag,
            "campaign_directory": str(_campaign_dir(measurement["id"], args.tag)),
            "options": options, "logical_jobs": len(logical),
            "shard_jobs": len(jobs), "jobs": jobs}


def _runtime(measurement: Mapping[str, Any]) -> dict[str, Any]:
    if measurement["process_kind"] in {"polarized_sidis", "unpolarized_sidis"}:
        return experimental.preflight_runtime(measurement)
    tools: dict[str, str] = {}
    for name in (
        "Herwig", "rivet", "rivet-build", "rivet-mkhtml", "rivet-config",
        "lhapdf", "lhapdf-config",
    ):
        executable = shutil.which(name)
        if not executable:
            raise CampaignError(f"Required executable {name!r} is not active; load herwig/pol")
        tools[name] = str(Path(executable).resolve())
    prefix = Path(tools["Herwig"]).parent.parent.resolve()
    herwig_core = experimental._find_runtime_library(
        prefix, "lib/Herwig", "Herwig.so*"
    )
    hwmedis = experimental._find_runtime_library(
        prefix, "lib/Herwig", "HwMEDIS*.so*"
    )
    hwmehadron = experimental._find_runtime_library(
        prefix, "lib/Herwig", "HwMEHadron*.so*"
    )
    hwshower = experimental._find_runtime_library(
        prefix, "lib/Herwig", "HwShower*.so*"
    )
    fixed_target = experimental._find_runtime_library(
        prefix, "lib/ThePEG", "FixedTargetLuminosity*.so*"
    )
    herwig_repository = prefix / "share" / "Herwig" / "HerwigDefaults.rpo"
    for axis in ("unpolarized", "polarized"):
        ensemble = measurement["pdf_ensembles"][axis]
        experimental._command_output([tools["lhapdf"], "show", ensemble["set"]])
    compiler = ""
    if sys.platform == "darwin":
        candidates = sorted(Path("/opt/homebrew/bin").glob("g++-[0-9]*"), reverse=True)
        compiler = str(candidates[0]) if candidates else ""
    pdf_sets = [
        measurement["pdf_ensembles"][name]["set"]
        for name in ("unpolarized", "polarized")
    ]
    try:
        file_provenance = provenance.runtime_record(
            repository=DISPOL_ROOT,
            tools=tools,
            herwig_prefix=prefix,
            artifact_paths={
                "Herwig": Path(tools["Herwig"]),
                "HerwigCore": herwig_core,
                "HerwigDefaults.rpo": herwig_repository,
                "HwMEDIS": hwmedis,
                "HwMEHadron": hwmehadron,
                "HwShower": hwshower,
                "FixedTargetLuminosity": fixed_target,
                "Rivet": Path(tools["rivet"]),
            },
            pdf_sets=pdf_sets,
            lhapdf_data_directory=Path(
                experimental._command_output(
                    [tools["lhapdf-config"], "--datadir"]
                )
            ),
        )
    except provenance.ProvenanceError as exc:
        raise CampaignError(str(exc)) from exc

    return {
        "checked_at": experimental.utc_now(), "tools": tools,
        "herwig_prefix": str(prefix),
        "herwig_version": experimental._command_output([tools["Herwig"], "--version"]),
        "rivet_version": experimental._command_output([tools["rivet"], "--version"]),
        "rivet_data_directory": experimental._command_output([tools["rivet-config"], "--datadir"]),
        "herwig_core_library": str(herwig_core),
        "hwmedis_library": str(hwmedis),
        "hwmehadron_library": str(hwmehadron),
        "hwshower_library": str(hwshower),
        "fixed_target_library": str(fixed_target),
        "herwig_repository": str(herwig_repository.resolve()),
        "rivet_plugin_compiler": compiler,
        "pdf_sets": pdf_sets,
        "provenance": file_provenance,
        "environment": {key: os.environ.get(key, "") for key in
                        ("HERWIG_ENV", "RIVET_ANALYSIS_PATH", "RIVET_DATA_PATH",
                         "DYLD_LIBRARY_PATH", "LD_LIBRARY_PATH")},
    }


def _signature(
    measurement: Mapping[str, Any], *, plot_bytes: bytes | None = None
) -> str:
    digest = hashlib.sha256(experimental.canonical_json_bytes(
        {key: value for key, value in measurement.items() if not key.startswith("_")}))
    plot_path = DISPOL_ROOT / measurement["analysis"]["plot"]
    files: list[Path] = []
    for spec in experimental.analysis_specs(measurement):
        files.extend(
            DISPOL_ROOT / str(spec[key])
            for key in ("source", "info", "plot")
        )
        files.extend(
            DISPOL_ROOT / str(path)
            for path in spec.get("support_files", [])
        )
    files.append(DISPOL_ROOT/measurement["reference"]["snapshot"])
    for key in ("source_manifest", "raw_snapshot"):
        if measurement["reference"].get(key):
            files.append(DISPOL_ROOT/str(measurement["reference"][key]))
    card_dir = DISPOL_ROOT/measurement["cards"]["directory"]
    files.extend(sorted(card_dir.glob("*.in")))
    for path in files:
        digest.update(str(path.relative_to(DISPOL_ROOT)).encode())
        digest.update(
            plot_bytes if plot_bytes is not None and path == plot_path
            else path.read_bytes()
        )
    return digest.hexdigest()


def _assert_manifest_signature_current(
    manifest: Mapping[str, Any], measurement: Mapping[str, Any], *,
    allow_plot_metadata_refresh: bool = False,
) -> dict[str, Any] | None:
    """Refuse to reinterpret shards produced by another measurement definition."""

    recorded = manifest.get("configuration", {}).get("measurement_signature")
    current = _signature(measurement)
    if recorded != current:
        if allow_plot_metadata_refresh:
            try:
                return experimental.authorize_plot_metadata_refresh(
                    manifest,
                    measurement,
                    current_signature=current,
                    signature_with_plot_bytes=lambda payload: _signature(
                        measurement, plot_bytes=payload
                    ),
                )
            except experimental.CampaignError as exc:
                raise CampaignError(str(exc)) from exc
        raise CampaignError(
            f"{measurement['id']} campaign products were generated with a "
            "different analysis, reference, card, or metadata signature. "
            "Postprocessing and plotting are intentionally refused; use a "
            "new immutable tag and rerun the event campaign."
        )
    return None


def _card_text(measurement: Mapping[str, Any], job: Mapping[str, Any]) -> str:
    cards = measurement["cards"]
    source_stem = cards["stem_pattern"].format(
        channel=job["channel"], helicity=job["helicity"],
        target=job.get("target_component", "none"),
        contribution=job["contribution"], order=job["contribution"],
    )
    source = DISPOL_ROOT/cards["directory"]/f"{source_stem}.in"
    if not source.is_file() and job["helicity"] == "00":
        source_stem = cards["stem_pattern"].format(
            channel=job["channel"], helicity="PP",
            target=job.get("target_component", "none"),
            contribution=job["contribution"], order=job["contribution"],
        )
        source = DISPOL_ROOT/cards["directory"]/f"{source_stem}.in"
    if not source.is_file():
        raise CampaignError(f"Missing base card {source}")
    text = source.read_text(encoding="utf-8")
    polarization = cards["helicities"].get(job["helicity"], [0, 0])
    text = re.sub(r"(FirstLongitudinalPolarization )[-0-9]+", rf"\g<1>{polarization[0]}", text)
    text = re.sub(r"(SecondLongitudinalPolarization )[-0-9]+", rf"\g<1>{polarization[1]}", text)
    overrides = [
        f"set /Herwig/Partons/HardLOPDF:Member {job['unpolarized_pdf_member']}",
        f"set /Herwig/Partons/HardNLOPDF:Member {job['unpolarized_pdf_member']}",
        f"set /Herwig/Partons/ShowerLOPDF:Member {job['unpolarized_pdf_member']}",
        f"set /Herwig/Partons/ShowerNLOPDF:Member {job['unpolarized_pdf_member']}",
        f"set /Herwig/Partons/MPIPDF:Member {job['unpolarized_pdf_member']}",
        f"set {cards['polarized_pdf_object']}:Member {job['polarized_pdf_member']}",
    ]
    if cards.get("unpolarized_pdf_object"):
        overrides.append(
            f"set {cards['unpolarized_pdf_object']}:Member "
            f"{job['unpolarized_pdf_member']}"
        )
    family = measurement["families"][job["family"]]
    hard_spin = family.get("hard_process_spin")
    if hard_spin is not None:
        overrides.append("set /Herwig/Shower/ShowerHandler:HardProcessSpin "
                         + ("Yes" if hard_spin == "on" else "No"))
    shower_spin = family.get("shower_spin_correlations")
    if shower_spin is not None:
        overrides.append(
            "set /Herwig/Shower/ShowerHandler:SpinCorrelations "
            + ("Yes" if shower_spin == "on" else "No")
        )
    if not math.isclose(float(job["scale"]), 1.0):
        if measurement["process_kind"] in {"polarized_sidis", "unpolarized_sidis"}:
            overrides.extend(
                f"set {scale_object}:ScaleFactor {float(job['scale']):.8g}"
                for scale_object in cards["scale_objects"]
            )
        else:
            scale_object = measurement["channels"][job["channel"]]["scale_object"]
            overrides.append(
                f"set {scale_object}:ScalePreFactor "
                f"{float(job['scale'])**2:.8g}"
            )
    if job["mpi"] == "off":
        overrides.extend([
            "set /Herwig/Shower/ShowerHandler:MPIHandler NULL",
            "set /Herwig/Shower/PowhegShowerHandler:MPIHandler NULL",
            "set /Herwig/DipoleShower/DipoleShowerHandler:MPIHandler NULL",
        ])
    else:
        overrides.extend([
            "set /Herwig/Shower/ShowerHandler:MPIHandler /Herwig/UnderlyingEvent/MPIHandler",
            "set /Herwig/Shower/PowhegShowerHandler:MPIHandler /Herwig/UnderlyingEvent/MPIHandler",
            "set /Herwig/DipoleShower/DipoleShowerHandler:MPIHandler /Herwig/UnderlyingEvent/MPIHandler",
        ])
    if measurement["process_kind"] == "polarized_pp_jets":
        jet_kt_min_gev = float(job["jet_kt_min_gev"])
        overrides.append(
            f"set /Herwig/Cuts/JetKtCut:MinKT {jet_kt_min_gev:.1f}*GeV"
        )
    specs = experimental.analysis_specs(measurement)
    for spec in specs:
        analysis_name = str(spec["name"])
        analysis_instance = _analysis_instance_for_spec(
            measurement, str(job["family"]), spec
        )
        text, analysis_replacements = re.subn(
            rf"(insert\s+/Herwig/Analysis/Rivet:Analyses\s+\d+\s+)"
            rf"{re.escape(analysis_name)}(?::\S+)?",
            rf"\g<1>{analysis_instance}",
            text,
        )
        if (
            (len(specs) > 1 or _analysis_options(
                measurement, str(job["family"]), spec
            ))
            and analysis_replacements != 1
        ):
            common_path = (
                DISPOL_ROOT / cards["directory"] / cards["common"]
            )
            common_text = common_path.read_text(encoding="utf-8")
            inherited = re.search(
                rf"insert\s+/Herwig/Analysis/Rivet:Analyses\s+\d+\s+"
                rf"{re.escape(analysis_name)}(?::\S+)?(?:\s|$)",
                common_text,
            )
            if analysis_replacements != 0 or inherited is None:
                raise CampaignError(
                    f"Expected one {analysis_name} Rivet insertion in "
                    f"{source} or its common card {common_path}"
                )
    saverun = f"saverun {job['stem']} EventGenerator"
    text, replacements = re.subn(
        r"saverun\s+\S+\s+EventGenerator",
        "\n".join(overrides+[saverun]),
        text,
    )
    if replacements != 1:
        raise CampaignError(f"Expected one saverun command in {source}")
    return text


def _manifest_configuration(measurement: Mapping[str, Any], args: argparse.Namespace,
                            plan: Mapping[str, Any]) -> dict[str, Any]:
    configuration = {
        "measurement": measurement["id"], "tag": args.tag,
        "measurement_signature": _signature(measurement), **plan["options"],
    }
    if measurement["analysis"].get("companions"):
        configuration["analysis_instances"] = {
            family_id: _analysis_instances(measurement, family_id)
            for family_id in plan["options"]["families"]
        }
    if "comparison_pair" in measurement:
        configuration["comparison_pair"] = list(_comparison_pair(measurement))
        configuration["shower_spin_policy"] = {
            family: {key: measurement["families"][family][key]
                     for key in ("hard_process_spin", "shower_spin_correlations")}
            for family in plan["options"]["families"]
        }
    return configuration


def _preflight_hard_process_spin(
    measurement: Mapping[str, Any], runtime: dict[str, Any]
) -> None:
    if not any("hard_process_spin" in family
               for family in measurement["families"].values()):
        return
    with tempfile.TemporaryDirectory(prefix="herwig-hard-spin-preflight-") as directory:
        card = Path(directory) / "interface.in"
        card.write_text(
            "library HwMEHadron.so\n"
            "set /Herwig/Shower/ShowerHandler:HardProcessSpin Yes\n"
            "set /Herwig/Shower/ShowerHandler:HardProcessSpin No\n"
            "get /Herwig/Shower/ShowerHandler:HardProcessSpin\n",
            encoding="utf-8",
        )
        environment = dict(os.environ, LD_DEBUG="libs", DYLD_PRINT_LIBRARIES="1")
        result = subprocess.run(
            [runtime["tools"]["Herwig"], "read", str(card)], cwd=directory,
            capture_output=True, text=True, check=False, env=environment,
        )
        output = result.stdout + result.stderr
        if result.returncode != 0 or not re.search(r"\bNo\b", output):
            raise CampaignError(
                "The active Herwig does not provide the required HardProcessSpin "
                "interface. Load the validated hard-spin runtime.\n" + output
            )
        loaded = set(re.findall(r"calling init: (/[^\n]+)", output))
        loaded.update(re.findall(r"dyld\[\d+\]: <[^>]+> (/[^\n]+)", output))
        resolved = {Path(path.strip()).resolve() for path in loaded}
        for label, key in (("HwShower", "hwshower_library"),
                           ("HwMEHadron", "hwmehadron_library")):
            if key not in runtime or Path(runtime[key]).resolve() not in resolved:
                raise CampaignError(f"Could not verify the actually loaded {label} library")
        runtime["loaded_library_fingerprints"] = {
            str(path): provenance.file_record(path)
            for path in sorted(resolved)
            if path.is_file() and ("Herwig" in str(path) or "ThePEG" in str(path))
        }


def prepare_pp(args: argparse.Namespace, measurement: Mapping[str, Any]) -> Path:
    plan = _plan(measurement, args)
    campaign_dir = _campaign_dir(measurement["id"], args.tag)
    if args.dry_run:
        print(json.dumps(plan, indent=2, sort_keys=True))
        return campaign_dir
    runtime = _runtime(measurement)
    _preflight_hard_process_spin(measurement, runtime)
    configuration = _manifest_configuration(measurement, args, plan)
    manifest_path = campaign_dir/experimental.MANIFEST_NAME
    existing = _load_json(manifest_path) if manifest_path.exists() else None
    if existing is not None and existing.get("configuration") != configuration:
        raise CampaignError("An incompatible manifest exists; use a new immutable tag")
    for directory in ("build", "cards", "runs", "yoda", "logs", "work", "postprocess", "plots"):
        (campaign_dir/directory).mkdir(parents=True, exist_ok=True)
    snapshot = _measurement_snapshot(measurement)
    if measurement["reference"].get("kind") != "internal_observable_definition":
        reference_yoda_name = measurement["analysis"].get("reference_yoda")
        if not reference_yoda_name:
            raise CampaignError(
                f"{measurement['id']} external measurement has no reference YODA"
            )
        reference_yoda = DISPOL_ROOT / str(reference_yoda_name)
        if not reference_yoda.is_file() or reference_yoda.stat().st_size == 0:
            write_reference_yoda(str(measurement["id"]), snapshot)
    common_source = DISPOL_ROOT/measurement["cards"]["directory"]/measurement["cards"]["common"]
    shutil.copy2(common_source, campaign_dir/"cards"/common_source.name)
    logical_cards: dict[str, Mapping[str, Any]] = {}
    for job in plan["jobs"]:
        logical_cards.setdefault(str(job["stem"]), job)
    for stem, job in logical_cards.items():
        destination = campaign_dir/"cards"/f"{stem}.in"
        generated = _card_text(measurement, job)
        if destination.exists() and destination.read_text(encoding="utf-8") != generated:
            raise CampaignError(f"Refusing to overwrite changed generated card {destination}")
        if not destination.exists():
            experimental.atomic_write_text(destination, generated)
    plugin = experimental.build_rivet_plugin(measurement, campaign_dir, runtime)
    environment = experimental.analysis_environment(measurement, campaign_dir, runtime)
    manifest = existing or {
        "manifest_version": 2, "measurement": measurement["id"], "tag": args.tag,
        "created_at": experimental.utc_now(), "configuration": configuration,
        "jobs": plan["jobs"], "history": [],
    }
    manifest.update({"status": "preparing", "updated_at": experimental.utc_now(),
                     "runtime": runtime, "plugin": str(plugin.relative_to(campaign_dir)),
                     "plugin_provenance": provenance.file_record(plugin),
                     "reference_snapshot": {"path": measurement["reference"]["snapshot"],
                     "sha256": experimental.sha256_file(DISPOL_ROOT/measurement["reference"]["snapshot"])}})
    manifest["history"].append({"at": experimental.utc_now(), "action": "prepare"})
    experimental.atomic_write_json(manifest_path, manifest)
    for stem in sorted(logical_cards):
        destination = campaign_dir/"runs"/f"{stem}.run"
        if destination.is_file() and destination.stat().st_size > 0:
            continue
        experimental._run_logged([runtime["tools"]["Herwig"], "read", f"{stem}.in"],
                                 campaign_dir/"cards", environment,
                                 campaign_dir/"logs"/f"read-{stem}.log")
        generated = campaign_dir/"cards"/f"{stem}.run"
        if not generated.is_file() or generated.stat().st_size == 0:
            raise CampaignError(f"Herwig read did not create {generated}")
        os.replace(generated, destination)
    prepared_records = {
        str(path.relative_to(campaign_dir)): provenance.file_record(path)
        for path in sorted(
            [plugin]
            + [path for path in (campaign_dir / "cards").rglob("*.in")]
            + [path for path in (campaign_dir / "runs").glob("*.run")]
        )
    }
    manifest["prepared_artifacts"] = {
        "files": prepared_records,
        "inventory_sha256": provenance.inventory_digest(prepared_records),
    }
    manifest["status"] = "prepared"
    manifest["updated_at"] = experimental.utc_now()
    manifest["history"].append({"at": experimental.utc_now(), "action": "prepared"})
    experimental.atomic_write_json(manifest_path, manifest)
    experimental.atomic_write_json(campaign_dir/"resolved-measurement.json",
        {key: value for key, value in measurement.items() if not key.startswith("_")})
    print(f"Prepared {measurement['id']} at {campaign_dir}")
    return campaign_dir


def run_pp(args: argparse.Namespace, measurement: Mapping[str, Any]) -> Path:
    plan = _plan(measurement, args)
    campaign_dir = _campaign_dir(measurement["id"], args.tag)
    manifest_path = campaign_dir/experimental.MANIFEST_NAME
    if args.dry_run and not manifest_path.exists():
        print(json.dumps(plan, indent=2, sort_keys=True))
        return campaign_dir
    if not manifest_path.exists():
        prepare_pp(args, measurement)
    manifest = _load_json(manifest_path)
    expected = _manifest_configuration(measurement, args, plan)
    if manifest.get("configuration") != expected:
        raise CampaignError("Campaign options differ from the immutable manifest")
    scheduled, blocked = experimental.pending_jobs(
        manifest, campaign_dir, bool(args.recover_failed))
    if blocked:
        raise CampaignError("Failed shards require --recover-failed: " +
                            ", ".join(str(job["id"]) for job in blocked))
    if args.dry_run:
        print(json.dumps({"scheduled": scheduled,
                          "already_complete": len(manifest["jobs"])-len(scheduled)},
                         indent=2, sort_keys=True))
        return campaign_dir
    if not scheduled:
        print(f"Campaign {measurement['id']}/{args.tag} is already complete")
        return campaign_dir
    runtime = manifest.get("runtime") or _runtime(measurement)
    started = float(getattr(args, "_tracker_started_at", time.time()))
    max_workers = max(1, int(plan["options"]["jobs"]))
    max_listed = max(0, int(args.max_listed))
    interval = float(args.progress_interval)
    active: dict[str, experimental.JobActivity] = {}
    active_lock = threading.Lock()
    failures: list[experimental.JobResult] = []
    by_id = {str(job["id"]): job for job in manifest["jobs"]}
    manifest["status"] = "running"
    manifest["updated_at"] = experimental.utc_now()
    manifest["history"].append({"at": experimental.utc_now(), "action": "campaign",
                                "scheduled": len(scheduled)})
    experimental.atomic_write_json(manifest_path, manifest)

    def snapshot(remaining: set[concurrent.futures.Future[experimental.JobResult]],
                 futures: Mapping[concurrent.futures.Future[experimental.JobResult], Mapping[str, Any]],
                 phase: str = "running-herwig", message: str | None = None) -> dict[str, Any]:
        with active_lock:
            running = list(active.values())
        active_ids = {str(item.job["id"]) for item in running}
        pending = [futures[future] for future in remaining
                   if str(futures[future]["id"]) not in active_ids]
        herwig = experimental.build_herwig_progress_payload(
            manifest["jobs"], pending, running, max_listed, started)
        payload = experimental.build_campaign_monitor_payload(
            measurement=measurement, tag=args.tag, phase=phase,
            started_at=started, message=message, herwig=herwig)
        experimental.write_campaign_monitor_files(campaign_dir, payload)
        return payload

    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(experimental._run_tracked_campaign_job, job, measurement,
                            campaign_dir, runtime, args.tag, active, active_lock): job
            for job in scheduled
        }
        remaining = set(futures)
        payload = snapshot(remaining, futures)
        experimental.emit_progress(experimental.progress_lines(args.tag, payload["herwig"]),
                                   sys.stdout.isatty())
        last_refresh = time.monotonic()
        while remaining:
            timeout = None if interval < 0 else max(0.05, max(1.0, interval)-
                                                     (time.monotonic()-last_refresh))
            done, _ = concurrent.futures.wait(
                remaining, timeout=timeout,
                return_when=concurrent.futures.FIRST_COMPLETED)
            for future in done:
                spec = futures[future]
                try:
                    result = future.result()
                except Exception as exc:
                    result = experimental.JobResult(str(spec["id"]), False, 1, 0,
                                                    f"Campaign worker failed: {exc}")
                remaining.remove(future)
                job = by_id[result.job_id]
                job.update({"finished_at": experimental.utc_now(),
                            "returncode": result.returncode, "message": result.message,
                            "output_size": result.output_size,
                            "status": "success" if result.success else "failed"})
                if not result.success:
                    failures.append(result)
                manifest["updated_at"] = experimental.utc_now()
                experimental.atomic_write_json(manifest_path, manifest)
            if interval >= 0 and remaining and time.monotonic()-last_refresh >= max(1.0, interval):
                payload = snapshot(remaining, futures)
                experimental.emit_progress(
                    experimental.progress_lines(args.tag, payload["herwig"]),
                    sys.stdout.isatty())
                last_refresh = time.monotonic()
    manifest["status"] = "failed" if failures else "complete"
    manifest["updated_at"] = experimental.utc_now()
    manifest["history"].append({"at": experimental.utc_now(),
                                "action": "campaign-finished",
                                "failures": len(failures)})
    experimental.atomic_write_json(manifest_path, manifest)
    final = snapshot(set(), {}, "campaign-failed" if failures else "campaign-complete",
                     f"{len(failures)} shard(s) failed" if failures else
                     f"Completed {len(scheduled)} shard job(s).")
    experimental.emit_progress(experimental.progress_lines(args.tag, final["herwig"]),
                               sys.stdout.isatty())
    if failures:
        raise CampaignError("One or more shards failed; inspect logs and use --recover-failed")
    return campaign_dir


def _logical_groups(manifest: Mapping[str, Any], campaign_dir: Path) -> dict[tuple[Any, ...], dict[str, list[Mapping[str, Any]]]]:
    mutable = dict(manifest)
    mutable["jobs"] = [dict(job) for job in manifest["jobs"]]
    experimental._reconcile_manifest_outputs(mutable, campaign_dir)
    incomplete = [str(job["id"]) for job in mutable["jobs"]
                  if job.get("status") != "success" or
                  not experimental._nonempty(campaign_dir/str(job["output_yoda"]))]
    if incomplete:
        raise CampaignError("Refusing postprocessing; incomplete components: " +
                            ", ".join(incomplete))
    groups: dict[tuple[Any, ...], dict[str, list[Mapping[str, Any]]]] = {}
    for job in mutable["jobs"]:
        key = (job["family"], job["channel"], int(job["polarized_pdf_member"]),
               int(job["unpolarized_pdf_member"]), float(job["scale"]), job["mpi"])
        groups.setdefault(key, {}).setdefault(str(job["helicity"]), []).append(job)
    for helicities in groups.values():
        for jobs in helicities.values():
            jobs.sort(key=lambda job: int(job["shard"]))
    return groups


def _load_series(jobs: Sequence[Mapping[str, Any]], campaign_dir: Path,
                 analysis: str, object_name: str) -> experimental.BinSeries:
    return _load_series_many(
        jobs, campaign_dir, {"histogram": (analysis, object_name)}
    )["histogram"]


def _load_series_many(
    jobs: Sequence[Mapping[str, Any]],
    campaign_dir: Path,
    requested: Mapping[str, tuple[str, str]],
) -> dict[str, experimental.BinSeries]:
    """Combine many shard objects after reading every shard YODA only once."""

    if not jobs or not requested:
        raise CampaignError("Multi-object shard loading requires jobs and objects")
    analysis_instances: set[str] = set()
    instances_by_analysis: dict[str, str] = {}
    source_analyses = {analysis for analysis, _ in requested.values()}
    for analysis in source_analyses:
        analysis_instances.clear()
        for job in jobs:
            configured = job.get("analysis_instances")
            if isinstance(configured, Mapping):
                instance = configured.get(analysis)
                if instance is None:
                    raise CampaignError(
                        f"Job {job.get('id')} has no Rivet instance for {analysis}"
                    )
                analysis_instances.add(str(instance))
            else:
                primary = str(job.get("analysis_instance", analysis))
                primary_name = primary.split(":", 1)[0]
                if analysis != primary_name:
                    raise CampaignError(
                        f"Legacy job {job.get('id')} has no companion analysis "
                        f"instance for {analysis}"
                    )
                analysis_instances.add(primary)
        if len(analysis_instances) != 1:
            raise CampaignError(
                f"Shards disagree on the {analysis} Rivet instance: "
                f"{sorted(analysis_instances)}"
            )
        instances_by_analysis[analysis] = next(iter(analysis_instances))

    object_paths = {
        label: f"/{instances_by_analysis[analysis]}/{object_name}"
        for label, (analysis, object_name) in requested.items()
    }
    shards: dict[str, list[experimental.BinSeries]] = {
        label: [] for label in requested
    }
    for job in jobs:
        loaded = experimental.read_histogram_series_many(
            campaign_dir / str(job["output_yoda"]), object_paths
        )
        for label, item in loaded.items():
            shards[label].append(item)
    event_counts = [int(job["events"]) for job in jobs]
    return {
        label: experimental.combine_shard_series(items, event_counts)
        for label, items in shards.items()
    }


def _logical_sidis_groups(
    manifest: Mapping[str, Any],
    campaign_dir: Path,
    measurement: Mapping[str, Any],
) -> dict[
    tuple[Any, ...],
    dict[str, dict[str, list[Mapping[str, Any]]]],
]:
    """Group complete SIDIS shards by target, helicity, and NLO component."""

    mutable = dict(manifest)
    mutable["jobs"] = [dict(job) for job in manifest["jobs"]]
    experimental._reconcile_manifest_outputs(mutable, campaign_dir)
    incomplete = [
        str(job["id"])
        for job in mutable["jobs"]
        if job.get("status") != "success"
        or not experimental._nonempty(
            campaign_dir / str(job["output_yoda"])
        )
    ]
    if incomplete:
        raise CampaignError(
            "Refusing postprocessing; incomplete components: "
            + ", ".join(incomplete)
        )
    groups: dict[
        tuple[Any, ...],
        dict[str, dict[str, list[Mapping[str, Any]]]],
    ] = {}
    for job in mutable["jobs"]:
        key = (
            job["family"],
            job["channel"],
            job["target_component"],
            int(job["polarized_pdf_member"]),
            int(job["unpolarized_pdf_member"]),
            float(job["scale"]),
            job["mpi"],
        )
        groups.setdefault(key, {}).setdefault(
            str(job["helicity"]), {}
        ).setdefault(str(job["contribution"]), []).append(job)
    for key, helicities in groups.items():
        family = str(key[0])
        expected_helicities = set(measurement["families"][family]["helicities"])
        if set(helicities) != expected_helicities:
            raise CampaignError(
                f"SIDIS postprocessing requires {sorted(expected_helicities)} "
                f"for family {family}"
            )
        for contributions in helicities.values():
            if set(contributions) != {"POSNLO", "NEGNLO"}:
                raise CampaignError(
                    "SIDIS postprocessing requires POSNLO and NEGNLO"
                )
            for jobs in contributions.values():
                jobs.sort(key=lambda job: int(job["shard"]))
    return groups


def _load_sidis_object(
    group: Mapping[str, Mapping[str, Sequence[Mapping[str, Any]]]],
    campaign_dir: Path,
    analysis: str,
    object_name: str,
    cache: dict[tuple[str, str, str], experimental.BinSeries] | None = None,
    helicities: Sequence[str] | None = None,
) -> dict[str, experimental.BinSeries]:
    """Load shards and add normalized NLO components for every helicity."""

    output: dict[str, experimental.BinSeries] = {}
    selected_helicities = tuple(helicities or group.keys())
    if set(selected_helicities) != set(group):
        raise CampaignError("SIDIS raw-object helicities do not match the job group")
    for helicity in selected_helicities:
        by_contribution: dict[str, experimental.BinSeries] = {}
        for contribution in ("POSNLO", "NEGNLO"):
            jobs = group[helicity][contribution]
            cache_key = (
                str(jobs[0]["id"]).rsplit("-s", 1)[0],
                contribution,
                object_name,
            )
            if cache is not None and cache_key in cache:
                by_contribution[contribution] = cache[cache_key]
                continue
            loaded = _load_series(
                jobs, campaign_dir, analysis, object_name
            )
            if cache is not None:
                cache[cache_key] = loaded
            by_contribution[contribution] = loaded
        output[helicity] = sidis.add_nlo_components(by_contribution)
    return output


def _sidis_variation_key(group_key: tuple[Any, ...]) -> tuple[Any, ...]:
    family, channel, _target, polarized, unpolarized, scale, mpi = group_key
    return family, channel, polarized, unpolarized, scale, mpi


def _sidis_component_samples(
    groups: Mapping[
        tuple[Any, ...],
        Mapping[str, Mapping[str, Sequence[Mapping[str, Any]]]],
    ],
    variation: tuple[Any, ...],
    campaign_dir: Path,
    analysis: str,
    object_name: str,
    target: str,
    cache: dict[tuple[str, str, str], experimental.BinSeries] | None = None,
) -> dict[str, experimental.BinSeries]:
    """Load the proton/neutron target components needed by one observable."""

    family, channel, polarized, unpolarized, scale, mpi = variation
    components = ("P",) if target == "proton" else ("P", "N")
    samples: dict[str, experimental.BinSeries] = {}
    for component in components:
        key = (
            family,
            channel,
            component,
            polarized,
            unpolarized,
            scale,
            mpi,
        )
        if key not in groups:
            raise CampaignError(
                f"Missing SIDIS target component {component} for "
                f"{_variation_id(variation)}"
            )
        component_samples = _load_sidis_object(
            groups[key], campaign_dir, analysis, object_name, cache
        )
        samples.update(
            {
                f"{component}:{helicity}": series
                for helicity, series in component_samples.items()
            }
        )
    return samples


def _linear_value(samples: Mapping[str, experimental.BinSeries],
                  coefficients: Mapping[str, float], index: int) -> tuple[float, float]:
    value = sum(coefficients[label]*samples[label].values[index] for label in coefficients)
    variance = sum(coefficients[label]**2*samples[label].variances[index]
                   for label in coefficients)
    return value, variance


def _ratio_at(samples: Mapping[str, experimental.BinSeries],
              numerator: Mapping[str, float], denominator: Mapping[str, float],
              index: int) -> tuple[float | None, float | None]:
    num, numvar = _linear_value(samples, numerator, index)
    den, denvar = _linear_value(samples, denominator, index)
    covariance = sum(numerator[label]*denominator[label]*samples[label].variances[index]
                     for label in numerator)
    return experimental.ratio_with_covariance(num, numvar, den, denvar, covariance)


def _ratio_arrays(samples: Mapping[str, experimental.BinSeries],
                  numerator: Mapping[str, float], denominator: Mapping[str, float]) -> tuple[list[float | None], list[float | None]]:
    values, errors = [], []
    for index in range(len(next(iter(samples.values())).values)):
        value, error = _ratio_at(samples, numerator, denominator, index)
        values.append(value); errors.append(error)
    return values, errors


DENOMINATOR = {"PP": 1.0, "PM": 1.0, "MP": 1.0, "MM": 1.0}
AL_A_NUMERATOR = {"PP": 1.0, "PM": 1.0, "MP": -1.0, "MM": -1.0}
AL_B_NUMERATOR = {"PP": 1.0, "PM": -1.0, "MP": 1.0, "MM": -1.0}
ALL_NUMERATOR = {"PP": 1.0, "PM": -1.0, "MP": -1.0, "MM": 1.0}
DELTA_SIGMA_LL_COEFFICIENTS = {
    "PP": 0.25, "PM": -0.25, "MP": -0.25, "MM": 0.25,
}


def _fold_star_all(samples: Mapping[str, experimental.BinSeries]) -> Mapping[str, experimental.BinSeries]:
    swap = {"PP": "PP", "PM": "MP", "MP": "PM", "MM": "MM"}
    folded: dict[str, experimental.BinSeries] = {}
    for helicity in DENOMINATOR:
        values, variances = [], []
        for abs_index in range(3):
            positive = 3+abs_index
            negative = 2-abs_index
            values.append(samples[helicity].values[positive] +
                          samples[swap[helicity]].values[negative])
            variances.append(samples[helicity].variances[positive] +
                             samples[swap[helicity]].variances[negative])
        folded[helicity] = experimental.BinSeries([0.0, 0.5, 1.0, 1.5],
                                                   values, variances)
    return folded


def _beam_exchange_closures(samples: Mapping[str, experimental.BinSeries]) -> dict[str, dict[str, Any]]:
    swap = {"PP": "PP", "PM": "MP", "MP": "PM", "MM": "MM"}
    edges = next(iter(samples.values())).edges
    outputs: dict[str, dict[str, Any]] = {}
    for helicity in DENOMINATOR:
        values, errors = [], []
        for index in range(len(edges)-1):
            mirror = len(edges)-2-index
            value, error = experimental.parity_residual(
                samples[helicity].values[index], samples[helicity].variances[index],
                samples[swap[helicity]].values[mirror],
                samples[swap[helicity]].variances[mirror])
            values.append(value); errors.append(error)
        outputs[helicity] = {"edges": list(edges), "values": values, "errors": errors}
    return outputs


def _star_prediction(channel: str, samples: Mapping[str, experimental.BinSeries]) -> dict[str, dict[str, Any]]:
    edges = list(next(iter(samples.values())).edges)
    output: dict[str, dict[str, Any]] = {}
    if channel in {"Wplus", "Wminus"}:
        al_a, err_a = _ratio_arrays(samples, AL_A_NUMERATOR, DENOMINATOR)
        al_b, err_b = _ratio_arrays(samples, AL_B_NUMERATOR, DENOMINATOR)
        values, errors = [], []
        for index in range(len(edges)-1):
            mirror = len(edges)-2-index
            if al_a[index] is None or al_b[mirror] is None:
                values.append(None); errors.append(None)
            else:
                values.append(0.5*(float(al_a[index])+float(al_b[mirror])))
                errors.append(0.5*math.hypot(float(err_a[index]), float(err_b[mirror])))
        output[f"{channel}_AL"] = {"edges": edges, "values": values, "errors": errors}
        folded = _fold_star_all(samples)
        all_values, all_errors = _ratio_arrays(folded, ALL_NUMERATOR, DENOMINATOR)
        output[f"{channel}_ALL"] = {"edges": [0.0,0.5,1.0,1.5],
                                     "values": all_values, "errors": all_errors}
        for helicity, closure in _beam_exchange_closures(samples).items():
            output[f"BeamExchange_{channel}_{helicity}"] = closure
    elif channel == "Zgamma":
        # Averaging the two beam asymmetries in a symmetric integrated
        # acceptance reduces algebraically to (PP-MM)/sum.
        numerator = {"PP": 1.0, "PM": 0.0, "MP": 0.0, "MM": -1.0}
        values, errors = _ratio_arrays(samples, numerator, DENOMINATOR)
        output["Zgamma_AL"] = {"edges": edges, "values": values, "errors": errors}
    else:
        raise CampaignError(f"Unknown STAR channel {channel}")
    return output


def _phenix_prediction(samples_by_object: Mapping[str, Mapping[str, experimental.BinSeries]]) -> dict[str, dict[str, Any]]:
    output: dict[str, dict[str, Any]] = {}
    for observable in ("inclusive_cross_section", "isolated_cross_section"):
        samples = samples_by_object[observable]
        series = experimental.linear_combine_series(samples,
            {label: 0.25 for label in DENOMINATOR})
        output[observable] = {
            "edges": list(series.edges), "values": list(series.values),
            "errors": [math.sqrt(max(0.0, value)) for value in series.variances],
        }
    samples = samples_by_object["isolated_all"]
    values, errors = _ratio_arrays(samples, ALL_NUMERATOR, DENOMINATOR)
    edges = list(next(iter(samples.values())).edges)
    output["isolated_all"] = {"edges": edges, "values": values, "errors": errors}
    for name, numerator in (("SingleSpinA", AL_A_NUMERATOR),
                            ("SingleSpinB", AL_B_NUMERATOR)):
        vals, errs = _ratio_arrays(samples, numerator, DENOMINATOR)
        output[f"{name}_isolated_all"] = {"edges": edges, "values": vals, "errors": errs}
    for first, second, name in (("PP", "MM", "Parity_PP_MM"),
                                ("PM", "MP", "Parity_PM_MP")):
        vals, errs = [], []
        for index in range(len(edges)-1):
            value, error = experimental.parity_residual(
                samples[first].values[index], samples[first].variances[index],
                samples[second].values[index], samples[second].variances[index])
            vals.append(value); errs.append(error)
        output[f"{name}_isolated_all"] = {"edges": edges, "values": vals, "errors": errs}
    return output


def _star_jet_prediction(
    samples_by_object: Mapping[str, Mapping[str, experimental.BinSeries]]
) -> dict[str, dict[str, Any]]:
    output: dict[str, dict[str, Any]] = {}
    for observable, samples in samples_by_object.items():
        sigma_uu = experimental.linear_combine_series(
            samples, {label: 0.25 for label in DENOMINATOR}
        )
        output[f"SigmaUU_{observable}"] = {
            "edges": list(sigma_uu.edges),
            "values": list(sigma_uu.values),
            "errors": [
                math.sqrt(max(0.0, variance))
                for variance in sigma_uu.variances
            ],
        }
        values, errors = _ratio_arrays(samples, ALL_NUMERATOR, DENOMINATOR)
        edges = list(next(iter(samples.values())).edges)
        output[observable] = {"edges": edges, "values": values, "errors": errors}
        for first, second, label in (
            ("PP", "MM", "Parity_PP_MM"),
            ("PM", "MP", "Parity_PM_MP"),
        ):
            closure_values, closure_errors = [], []
            for index in range(len(edges)-1):
                value, error = experimental.parity_residual(
                    samples[first].values[index], samples[first].variances[index],
                    samples[second].values[index], samples[second].variances[index],
                )
                closure_values.append(value)
                closure_errors.append(error)
            output[f"{label}_{observable}"] = {
                "edges": edges, "values": closure_values,
                "errors": closure_errors,
            }
        for label, numerator in (
            ("SingleSpinA", AL_A_NUMERATOR),
            ("SingleSpinB", AL_B_NUMERATOR),
        ):
            closure_values, closure_errors = _ratio_arrays(
                samples, numerator, DENOMINATOR
            )
            output[f"{label}_{observable}"] = {
                "edges": edges, "values": closure_values,
                "errors": closure_errors,
            }
    return output


def _mc_poldijets_prediction(
    samples_by_object: Mapping[str, Mapping[str, experimental.BinSeries]],
    helicity_resolved_observables: Sequence[str] = (),
) -> dict[str, dict[str, Any]]:
    """Construct cross sections and spin observables for loose MC dijets."""

    output: dict[str, dict[str, Any]] = {}
    resolved = set(helicity_resolved_observables)
    unknown = sorted(resolved - set(samples_by_object))
    if unknown:
        raise CampaignError(
            "Unknown helicity-resolved MC_POLDIJETS observables: "
            + ", ".join(unknown)
        )
    for observable, samples in samples_by_object.items():
        edges = list(next(iter(samples.values())).edges)
        sigma_uu = experimental.linear_combine_series(
            samples, {label: 0.25 for label in DENOMINATOR}
        )
        delta_sigma_ll = experimental.linear_combine_series(
            samples, DELTA_SIGMA_LL_COEFFICIENTS
        )
        output[f"SigmaUU_{observable}"] = {
            "edges": edges,
            "values": list(sigma_uu.values),
            "errors": [
                math.sqrt(max(0.0, variance))
                for variance in sigma_uu.variances
            ],
        }
        output[f"DeltaSigmaLL_{observable}"] = {
            "edges": edges,
            "values": list(delta_sigma_ll.values),
            "errors": [
                math.sqrt(max(0.0, variance))
                for variance in delta_sigma_ll.variances
            ],
        }
        values, errors = _ratio_arrays(
            samples, ALL_NUMERATOR, DENOMINATOR
        )
        output[f"ALL_{observable}"] = {
            "edges": edges, "values": values, "errors": errors,
        }
        if observable in resolved:
            for helicity in DENOMINATOR:
                series = samples[helicity]
                output[f"Sigma{helicity}_{observable}"] = {
                    "edges": edges,
                    "values": list(series.values),
                    "errors": [
                        math.sqrt(max(0.0, variance))
                        for variance in series.variances
                    ],
                }

        for label, numerator in (
            ("SingleSpinA", AL_A_NUMERATOR),
            ("SingleSpinB", AL_B_NUMERATOR),
        ):
            closure_values, closure_errors = _ratio_arrays(
                samples, numerator, DENOMINATOR
            )
            output[f"{label}_{observable}"] = {
                "edges": edges,
                "values": closure_values,
                "errors": closure_errors,
            }
        for first, second, label in (
            ("PP", "MM", "Parity_PP_MM"),
            ("PM", "MP", "Parity_PM_MP"),
        ):
            closure_values, closure_errors = [], []
            for index in range(len(edges) - 1):
                value, error = experimental.parity_residual(
                    samples[first].values[index],
                    samples[first].variances[index],
                    samples[second].values[index],
                    samples[second].variances[index],
                )
                closure_values.append(value)
                closure_errors.append(error)
            output[f"{label}_{observable}"] = {
                "edges": edges,
                "values": closure_values,
                "errors": closure_errors,
            }
    return output


def _integrated_series(series: experimental.BinSeries) -> tuple[float, float]:
    value = 0.0
    variance = 0.0
    for index, content in enumerate(series.values):
        width = series.edges[index + 1] - series.edges[index]
        value += content * width
        variance += series.variances[index] * width * width
    return value, variance


def _normalized_shape(series: experimental.BinSeries) -> dict[str, Any]:
    denominator, denominator_variance = _integrated_series(series)
    values: list[float | None] = []
    errors: list[float | None] = []
    for index, numerator in enumerate(series.values):
        width = series.edges[index + 1] - series.edges[index]
        value, error = experimental.ratio_with_covariance(
            numerator,
            series.variances[index],
            denominator,
            denominator_variance,
            series.variances[index] * width,
        )
        values.append(value)
        errors.append(error)
    return {"edges": list(series.edges), "values": values, "errors": errors}


def _angular_moment(
    samples: Mapping[str, experimental.BinSeries],
    numerator_coefficients: Mapping[str, float],
    trigonometric: str,
) -> tuple[float | None, float | None]:
    """Return an exact binned A2/B2 ratio with shared-sample covariance."""

    first = next(iter(samples.values()))
    numerator = 0.0
    denominator = 0.0
    numerator_variance = 0.0
    denominator_variance = 0.0
    covariance = 0.0
    for helicity, series in samples.items():
        if not experimental._same_edges(first.edges, series.edges):
            raise CampaignError("Angular moment inputs have different bin edges")
        for index, content in enumerate(series.values):
            low, high = series.edges[index:index + 2]
            centre = 0.5 * (low + high)
            width = high - low
            harmonic = (
                math.cos(2.0 * centre)
                if trigonometric == "cos"
                else math.sin(2.0 * centre)
            )
            numerator_weight = (
                2.0 * float(numerator_coefficients[helicity])
                * harmonic * width
            )
            denominator_weight = 0.25 * width
            variance = series.variances[index]
            numerator += numerator_weight * content
            denominator += denominator_weight * content
            numerator_variance += numerator_weight**2 * variance
            denominator_variance += denominator_weight**2 * variance
            covariance += numerator_weight * denominator_weight * variance
    return experimental.ratio_with_covariance(
        numerator, numerator_variance, denominator, denominator_variance,
        covariance,
    )


def _conditional_angular_moment(
    samples: Mapping[str, experimental.BinSeries],
    numerator_coefficients: Mapping[str, float],
    denominator_coefficients: Mapping[str, float],
    trigonometric: str,
    fraction_edges: Sequence[float],
    angle_bins: int,
) -> dict[str, Any]:
    """Project a flattened (fraction, angle) histogram with exact covariance."""

    if angle_bins <= 0 or len(fraction_edges) < 2:
        raise CampaignError("Conditional angular binning is empty")
    first = next(iter(samples.values()))
    fraction_bins = len(fraction_edges) - 1
    expected_bins = fraction_bins * angle_bins
    if len(first.values) != expected_bins:
        raise CampaignError(
            "Conditional angular histogram has "
            f"{len(first.values)} bins, expected {expected_bins}"
        )
    for series in samples.values():
        if not experimental._same_edges(first.edges, series.edges):
            raise CampaignError(
                "Conditional angular inputs have different flattened edges"
            )

    values: list[float | None] = []
    errors: list[float | None] = []
    for fraction_index in range(fraction_bins):
        numerator = 0.0
        denominator = 0.0
        numerator_variance = 0.0
        denominator_variance = 0.0
        covariance = 0.0
        for helicity, series in samples.items():
            numerator_coefficient = float(
                numerator_coefficients.get(helicity, 0.0)
            )
            denominator_coefficient = float(
                denominator_coefficients.get(helicity, 0.0)
            )
            for angle_index in range(angle_bins):
                flat_index = fraction_index * angle_bins + angle_index
                centre = -math.pi + (
                    float(angle_index) + 0.5
                ) * 2.0 * math.pi / float(angle_bins)
                harmonic = (
                    math.cos(2.0 * centre)
                    if trigonometric == "cos"
                    else math.sin(2.0 * centre)
                )
                numerator_weight = 2.0 * numerator_coefficient * harmonic
                denominator_weight = denominator_coefficient
                content = series.values[flat_index]
                variance = series.variances[flat_index]
                numerator += numerator_weight * content
                denominator += denominator_weight * content
                numerator_variance += numerator_weight**2 * variance
                denominator_variance += denominator_weight**2 * variance
                covariance += (
                    numerator_weight * denominator_weight * variance
                )
        value, error = experimental.ratio_with_covariance(
            numerator,
            numerator_variance,
            denominator,
            denominator_variance,
            covariance,
        )
        values.append(value)
        errors.append(error)
    return {
        "edges": [float(value) for value in fraction_edges],
        "values": values,
        "errors": errors,
    }


def _nested_ratio_prediction(
    numerator: experimental.BinSeries,
    denominator: experimental.BinSeries,
) -> dict[str, Any]:
    if not experimental._same_edges(numerator.edges, denominator.edges):
        raise CampaignError("Nested rate scans have different bin edges")
    values: list[float | None] = []
    errors: list[float | None] = []
    for index in range(len(numerator.values)):
        value, error = experimental.ratio_with_covariance(
            numerator.values[index], numerator.variances[index],
            denominator.values[index], denominator.variances[index],
            # The numerator is an event subset of the denominator.
            numerator.variances[index],
        )
        values.append(value)
        errors.append(error)
    return {"edges": list(numerator.edges), "values": values, "errors": errors}


def _one_minus_prediction(prediction: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "edges": list(prediction["edges"]),
        "values": [None if value is None else 1.0 - value
                   for value in prediction["values"]],
        "errors": list(prediction["errors"]),
    }


def _mc_poljetshapes_prediction(
    samples_by_object: Mapping[str, Mapping[str, experimental.BinSeries]],
    angular_observables: Sequence[str],
    conditional_angular_observables: Mapping[str, Mapping[str, Any]] | None = None,
    support_observables: Sequence[str] = (),
) -> dict[str, dict[str, Any]]:
    """Build helicity spectra, normalized shapes, moments, and rate scans."""

    output = _mc_poldijets_prediction(
        samples_by_object, tuple(samples_by_object)
    )
    angular = set(angular_observables)
    unknown = sorted(angular - set(samples_by_object))
    if unknown:
        raise CampaignError(
            "Unknown MC_POLJETSHAPES angular observables: "
            + ", ".join(unknown)
        )
    uu_coefficients = {label: 0.25 for label in DENOMINATOR}
    for observable in sorted(angular):
        samples = samples_by_object[observable]
        sigma_uu = experimental.linear_combine_series(samples, uu_coefficients)
        output[f"ShapeUU_{observable}"] = _normalized_shape(sigma_uu)
        for helicity in DENOMINATOR:
            output[f"Shape{helicity}_{observable}"] = _normalized_shape(
                samples[helicity]
            )
        for prefix, coefficients in (
            ("A2UU", uu_coefficients),
            ("A2LL", DELTA_SIGMA_LL_COEFFICIENTS),
        ):
            value, error = _angular_moment(samples, coefficients, "cos")
            output[f"{prefix}_{observable}"] = {
                "edges": [0.0, 1.0], "values": [value], "errors": [error],
            }
        for prefix, coefficients in (
            ("B2UU", uu_coefficients),
            ("B2LL", DELTA_SIGMA_LL_COEFFICIENTS),
        ):
            value, error = _angular_moment(samples, coefficients, "sin")
            output[f"{prefix}_{observable}"] = {
                "edges": [0.0, 1.0], "values": [value], "errors": [error],
            }

    required_rates = {
        "dijet_threshold_denominator", "ge3_threshold", "ge4_threshold"
    }
    if required_rates <= set(samples_by_object):
        denominator_samples = samples_by_object["dijet_threshold_denominator"]
        ge3_samples = samples_by_object["ge3_threshold"]
        ge4_samples = samples_by_object["ge4_threshold"]
        labels: dict[str, tuple[experimental.BinSeries,
                               experimental.BinSeries,
                               experimental.BinSeries]] = {}
        for helicity in DENOMINATOR:
            labels[helicity] = (
                denominator_samples[helicity], ge3_samples[helicity],
                ge4_samples[helicity],
            )
        labels["UU"] = (
            experimental.linear_combine_series(
                denominator_samples, uu_coefficients
            ),
            experimental.linear_combine_series(ge3_samples, uu_coefficients),
            experimental.linear_combine_series(ge4_samples, uu_coefficients),
        )
        for label, (denominator, ge3, ge4) in labels.items():
            r32 = _nested_ratio_prediction(ge3, denominator)
            r43 = _nested_ratio_prediction(ge4, ge3)
            output[f"R32_{label}"] = r32
            output[f"R43_{label}"] = r43
            output[f"ThirdJetVeto_{label}"] = _one_minus_prediction(r32)

    support = set(support_observables)
    for output_name in list(output):
        if any(output_name.endswith("_" + name) for name in support):
            del output[output_name]

    uu_coefficients = {label: 0.25 for label in DENOMINATOR}
    for output_name, definition in sorted(
        (conditional_angular_observables or {}).items()
    ):
        raw_observable = str(definition["raw_observable"])
        if raw_observable not in samples_by_object:
            raise CampaignError(
                f"Missing conditional angular input {raw_observable}"
            )
        fraction_edges = [
            float(value) for value in definition["fraction_edges"]
        ]
        angle_bins = int(definition["angle_bins"])
        samples = samples_by_object[raw_observable]
        for label, numerator_coefficients, denominator_coefficients in (
            ("UU", uu_coefficients, uu_coefficients),
            ("LL", DELTA_SIGMA_LL_COEFFICIENTS, uu_coefficients),
            *(
                (
                    helicity,
                    {name: 1.0 if name == helicity else 0.0
                     for name in DENOMINATOR},
                    {name: 1.0 if name == helicity else 0.0
                     for name in DENOMINATOR},
                )
                for helicity in DENOMINATOR
            ),
        ):
            output[f"C2{label}_{output_name}"] = _conditional_angular_moment(
                samples,
                numerator_coefficients,
                denominator_coefficients,
                "cos",
                fraction_edges,
                angle_bins,
            )
            output[f"S2{label}_{output_name}"] = _conditional_angular_moment(
                samples,
                numerator_coefficients,
                denominator_coefficients,
                "sin",
                fraction_edges,
                angle_bins,
            )
    return output


def _independent_difference(
    spin_on: Mapping[str, Any], spin_off: Mapping[str, Any]
) -> dict[str, Any]:
    if not experimental._same_edges(spin_on["edges"], spin_off["edges"]):
        raise CampaignError("Spin-on/off predictions have different bin edges")
    values: list[float | None] = []
    errors: list[float | None] = []
    for on_value, on_error, off_value, off_error in zip(
        spin_on["values"], spin_on["errors"],
        spin_off["values"], spin_off["errors"],
    ):
        if None in (on_value, on_error, off_value, off_error):
            values.append(None)
            errors.append(None)
            continue
        values.append(float(on_value) - float(off_value))
        errors.append(math.hypot(float(on_error), float(off_error)))
    return {"edges": list(spin_on["edges"]), "values": values, "errors": errors}


def _independent_ratio(
    spin_on: Mapping[str, Any], spin_off: Mapping[str, Any],
    minimum_denominator_sigma: float = 0.0,
) -> dict[str, Any]:
    """Return spin-on/spin-off with both independent MC errors propagated."""

    if not experimental._same_edges(spin_on["edges"], spin_off["edges"]):
        raise CampaignError("Spin-on/off predictions have different bin edges")
    values: list[float | None] = []
    errors: list[float | None] = []
    for on_value, on_error, off_value, off_error in zip(
        spin_on["values"], spin_on["errors"],
        spin_off["values"], spin_off["errors"],
    ):
        if None in (on_value, on_error, off_value, off_error):
            values.append(None)
            errors.append(None)
            continue
        numerator = float(on_value)
        numerator_error = float(on_error)
        denominator = float(off_value)
        denominator_error = float(off_error)
        if (
            denominator == 0.0
            or (minimum_denominator_sigma > 0.0 and (
                numerator < 0.0 or denominator <= 0.0
                or denominator_error <= 0.0
                or denominator < minimum_denominator_sigma * denominator_error
            ))
            or not all(
                math.isfinite(value)
                for value in (
                    numerator, numerator_error,
                    denominator, denominator_error,
                )
            )
        ):
            values.append(None)
            errors.append(None)
            continue
        ratio = numerator / denominator
        variance = (
            (numerator_error / denominator) ** 2
            + (numerator * denominator_error / denominator**2) ** 2
        )
        if not math.isfinite(ratio) or not math.isfinite(variance):
            values.append(None)
            errors.append(None)
            continue
        values.append(ratio)
        errors.append(math.sqrt(max(0.0, variance)))
    return {"edges": list(spin_on["edges"]), "values": values, "errors": errors}


def _mc_poljetshapes_plot_ratios(
    measurement: Mapping[str, Any], summary: Mapping[str, Any]
) -> dict[str, dict[str, Any]]:
    """Build presentation-only positive-observable ratios from summary JSON.

    The campaign summary already contains the postprocessed values and MC
    errors for both independent families.  Reusing it lets a completed
    immutable campaign acquire ratio plots with a plot-only rerun: no shard or
    postprocessing product is reinterpreted or overwritten.
    """

    variations = summary.get("variations", {})
    if not isinstance(variations, Mapping):
        return {}
    ratios: dict[str, dict[str, Any]] = {}
    for channel_id in measurement["channels"]:
        nominal_key = _variation_id(
            (
                "nominal", channel_id, 0, 0, 1.0,
                measurement["families"]["nominal"]["mpi"],
            )
        )
        control_key = _variation_id(
            (
                _comparison_pair(measurement)[1], channel_id, 0, 0, 1.0,
                measurement["families"][_comparison_pair(measurement)[1]]["mpi"],
            )
        )
        nominal = variations.get(nominal_key)
        control = variations.get(control_key)
        if not isinstance(nominal, Mapping) or not isinstance(control, Mapping):
            continue
        for observable in sorted(set(nominal).intersection(control)):
            if not observable.startswith(MC_POLJETSHAPES_RATIO_PREFIXES):
                continue
            if observable in ratios:
                raise CampaignError(
                    f"Duplicate MC_POLJETSHAPES ratio observable {observable}"
                )
            ratios[observable] = _independent_ratio(
                nominal[observable], control[observable],
                float(measurement.get("comparison_support", {}).get(
                    "minimum_ratio_denominator_sigma", 0.0
                )),
            )
    return ratios


def _effective_entries(series: experimental.BinSeries) -> float:
    value, variance = _integrated_series(series)
    if variance <= 0.0 or not math.isfinite(value) or not math.isfinite(variance):
        return 0.0
    return value * value / variance


def _mask_spin_asymmetry_support(
    measurement: Mapping[str, Any],
    predictions: Mapping[tuple[Any, ...], dict[str, Any]],
    raw_samples: Mapping[tuple[Any, ...], Any],
) -> None:
    """Conservative MC-support masks; never censor an observed discrepancy."""
    minimum = float(measurement.get("comparison_support", {}).get(
        "minimum_effective_entries_per_helicity", 0.0
    ))
    if minimum <= 0.0:
        return
    for key, outputs in predictions.items():
        for name, prediction in outputs.items():
            match = re.match(r"^(ALL|[ABCS]2(?:UU|LL))_(.+)$", name)
            if not match:
                continue
            kind, observable = match.groups()
            conditional = measurement["channels"][key[1]].get(
                "conditional_angular_observables", {}
            ).get(observable) if kind.startswith(("C2", "S2")) else None
            raw_observable = conditional["raw_observable"] if conditional else observable
            samples = raw_samples.get(key, {}).get(raw_observable, {})
            mask = []
            for index in range(len(prediction["values"])):
                valid = set(DENOMINATOR) <= set(samples)
                for sample in samples.values():
                    selection = slice(None)
                    if kind == "ALL":
                        selection = slice(index, index + 1)
                    elif conditional:
                        angle_bins = int(conditional["angle_bins"])
                        selection = slice(index * angle_bins, (index + 1) * angle_bins)
                    value = sum(sample.values[selection])
                    variance = sum(sample.variances[selection])
                    neff = value * value / variance if variance > 0 else 0.0
                    valid = valid and value > 0 and math.isfinite(neff) and neff >= minimum
                valid = valid and None not in (
                    prediction["values"][index], prediction["errors"][index]
                )
                if valid:
                    valid = (math.isfinite(prediction["values"][index]) and
                             math.isfinite(prediction["errors"][index]) and
                             prediction["errors"][index] > 0)
                mask.append(bool(valid))
                if not valid:
                    prediction["values"][index] = None
                    prediction["errors"][index] = None
            prediction["support_mask"] = mask
            prediction["support_policy"] = f"at least {minimum:g} effective entries in each helicity"


def _mc_poljetshapes_assessment(
    measurement: Mapping[str, Any],
    predictions: Mapping[tuple[Any, ...], Mapping[str, Any]],
    raw_samples: Mapping[
        tuple[Any, ...], Mapping[str, Mapping[str, experimental.BinSeries]]
    ],
    manifest: Mapping[str, Any],
) -> tuple[dict[str, Any], dict[str, dict[str, Any]], list[dict[str, Any]]]:
    """Compare independent shower samples and project bounded production tiers."""

    if not {"nominal", _comparison_pair(measurement)[1]} <= set(
        manifest["configuration"]["families"]
    ):
        return {}, {}, []
    comparison: dict[str, dict[str, Any]] = {}
    ranking: list[dict[str, Any]] = []
    moment_differences: dict[str, Any] = {}
    baseline_effective: dict[str, dict[str, dict[str, float]]] = {}
    rare_effective: dict[str, dict[str, dict[str, float]]] = {}
    projection = measurement.get("statistics_projection", {})
    baseline_names = tuple(
        projection.get(
            "baseline_observables",
            ("dpsi12_j1_loose", "dpsi12_j2_loose"),
        )
    )
    rare_names = tuple(
        projection.get(
            "rare_observables",
            ("interjet_dpsi11_kt05", "bz_angle"),
        )
    )
    if "comparison_pair" in measurement:
        rare_names = tuple(dict.fromkeys((*rare_names, "interjet_dpsi11_kt05",
                                          "interjet_dpsi11_kt10", "bz_angle")))

    for channel_id, channel in measurement["channels"].items():
        nominal_key = ("nominal", channel_id, 0, 0, 1.0,
                       measurement["families"]["nominal"]["mpi"])
        control_key = (_comparison_pair(measurement)[1], channel_id, 0, 0, 1.0,
                       measurement["families"][_comparison_pair(measurement)[1]]["mpi"])
        if nominal_key not in predictions or control_key not in predictions:
            continue
        nominal = predictions[nominal_key]
        control = predictions[control_key]
        for observable in sorted(set(nominal).intersection(control)):
            if observable.startswith(("SingleSpin", "Parity_")):
                continue
            difference = _independent_difference(
                nominal[observable], control[observable]
            )
            comparison[observable] = difference
            pulls = [
                float(value) / float(error)
                for value, error in zip(
                    difference["values"], difference["errors"]
                )
                if value is not None and error is not None and error > 0.0
            ]
            if pulls:
                ranking.append(
                    {
                        "observable": observable,
                        "bins": len(pulls),
                        "chi2_independent_samples": sum(
                            pull * pull for pull in pulls
                        ),
                        "quadrature_sensitivity": math.sqrt(
                            sum(pull * pull for pull in pulls)
                        ),
                        "maximum_absolute_pull": max(abs(pull) for pull in pulls),
                    }
                )
            if observable.startswith(("A2UU_", "A2LL_", "B2UU_", "B2LL_")):
                moment_differences[observable] = difference

        for family_id, key in (("nominal", nominal_key),
                               (_comparison_pair(measurement)[1], control_key)):
            baseline_effective.setdefault(family_id, {})
            rare_effective.setdefault(family_id, {})
            for helicity, series_map in (
                (helicity, {
                    name: raw_samples[key][name][helicity]
                    for name in (*baseline_names, *rare_names)
                    if name in raw_samples[key]
                })
                for helicity in DENOMINATOR
            ):
                baseline_effective[family_id][helicity] = {
                    name: _effective_entries(series_map[name])
                    for name in baseline_names if name in series_map
                }
                rare_effective[family_id][helicity] = {
                    name: _effective_entries(series_map[name])
                    for name in rare_names if name in series_map
                }

    ranking.sort(
        key=lambda entry: (
            -float(entry["quadrature_sensitivity"]),
            str(entry["observable"]),
        )
    )
    pilot_events = int(manifest["configuration"]["lo_events"])
    candidates = tuple(
        int(value)
        for value in projection.get(
            "candidate_events_per_helicity_family",
            (100_000_000, 250_000_000, 500_000_000),
        )
    )
    minimum_entries = float(
        projection.get("minimum_effective_baseline_entries", 250_000)
    )
    maximum_moment_error_target = float(
        projection.get("maximum_independent_on_off_A2_error", 0.002)
    )
    tier_assessments: list[dict[str, Any]] = []
    selected_events: int | None = None
    baseline_values = [
        value
        for families in baseline_effective.values()
        for helicities in families.values()
        for value in helicities.values()
    ]
    relevant_moments = [
        prediction
        for name, prediction in moment_differences.items()
        if name in {
            *(f"A2UU_{observable}" for observable in baseline_names),
            *(f"A2LL_{observable}" for observable in baseline_names),
        }
    ]
    for candidate in candidates:
        scale = float(candidate)/pilot_events
        minimum_effective = (
            min(baseline_values) * scale if baseline_values else 0.0
        )
        projected_moment_errors = [
            float(prediction["errors"][0]) / math.sqrt(scale)
            for prediction in relevant_moments
            if prediction["errors"][0] is not None
        ]
        maximum_moment_error = (
            max(projected_moment_errors)
            if len(projected_moment_errors) == len(relevant_moments)
            and projected_moment_errors
            else None
        )
        passes = (
            minimum_effective >= minimum_entries
            and maximum_moment_error is not None
            and maximum_moment_error <= maximum_moment_error_target
        )
        tier_assessments.append(
            {
                "events_per_helicity_family": candidate,
                "projected_minimum_effective_baseline_entries": minimum_effective,
                "projected_maximum_independent_on_off_A2_error": maximum_moment_error,
                "passes": passes,
            }
        )
        if passes and selected_events is None:
            selected_events = candidate

    bounded_recommendation = selected_events or candidates[-1]
    tier_millions = bounded_recommendation // 1_000_000
    tier_label = (
        f"{bounded_recommendation // 1_000_000_000}b"
        if bounded_recommendation % 1_000_000_000 == 0
        else f"{tier_millions}m"
    )
    events_per_shard = int(projection.get("events_per_shard", 500_000))
    shards = bounded_recommendation // events_per_shard
    production_seed_base = int(
        projection.get("production_seed_base", 8307000)
    )
    production_tag_date = str(
        projection.get("production_tag_date", "20260829")
    )
    production_tag_stem = str(
        projection.get("production_tag_stem", "mc_poljetshapes_spin")
    )
    command = (
        f"python3 {measurement.get('runner', 'scripts/run_mc_poljetshapes_campaign.py')} full "
        f"--tag {production_tag_stem}_{tier_label}_"
        f"{production_tag_date}_v1 "
        f"--families nominal,{_comparison_pair(measurement)[1]} "
        f"--lo-events {bounded_recommendation} --shards {shards} --jobs 100 "
        f"--seed-base {production_seed_base} "
        "--plot-comparisons --include-diagnostics"
    )
    assessment = {
        "pilot_events_per_helicity_family": pilot_events,
        "baseline_effective_entries": baseline_effective,
        "rare_effective_entries": rare_effective,
        "moment_differences": moment_differences,
        "tiers": tier_assessments,
        "selected_events_per_helicity_family": selected_events,
        "bounded_recommendation_events_per_helicity_family": bounded_recommendation,
        "selection_status": (
            "criteria_satisfied" if selected_events is not None
            else "no_bounded_tier_satisfies_all_criteria"
        ),
        "production_command": command,
        "criteria": {
            "minimum_effective_baseline_entries_per_helicity_family": (
                minimum_entries
            ),
            "maximum_independent_on_off_A2UU_or_A2LL_error": (
                maximum_moment_error_target
            ),
            "candidate_tiers_events": list(candidates),
        },
        "note": (
            "Rare channels are reported separately and do not force an event "
            "tier beyond the largest configured bounded candidate."
        ),
    }
    return assessment, comparison, ranking


def _exclude_shard_series(
    combined: experimental.BinSeries,
    excluded: experimental.BinSeries,
    total_events: int,
    excluded_events: int,
) -> experimental.BinSeries:
    """Return the event-weighted combination with one shard removed."""

    if not experimental._same_edges(combined.edges, excluded.edges):
        raise CampaignError("Jackknife shard has different bin edges")
    retained_events = total_events - excluded_events
    if retained_events <= 0:
        raise CampaignError("Cannot jackknife the only shard")
    values = [
        (
            float(total_events) * combined.values[index]
            - float(excluded_events) * excluded.values[index]
        ) / float(retained_events)
        for index in range(len(combined.values))
    ]
    variances = [
        max(
            0.0,
            (
                float(total_events)**2 * combined.variances[index]
                - float(excluded_events)**2 * excluded.variances[index]
            ) / float(retained_events)**2,
        )
        for index in range(len(combined.values))
    ]
    return experimental.BinSeries(
        list(combined.edges), values, variances
    )


def _jackknife_matrix(
    replicates: Sequence[Sequence[float | None]],
) -> dict[str, Any]:
    """Return a delete-one-block covariance and correlation matrix."""

    try:
        import numpy as np
    except ImportError as exc:
        raise CampaignError(
            "NumPy is required for shard-block covariance"
        ) from exc
    if len(replicates) < 3:
        raise CampaignError("At least three shard blocks are required")
    array = np.asarray(
        [
            [np.nan if value is None else float(value) for value in row]
            for row in replicates
        ],
        dtype=float,
    )
    active = np.all(np.isfinite(array), axis=0)
    size = array.shape[1]
    covariance_full: list[list[float | None]] = [
        [None for _ in range(size)] for _ in range(size)
    ]
    correlation_full: list[list[float | None]] = [
        [None for _ in range(size)] for _ in range(size)
    ]
    errors: list[float | None] = [None for _ in range(size)]
    active_indices = np.flatnonzero(active)
    if active_indices.size:
        selected = array[:, active]
        centred = selected - np.mean(selected, axis=0)
        blocks = float(array.shape[0])
        covariance = (blocks - 1.0) / blocks * centred.T.dot(centred)
        diagonal = np.maximum(np.diag(covariance), 0.0)
        denominator = np.sqrt(np.outer(diagonal, diagonal))
        correlation = np.divide(
            covariance,
            denominator,
            out=np.zeros_like(covariance),
            where=denominator > 0.0,
        )
        for local_i, global_i in enumerate(active_indices):
            errors[int(global_i)] = math.sqrt(
                max(0.0, float(diagonal[local_i]))
            )
            for local_j, global_j in enumerate(active_indices):
                covariance_full[int(global_i)][int(global_j)] = float(
                    covariance[local_i, local_j]
                )
                correlation_full[int(global_i)][int(global_j)] = float(
                    correlation[local_i, local_j]
                )
    return {
        "active_bins": [bool(value) for value in active.tolist()],
        "errors": errors,
        "covariance": covariance_full,
        "correlation": correlation_full,
    }


def _mc_poljetshapes_shard_covariance(
    measurement: Mapping[str, Any],
    groups: Mapping[
        tuple[Any, ...], dict[str, list[Mapping[str, Any]]]
    ],
    campaign_dir: Path,
    raw_samples: Mapping[
        tuple[Any, ...], Mapping[str, Mapping[str, experimental.BinSeries]]
    ],
    predictions: Mapping[tuple[Any, ...], Mapping[str, Any]],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Estimate selected full-result covariances from delete-one shard blocks."""

    configuration = measurement.get("shard_block_covariance")
    if not isinstance(configuration, Mapping):
        return {}, []
    raw_names = tuple(str(value) for value in configuration["raw_observables"])
    output_names = tuple(
        str(value) for value in configuration["prediction_observables"]
    )
    payload: dict[str, Any] = {
        "method": "delete-one-common-shard jackknife",
        "block_definition": (
            "one equal-ordinal shard from each independent physical-helicity "
            "sample in a family"
        ),
        "raw_observables": list(raw_names),
        "prediction_observables": list(output_names),
        "families": {},
        "spin_on_minus_off": {},
    }
    csv_rows: list[dict[str, Any]] = []
    channel_id = next(iter(measurement["channels"]))
    channel = measurement["channels"][channel_id]
    central_keys: dict[str, tuple[Any, ...]] = {}

    for family_id in ("nominal", _comparison_pair(measurement)[1]):
        if family_id not in measurement["families"]:
            continue
        key = (
            family_id,
            channel_id,
            0,
            0,
            1.0,
            measurement["families"][family_id]["mpi"],
        )
        if key not in groups or key not in raw_samples or key not in predictions:
            continue
        central_keys[family_id] = key
        helicity_jobs = groups[key]
        jobs_by_helicity = {
            helicity: {int(job["shard"]): job for job in jobs}
            for helicity, jobs in helicity_jobs.items()
        }
        shard_sets = [set(jobs) for jobs in jobs_by_helicity.values()]
        if not shard_sets or any(shards != shard_sets[0] for shards in shard_sets[1:]):
            raise CampaignError(
                f"{family_id} physical helicities have different shard ordinals"
            )
        block_ids = sorted(shard_sets[0])
        if len(block_ids) < 3:
            payload["families"][family_id] = {
                "status": "insufficient_common_blocks",
                "blocks": len(block_ids),
            }
            continue

        event_vectors = {
            tuple(
                int(jobs_by_helicity[helicity][block]["events"])
                for helicity in DENOMINATOR
            )
            for block in block_ids
        }
        if (
            len(event_vectors) != 1
            or len(set(next(iter(event_vectors)))) != 1
        ):
            payload["families"][family_id] = {
                "status": "unequal_block_sizes_not_supported",
                "blocks": len(block_ids),
                "event_vectors": [list(values) for values in sorted(event_vectors)],
            }
            continue

        total_events = {
            helicity: sum(int(job["events"]) for job in jobs.values())
            for helicity, jobs in jobs_by_helicity.items()
        }
        shard_cache: dict[
            str, dict[str, dict[int, experimental.BinSeries]]
        ] = {
            raw_name: {helicity: {} for helicity in DENOMINATOR}
            for raw_name in raw_names
        }
        requested: dict[str, tuple[str, str]] = {}
        for raw_name in raw_names:
            if raw_name not in raw_samples[key]:
                raise CampaignError(
                    f"Shard covariance input {raw_name} is not configured"
                )
            object_spec = channel["raw_objects"][raw_name]
            source_analysis, object_name = _raw_object_source(
                measurement, object_spec
            )
            requested[raw_name] = (source_analysis, object_name)
        for helicity in DENOMINATOR:
            for block_id in block_ids:
                loaded = _load_series_many(
                    [jobs_by_helicity[helicity][block_id]],
                    campaign_dir,
                    requested,
                )
                for raw_name, item in loaded.items():
                    shard_cache[raw_name][helicity][block_id] = item

        replicates: dict[str, list[list[float | None]]] = {
            name: [] for name in output_names
        }
        for block_id in block_ids:
            leave_one_out: dict[
                str, dict[str, experimental.BinSeries]
            ] = {}
            for raw_name in raw_names:
                leave_one_out[raw_name] = {}
                for helicity in DENOMINATOR:
                    job = jobs_by_helicity[helicity][block_id]
                    leave_one_out[raw_name][helicity] = _exclude_shard_series(
                        raw_samples[key][raw_name][helicity],
                        shard_cache[raw_name][helicity][block_id],
                        total_events[helicity],
                        int(job["events"]),
                    )
            replicate_prediction = _mc_poljetshapes_prediction(
                leave_one_out,
                (),
                channel.get("conditional_angular_observables", {}),
                channel.get("support_observables", []),
            )
            for output_name in output_names:
                if output_name not in replicate_prediction:
                    raise CampaignError(
                        f"Shard covariance output {output_name} was not produced"
                    )
                replicates[output_name].append(
                    list(replicate_prediction[output_name]["values"])
                )

        family_payload: dict[str, Any] = {
            "status": "complete",
            "blocks": len(block_ids),
            "block_ids": block_ids,
            "equal_events_per_block": True,
            "observables": {},
        }
        for output_name in output_names:
            matrix = _jackknife_matrix(replicates[output_name])
            prediction = predictions[key][output_name]
            result = {
                "edges": list(prediction["edges"]),
                "values": list(prediction["values"]),
                **matrix,
            }
            family_payload["observables"][output_name] = result
            for index, row in enumerate(matrix["covariance"]):
                for column, covariance in enumerate(row):
                    csv_rows.append(
                        {
                            "sample": family_id,
                            "observable": output_name,
                            "bin_i": index,
                            "bin_j": column,
                            "covariance": covariance,
                            "correlation": matrix["correlation"][index][column],
                        }
                    )
        payload["families"][family_id] = family_payload

    if {"nominal", _comparison_pair(measurement)[1]} <= set(central_keys):
        nominal_family = payload["families"].get("nominal", {})
        control_family = payload["families"].get(_comparison_pair(measurement)[1], {})
        if (
            nominal_family.get("status") == "complete"
            and control_family.get("status") == "complete"
        ):
            try:
                import numpy as np
            except ImportError as exc:
                raise CampaignError(
                    "NumPy is required for shard-block covariance"
                ) from exc
            nominal_key = central_keys["nominal"]
            control_key = central_keys[_comparison_pair(measurement)[1]]
            for output_name in output_names:
                nominal_matrix = nominal_family["observables"][output_name]
                control_matrix = control_family["observables"][output_name]
                size = len(nominal_matrix["values"])
                covariance: list[list[float | None]] = [
                    [None for _ in range(size)] for _ in range(size)
                ]
                correlation: list[list[float | None]] = [
                    [None for _ in range(size)] for _ in range(size)
                ]
                active = [
                    bool(nominal_matrix["active_bins"][index])
                    and bool(control_matrix["active_bins"][index])
                    and predictions[nominal_key][output_name]["values"][index]
                    is not None
                    and predictions[control_key][output_name]["values"][index]
                    is not None
                    and math.isfinite(float(
                        predictions[nominal_key][output_name]["values"][index]
                    ))
                    and math.isfinite(float(
                        predictions[control_key][output_name]["values"][index]
                    ))
                    for index in range(size)
                ]
                active_indices = [
                    index for index, keep in enumerate(active) if keep
                ]
                if active_indices:
                    dense = np.asarray(
                        [
                            [
                                float(nominal_matrix["covariance"][i][j])
                                + float(control_matrix["covariance"][i][j])
                                for j in active_indices
                            ]
                            for i in active_indices
                        ],
                        dtype=float,
                    )
                    diagonal = np.maximum(np.diag(dense), 0.0)
                    denominator = np.sqrt(np.outer(diagonal, diagonal))
                    dense_correlation = np.divide(
                        dense,
                        denominator,
                        out=np.zeros_like(dense),
                        where=denominator > 0.0,
                    )
                    for local_i, global_i in enumerate(active_indices):
                        for local_j, global_j in enumerate(active_indices):
                            covariance[global_i][global_j] = float(
                                dense[local_i, local_j]
                            )
                            correlation[global_i][global_j] = float(
                                dense_correlation[local_i, local_j]
                            )
                    difference = np.asarray(
                        [
                            float(predictions[nominal_key][output_name]["values"][i])
                            - float(predictions[control_key][output_name]["values"][i])
                            for i in active_indices
                        ],
                        dtype=float,
                    )
                    inverse = np.linalg.pinv(
                        dense, hermitian=True, rcond=1.0e-12
                    )
                    chi2 = float(difference.dot(inverse).dot(difference))
                    rank = int(np.linalg.matrix_rank(dense))
                else:
                    chi2 = None
                    rank = 0
                difference_values = [
                    (
                        None
                        if predictions[nominal_key][output_name]["values"][index]
                        is None
                        or predictions[control_key][output_name]["values"][index]
                        is None
                        else float(
                            predictions[nominal_key][output_name]["values"][index]
                        )
                        - float(
                            predictions[control_key][output_name]["values"][index]
                        )
                    )
                    for index in range(size)
                ]
                payload["spin_on_minus_off"][output_name] = {
                    "edges": list(
                        predictions[nominal_key][output_name]["edges"]
                    ),
                    "values": difference_values,
                    "active_bins": active,
                    "covariance": covariance,
                    "correlation": correlation,
                    "chi2": chi2,
                    "rank": rank,
                    "interpretation": (
                        "independent-family covariance sum; chi2 uses the "
                        "Moore-Penrose inverse"
                    ),
                }
                for index, row in enumerate(covariance):
                    for column, value in enumerate(row):
                        csv_rows.append(
                            {
                                "sample": "spin_on_minus_off",
                                "observable": output_name,
                                "bin_i": index,
                                "bin_j": column,
                                "covariance": value,
                                "correlation": correlation[index][column],
                            }
                        )
    return payload, csv_rows


def _apply_star_display_binning(
    predictions: dict[str, dict[str, Any]],
    snapshot: Mapping[str, Any],
) -> None:
    """Separate event-selection bins from unfolded plotting coordinates.

    STAR quotes the asymmetry in detector-level event intervals and reports
    each result at a corrected parton-level representative coordinate.  The
    physical yields are therefore accumulated in the published event bins,
    while the final estimates are displayed in a contiguous Voronoi binning
    around those published coordinates.
    """

    for identifier, dataset in snapshot["datasets"].items():
        prediction = predictions.get(str(identifier))
        if prediction is None:
            continue
        analysis_edges = [float(value) for value in dataset["bin_edges"]]
        if not experimental._same_edges(
            prediction["edges"], analysis_edges
        ):
            raise CampaignError(
                f"STAR analysis bins disagree with the pinned reference for "
                f"{identifier}"
            )
        display_edges = [
            float(value) for value in dataset["display_bin_edges"]
        ]
        if len(display_edges) != len(prediction["values"]) + 1:
            raise CampaignError(
                f"STAR display-bin count changed for {identifier}"
            )
        prediction["analysis_edges"] = analysis_edges
        prediction["edges"] = display_edges


def _variation_id(key: tuple[Any, ...]) -> str:
    family, channel, polarized, unpolarized, scale, mpi = key
    return (f"{family}-{channel}-p{polarized:03d}-u{unpolarized:03d}-"
            f"mu{_scale_token(float(scale))}-mpi{mpi}")


def _compass_reference_entry(
    snapshot: Mapping[str, Any], observable: str
) -> tuple[Mapping[str, Any], Mapping[str, Any] | None] | None:
    dataset = snapshot.get("datasets", {}).get(observable)
    if dataset is not None:
        return dataset, None
    projection = snapshot.get("readable_projections", {}).get(observable)
    if projection is not None:
        return projection, None
    for candidate in snapshot.get("datasets", {}).values():
        for slice_spec in candidate.get("slices", []):
            if Path(str(slice_spec["rivet_path"])).name == observable:
                return candidate, slice_spec
    return None


def _reference_path(measurement: Mapping[str, Any], observable: str,
                    snapshot: Mapping[str, Any]) -> str | None:
    if measurement["postprocessor"] == "star_weak_bosons":
        if observable.startswith("Wplus_"):
            channel = snapshot["channels"]["Wplus"]
            return channel["rivet_path"] if observable.endswith("_AL") else channel["all"]["rivet_path"]
        if observable.startswith("Wminus_"):
            channel = snapshot["channels"]["Wminus"]
            return channel["rivet_path"] if observable.endswith("_AL") else channel["all"]["rivet_path"]
        if observable == "Zgamma_AL":
            return snapshot["channels"]["Zgamma"]["rivet_path"]
    if measurement["postprocessor"] == "phenix_prompt_photon" and observable in snapshot["channels"]:
        return snapshot["channels"][observable]["rivet_path"]
    if measurement["postprocessor"] == "star_jet_all":
        dataset = snapshot["datasets"].get(observable)
        return str(dataset["rivet_path"]) if dataset else None
    if measurement["postprocessor"] in {"mc_poldijets", "mc_poljetshapes"}:
        for prefix in (
            "ALL_", "DeltaSigmaLL_", "SigmaUU_",
            "SigmaPP_", "SigmaPM_", "SigmaMP_", "SigmaMM_",
        ):
            if observable.startswith(prefix):
                raw_observable = observable[len(prefix):]
                if raw_observable in snapshot.get("observables", {}):
                    return f"/{measurement['analysis']['name']}/{observable}"
        if measurement["postprocessor"] == "mc_poljetshapes":
            for prefix in (
                "ShapeUU_", "ShapePP_", "ShapePM_", "ShapeMP_", "ShapeMM_",
                "A2UU_", "A2LL_", "B2UU_", "B2LL_",
            ):
                if observable.startswith(prefix):
                    raw_observable = observable[len(prefix):]
                    definition = snapshot.get("observables", {}).get(
                        raw_observable, {}
                    )
                    if definition.get("angular"):
                        return f"/{measurement['analysis']['name']}/{observable}"
            conditional_match = re.match(
                r"^(C2|S2)(UU|LL|PP|PM|MP|MM)_(.+)$", observable
            )
            if conditional_match:
                conditional_name = conditional_match.group(3)
                if any(
                    conditional_name in channel.get(
                        "conditional_angular_observables", {}
                    )
                    for channel in measurement["channels"].values()
                ):
                    return f"/{measurement['analysis']['name']}/{observable}"
            if re.match(r"^(R32|R43|ThirdJetVeto)_(UU|PP|PM|MP|MM)$",
                        observable):
                return f"/{measurement['analysis']['name']}/{observable}"
        return None
    if measurement["postprocessor"] == "hermes_sidis":
        dataset = snapshot.get("datasets", {}).get(observable)
        if dataset and dataset.get("observable") == "A_parallel":
            return str(dataset["rivet_path"])
    if measurement["postprocessor"] in {
        "compass_sidis_a1", "compass_sidis_multiplicity",
        "hermes_sidis_multiplicity", "compass_sidis_charge_ratio",
        "compass_sidis_pt_slope", "sidis_azimuthal_moment",
    }:
        entry = _compass_reference_entry(snapshot, observable)
        if entry is None:
            return None
        dataset, slice_spec = entry
        if slice_spec is not None:
            return str(slice_spec["rivet_path"])
        return str(dataset.get("rivet_path") or dataset["flat_rivet_path"])
    return None


def _replica_sigma(predictions: Mapping[tuple[Any, ...], Mapping[str, Any]],
                   keys: Sequence[tuple[Any, ...]], observable: str,
                   central: Sequence[float | None]) -> list[float | None] | None:
    available = [predictions[key][observable]["values"] for key in keys
                 if key in predictions and observable in predictions[key]]
    if len(available) < 2:
        return None
    result: list[float | None] = []
    for index, centre in enumerate(central):
        values = [row[index] for row in available if row[index] is not None]
        if centre is None or len(values) < 2:
            result.append(None)
        else:
            mean = sum(float(value) for value in values) / len(values)
            result.append(
                math.sqrt(
                    sum((float(value) - mean) ** 2 for value in values)
                    / (len(values) - 1)
                )
            )
    return result


def aggregate_uncertainties(predictions: Mapping[tuple[Any, ...], Mapping[str, Any]],
                            measurement: Mapping[str, Any]) -> dict[str, Any]:
    bands: dict[str, Any] = {}
    snapshot = _measurement_snapshot(measurement)
    for channel in measurement["channels"]:
        nominal_mpi = measurement["families"]["nominal"]["mpi"]
        central_key = ("nominal", channel, 0, 0, 1.0, nominal_mpi)
        if central_key not in predictions:
            continue
        for observable, central_prediction in predictions[central_key].items():
            if _reference_path(measurement, observable, snapshot) is None:
                continue
            central = central_prediction["values"]
            polarized_keys = [("nominal", channel, member, 0, 1.0, nominal_mpi)
                              for member in range(1,101)]
            unpolarized_keys = [("nominal", channel, 0, member, 1.0, nominal_mpi)
                                for member in range(1,101)]
            pol = _replica_sigma(predictions, polarized_keys, observable, central)
            unpol = _replica_sigma(predictions, unpolarized_keys, observable, central)
            pdf: list[float | None] | None = None
            if pol is not None or unpol is not None:
                pdf = []
                for index in range(len(central)):
                    pieces = [array[index] for array in (pol,unpol)
                              if array is not None and array[index] is not None]
                    pdf.append(math.sqrt(sum(float(value)**2 for value in pieces))
                               if pieces else None)
            scale_rows = [predictions[key][observable]["values"]
                          for key in (("nominal",channel,0,0,0.5,nominal_mpi),
                                      central_key,
                                      ("nominal",channel,0,0,2.0,nominal_mpi))
                          if key in predictions and observable in predictions[key]]
            scale_low, scale_high = None, None
            if len(scale_rows) >= 2:
                scale_low = []; scale_high = []
                for index, centre in enumerate(central):
                    values = [row[index] for row in scale_rows if row[index] is not None]
                    if centre is None or not values:
                        scale_low.append(None); scale_high.append(None)
                    else:
                        scale_low.append(min(values)-float(centre))
                        scale_high.append(max(values)-float(centre))
            mpi_family = next(
                (family_id for family_id, family in measurement["families"].items()
                 if family_id != "nominal" and family.get("mpi") != nominal_mpi
                 and len(family.get("helicities", [])) == 4),
                None,
            )
            mpi_key = (
                (mpi_family, channel, 0, 0, 1.0,
                 measurement["families"][mpi_family]["mpi"])
                if mpi_family else None
            )
            mpi = None
            if mpi_key is not None and mpi_key in predictions and observable in predictions[mpi_key]:
                mpi = [None if centre is None or alternative is None else
                       float(alternative)-float(centre)
                       for centre, alternative in zip(central,
                           predictions[mpi_key][observable]["values"])]
            bands[observable] = {
                "pdf_68": pdf, "polarized_pdf_68": pol,
                "unpolarized_pdf_68": unpol,
                "scale_down": scale_low, "scale_up": scale_high,
                "mpi_shift": mpi,
                "monte_carlo": central_prediction["errors"],
            }
    return bands


def _estimate_with_bands(yoda: Any, prediction: Mapping[str, Any], path: str,
                         annotations: Mapping[str, Any], bands: Mapping[str, Any] | None = None) -> Any:
    result = experimental._estimate_from_values(
        yoda, prediction["edges"], path, prediction["values"], prediction["errors"],
        annotations)
    if not bands:
        return result
    for index, value in enumerate(prediction["values"], start=1):
        if value is None:
            continue
        bin_object = result.bin(index)
        pdf = bands.get("pdf_68")
        if pdf is not None and pdf[index-1] is not None:
            error = float(pdf[index-1])
            bin_object.setErr(-error, error, "pdf68")
        lower, upper = bands.get("scale_down"), bands.get("scale_up")
        if lower is not None and upper is not None and lower[index-1] is not None:
            bin_object.setErr(float(lower[index-1]), float(upper[index-1]), "hard_scale")
        mpi = bands.get("mpi_shift")
        if mpi is not None and mpi[index-1] is not None:
            shift = float(mpi[index-1])
            bin_object.setErr(min(0.0, shift), max(0.0, shift), "mpi")
    return result


def _reference_points(measurement: Mapping[str, Any], snapshot: Mapping[str, Any],
                      observable: str) -> Sequence[Mapping[str, Any]] | None:
    if measurement["postprocessor"] == "star_weak_bosons":
        if observable.startswith("Wplus_"):
            channel = snapshot["channels"]["Wplus"]
            return channel["points"] if observable.endswith("_AL") else channel["all"]["points"]
        if observable.startswith("Wminus_"):
            channel = snapshot["channels"]["Wminus"]
            return channel["points"] if observable.endswith("_AL") else channel["all"]["points"]
        if observable == "Zgamma_AL":
            return snapshot["channels"]["Zgamma"]["points"]
    if measurement["postprocessor"] == "phenix_prompt_photon":
        channel = snapshot["channels"].get(observable)
        return channel["points"] if channel else None
    if measurement["postprocessor"] == "star_jet_all":
        dataset = snapshot["datasets"].get(observable)
        return dataset["points"] if dataset else None
    if measurement["postprocessor"] == "hermes_sidis":
        dataset = snapshot.get("datasets", {}).get(observable)
        if dataset and dataset.get("observable") == "A_parallel":
            return dataset["points"]
    if measurement["postprocessor"] in {
        "compass_sidis_a1", "compass_sidis_multiplicity",
        "hermes_sidis_multiplicity", "compass_sidis_charge_ratio",
        "compass_sidis_pt_slope", "sidis_azimuthal_moment",
    }:
        entry = _compass_reference_entry(snapshot, observable)
        if entry is None:
            return None
        dataset, slice_spec = entry
        if slice_spec is None:
            return dataset.get("points")
        if not isinstance(dataset.get("points"), list):
            return None
        return [
            dataset["points"][int(index)]
            for index in slice_spec["flat_bins"]
        ]
    return None


def _comparison_mask(
    measurement: Mapping[str, Any],
    snapshot: Mapping[str, Any],
    observable: str,
    size: int,
) -> list[bool]:
    if measurement["postprocessor"] != "star_jet_all":
        return [True] * size
    try:
        return star_policy.primary_bin_mask(
            measurement, snapshot, observable, size
        )
    except star_policy.ComparisonPolicyError as exc:
        raise CampaignError(str(exc)) from exc


def _masked_bands(
    bands: Mapping[str, Any], mask: Sequence[bool]
) -> dict[str, Any]:
    result = copy.deepcopy(dict(bands))
    for key, values in list(result.items()):
        if isinstance(values, list) and len(values) == len(mask):
            result[key] = [value if keep else None for value, keep in zip(values, mask)]
    return result


def _pp_reference_overlay_points(
    measurement: Mapping[str, Any], snapshot: Mapping[str, Any], plot_stem: str
) -> list[dict[str, Any]] | None:
    def mean_coordinate(point: Mapping[str, Any], name: str) -> float:
        if name in point:
            return float(point[name])
        axis = name.removesuffix("_mean")
        return .5 * (
            float(point[f"{axis}_low"]) + float(point[f"{axis}_high"])
        )

    if plot_stem.startswith("Pull_"):
        return None
    if measurement["postprocessor"] == "star_weak_bosons":
        observables = (
            "Wplus_AL", "Wplus_ALL", "Wminus_AL", "Wminus_ALL", "Zgamma_AL"
        )
    elif measurement["postprocessor"] == "star_jet_all":
        observables = tuple(snapshot.get("datasets", {}))
    elif measurement["postprocessor"] == "hermes_sidis":
        observables = tuple(snapshot.get("datasets", {}))
    elif measurement["postprocessor"] in {
        "compass_sidis_a1", "compass_sidis_multiplicity",
        "hermes_sidis_multiplicity", "compass_sidis_charge_ratio",
        "compass_sidis_pt_slope", "sidis_azimuthal_moment",
    }:
        observables = tuple(snapshot.get("datasets", {})) + tuple(
            Path(str(slice_spec["rivet_path"])).name
            for dataset in snapshot.get("datasets", {}).values()
            for slice_spec in dataset.get("slices", [])
        ) + tuple(snapshot.get("readable_projections", {}))
    else:
        observables = tuple(snapshot.get("channels", {}))
    for observable in observables:
        reference_path = _reference_path(measurement, observable, snapshot)
        if not reference_path or Path(reference_path).name != plot_stem:
            continue
        points = _reference_points(measurement, snapshot, observable)
        if not points:
            return None
        coordinate = "pt_mean"
        if observable in {"Wplus_AL", "Wminus_AL"}:
            coordinate = "eta"
        elif observable in {"Wplus_ALL", "Wminus_ALL"}:
            coordinate = "abs_eta"
        elif observable == "Zgamma_AL":
            coordinate = "coordinate"
        elif measurement["postprocessor"] == "star_jet_all":
            coordinate = "coordinate"
        elif measurement["postprocessor"] == "hermes_sidis":
            dataset = snapshot["datasets"][observable]
            return [
                {
                    **dict(point),
                    "plot_x": (
                        float(point["x"])
                        if dataset["projection"] == "x"
                        else float(point["flat_bin"]) + 0.5
                    ),
                    "value": float(point["aparallel"]),
                    "stat": float(point["aparallel_stat"]),
                }
                for point in points
            ]
        elif measurement["postprocessor"] in {
            "compass_sidis_a1", "compass_sidis_multiplicity",
            "hermes_sidis_multiplicity", "compass_sidis_charge_ratio",
            "compass_sidis_pt_slope", "sidis_azimuthal_moment",
        }:
            entry = _compass_reference_entry(snapshot, observable)
            if entry is None:
                return None
            dataset, slice_spec = entry
            is_a1 = measurement["postprocessor"] == "compass_sidis_a1"
            plotted_coordinate = "z_mean"
            if dataset.get("integrated_projection"):
                plotted_coordinate = f"{dataset['axis']}_mean"
            elif slice_spec is None and dataset.get("projection") in {"x", "z", "pt"}:
                plotted_coordinate = f"{dataset['projection']}_mean"
            elif slice_spec is not None and slice_spec.get("plotted_dimension"):
                plotted_coordinate = f"{slice_spec['plotted_dimension']}_mean"
            elif slice_spec is not None and dataset.get("density_widths") == ["z", "pt2"]:
                plotted_coordinate = "pt2_mean"
            elif slice_spec is not None and dataset.get("density_widths") == ["z", "phperp"]:
                plotted_coordinate = "phperp_mean"
            return [
                {
                    **dict(point),
                    "plot_x": (
                        float(point["x_mean"])
                        if is_a1
                        else (
                            mean_coordinate(point, plotted_coordinate)
                            if (
                                slice_spec is not None
                                or dataset.get("integrated_projection")
                                or dataset.get("projection") in {"x", "z", "pt"}
                            )
                            else float(point["flat_bin"]) + 0.5
                        )
                    ),
                    "value": float(point["a1"] if is_a1 else point["value"]),
                    "stat": float(point["stat"]),
                }
                for point in points
            ]
        mask = _comparison_mask(
            measurement, snapshot, observable, len(points)
        )
        return [
            {**dict(point), "plot_x": float(point[coordinate])}
            for point, keep in zip(points, mask) if keep
        ]
    return None


def _pp_theory_uncertainty_bands(
    measurement: Mapping[str, Any],
    snapshot: Mapping[str, Any],
    summary: Mapping[str, Any],
    plot_stem: str,
) -> Mapping[str, Any] | None:
    uncertainties = summary.get("uncertainties", {})
    if not isinstance(uncertainties, Mapping):
        return None
    for observable, bands in uncertainties.items():
        reference_path = _reference_path(measurement, str(observable), snapshot)
        if reference_path and Path(reference_path).name == plot_stem:
            if not isinstance(bands, Mapping):
                return None
            values = bands.get("monte_carlo", [])
            mask = _comparison_mask(
                measurement, snapshot, str(observable), len(values)
            )
            return _masked_bands(bands, mask)
    return None


def _pulls(
    prediction: Mapping[str, Any],
    points: Sequence[Mapping[str, Any]],
    mask: Sequence[bool] | None = None,
) -> tuple[list[float | None], float, int]:
    values: list[float | None] = []
    chi2 = 0.0
    count = 0
    selected = list(mask) if mask is not None else [True] * len(points)
    if len(selected) != len(points):
        raise CampaignError("pull mask and reference points differ in size")
    for theory, theory_error, point, keep in zip(
        prediction["values"], prediction["errors"], points, selected
    ):
        if not keep:
            values.append(None); continue
        if theory is None or theory_error is None:
            values.append(None); continue
        data_error = math.hypot(float(point["stat"]), float(point["systematic_combined"]))
        denominator = math.hypot(float(theory_error), data_error)
        if denominator <= 0.0:
            values.append(None); continue
        pull = (float(theory)-float(point["value"]))/denominator
        values.append(pull); chi2 += pull*pull; count += 1
    return values, chi2, count


def _star_correlated_goodness_of_fit(
    prediction_set: Mapping[str, Mapping[str, Any]],
    snapshot: Mapping[str, Any],
    measurement: Mapping[str, Any] | None = None,
    *,
    require_complete: bool = True,
) -> dict[str, Any]:
    try:
        import numpy as np
    except ImportError as exc:
        raise CampaignError(
            "NumPy is required for correlated STAR goodness-of-fit calculations"
        ) from exc

    theory: list[float] = []
    data: list[float] = []
    theory_variance: list[float] = []
    retained_indices: list[int] = []
    retained_labels: list[str] = []
    ordering = list(snapshot["primary_covariance"]["ordering"])
    for covariance_index, label in enumerate(ordering):
        dataset_id, point_token = str(label).rsplit(":", 1)
        prediction = prediction_set.get(dataset_id)
        point = snapshot["datasets"][dataset_id]["points"][int(point_token)-1]
        if prediction is None:
            continue
        mask = _comparison_mask(
            measurement or {"postprocessor": "star_jet_all"},
            snapshot,
            dataset_id,
            len(prediction["values"]),
        )
        if not mask[int(point_token)-1]:
            continue
        value = prediction["values"][int(point_token)-1]
        error = prediction["errors"][int(point_token)-1]
        if value is None or error is None:
            continue
        theory.append(float(value))
        data.append(float(point["value"]))
        theory_variance.append(float(error)**2)
        retained_indices.append(covariance_index)
        retained_labels.append(str(label))
    expected_points: int | None = None
    if measurement is not None:
        configured_points = measurement.get(
            "comparison_policy", {}
        ).get("covariance_points")
        if configured_points is not None:
            expected_points = int(configured_points)
        if (
            require_complete
            and expected_points is not None
            and len(retained_indices) != expected_points
        ):
            raise CampaignError(
                f"STAR primary covariance retained {len(retained_indices)} points, "
                f"expected {expected_points}"
            )
        if (
            not require_complete
            and expected_points is not None
            and len(retained_indices) != expected_points
        ):
            return {
                "points": len(retained_indices),
                "expected_points": expected_points,
                "coverage_status": "partial smoke sample",
                "ordering": retained_labels,
                "status": (
                    "goodness of fit not evaluated for partial smoke coverage"
                ),
            }
    if not retained_indices:
        return {
            "points": 0,
            "expected_points": expected_points,
            "coverage_status": (
                "complete"
                if expected_points in {None, 0}
                else "partial smoke sample"
            ),
            "status": "no finite theory bins",
        }

    published = np.asarray(snapshot["primary_covariance"]["covariance"], dtype=float)
    covariance = published[np.ix_(retained_indices, retained_indices)]
    covariance = covariance + np.diag(np.asarray(theory_variance, dtype=float))
    inverse = np.linalg.pinv(covariance, hermitian=True, rcond=1.0e-12)
    residual = np.asarray(theory)-np.asarray(data)
    unprofiled = float(residual @ inverse @ residual)

    global_uncertainties = snapshot["global_uncertainties"]
    lumi = float(global_uncertainties["relative_luminosity_absolute"])
    polarization = float(global_uncertainties["polarization_relative"])
    nuisance_vectors = np.column_stack((
        np.full(len(data), lumi, dtype=float),
        np.asarray(data, dtype=float)*polarization,
    ))
    normal = np.eye(2)+nuisance_vectors.T @ inverse @ nuisance_vectors
    rhs = nuisance_vectors.T @ inverse @ residual
    nuisance = np.linalg.solve(normal, rhs)
    profiled_residual = residual-nuisance_vectors @ nuisance
    profiled = float(
        profiled_residual @ inverse @ profiled_residual+nuisance @ nuisance
    )
    try:
        cholesky = np.linalg.cholesky(covariance)
        decorrelated = np.linalg.solve(cholesky, profiled_residual)
    except np.linalg.LinAlgError:
        decorrelated = np.full(len(data), np.nan)
    return {
        "points": len(data),
        "expected_points": expected_points,
        "coverage_status": (
            "complete"
            if expected_points is None or len(data) == expected_points
            else "partial smoke sample"
        ),
        "ordering": retained_labels,
        "chi2_correlated_without_global_nuisances": unprofiled,
        "chi2_profiled": profiled,
        "nuisance_pulls": {
            "relative_luminosity": float(nuisance[0]),
            "polarization": float(nuisance[1]),
        },
        "decorrelated_pulls": [
            None if not math.isfinite(float(value)) else float(value)
            for value in decorrelated
        ],
        "covariance": (
            "published point-to-point covariance plus diagonal Monte Carlo "
            "statistical variance"
        ),
    }


def _unpolarized_closure(physical: Mapping[str, experimental.BinSeries],
                         unpolarized: experimental.BinSeries) -> dict[str, Any]:
    sigma_uu = experimental.linear_combine_series(
        physical, {label: 0.25 for label in DENOMINATOR})
    values, errors = [], []
    for index in range(len(sigma_uu.values)):
        value, error = experimental.parity_residual(
            sigma_uu.values[index], sigma_uu.variances[index],
            unpolarized.values[index], unpolarized.variances[index])
        values.append(value); errors.append(error)
    return {"edges": list(sigma_uu.edges), "values": values, "errors": errors}


def _sidis_charge_objects(dataset: Mapping[str, Any]) -> tuple[str, str, str]:
    target = str(dataset["target"])
    species = str(dataset["species"])
    if species == "pi_charge_difference":
        positive, negative, covariance = "piplus", "piminus", "pi"
    elif species == "k_charge_difference":
        positive, negative, covariance = "kplus", "kminus", "k"
    elif species == "h_charge_difference":
        positive, negative, covariance = "hplus", "hminus", "h"
    else:
        raise CampaignError(f"Unknown HERMES charge-difference species {species}")
    return (
        f"Yield_{target}_{positive}_x",
        f"Yield_{target}_{negative}_x",
        f"CovarianceProxy_{target}_{covariance}_charge_difference",
    )


def _sidis_prediction_sets(
    groups: Mapping[
        tuple[Any, ...],
        Mapping[str, Mapping[str, Sequence[Mapping[str, Any]]]],
    ],
    campaign_dir: Path,
    measurement: Mapping[str, Any],
    snapshot: Mapping[str, Any],
    include_diagnostics: bool,
) -> dict[tuple[Any, ...], dict[str, Any]]:
    """Construct every normalized-yield SIDIS variation."""

    analysis = str(measurement["analysis"]["name"])
    variations = sorted(
        {_sidis_variation_key(key) for key in groups},
        key=lambda item: tuple(str(value) for value in item),
    )
    cache: dict[tuple[str, str, str], experimental.BinSeries] = {}
    output: dict[tuple[Any, ...], dict[str, Any]] = {}
    for variation in variations:
        prediction_set: dict[str, Any] = {}
        for dataset_id, dataset in snapshot["datasets"].items():
            if dataset["projection"] == "difference":
                if not include_diagnostics:
                    continue
                positive_name, negative_name, covariance_name = (
                    _sidis_charge_objects(dataset)
                )
                positive = _sidis_component_samples(
                    groups,
                    variation,
                    campaign_dir,
                    analysis,
                    positive_name,
                    str(dataset["target"]),
                    cache,
                )
                negative = _sidis_component_samples(
                    groups,
                    variation,
                    campaign_dir,
                    analysis,
                    negative_name,
                    str(dataset["target"]),
                    cache,
                )
                covariance = _sidis_component_samples(
                    groups,
                    variation,
                    campaign_dir,
                    analysis,
                    covariance_name,
                    str(dataset["target"]),
                    cache,
                )
                differences = {
                    label: sidis.difference_series(
                        positive[label], negative[label], covariance[label]
                    )
                    for label in positive
                }
                result = sidis.asymmetry(
                    differences, str(dataset["target"])
                )
                observable = f"AParallelChargeDifference_{dataset_id}"
            else:
                object_name = (
                    f"Yield_{dataset['target']}_{dataset['species']}_"
                    f"{dataset['projection']}"
                )
                samples = _sidis_component_samples(
                    groups,
                    variation,
                    campaign_dir,
                    analysis,
                    object_name,
                    str(dataset["target"]),
                    cache,
                )
                result = sidis.asymmetry(samples, str(dataset["target"]))
                observable = str(dataset_id)
            prediction_set[observable] = {
                "edges": result["edges"],
                "values": result["values"],
                "errors": result["errors"],
            }
            if include_diagnostics:
                for label, closure in result["closures"].items():
                    prediction_set[f"{label}_{dataset_id}"] = closure

        if include_diagnostics:
            for diagnostic_id, diagnostic in snapshot["diagnostics"].items():
                prefix = (
                    f"{diagnostic['target']}_{diagnostic['species']}"
                )
                numerators = _sidis_component_samples(
                    groups,
                    variation,
                    campaign_dir,
                    analysis,
                    f"CosPhiNumerator_{prefix}",
                    str(diagnostic["target"]),
                    cache,
                )
                denominators = _sidis_component_samples(
                    groups,
                    variation,
                    campaign_dir,
                    analysis,
                    f"CosPhiDenominator_{prefix}",
                    str(diagnostic["target"]),
                    cache,
                )
                positive = _sidis_component_samples(
                    groups,
                    variation,
                    campaign_dir,
                    analysis,
                    f"CosPhiCovariancePositive_{prefix}",
                    str(diagnostic["target"]),
                    cache,
                )
                negative = _sidis_component_samples(
                    groups,
                    variation,
                    campaign_dir,
                    analysis,
                    f"CosPhiCovarianceNegative_{prefix}",
                    str(diagnostic["target"]),
                    cache,
                )
                covariance = {
                    label: [
                        positive[label].variances[index]
                        - negative[label].variances[index]
                        for index in range(len(positive[label].values))
                    ]
                    for label in numerators
                }
                prediction_set[f"UnpolarizedCosPhi_{diagnostic_id}"] = (
                    sidis.azimuthal_moment(
                        numerators,
                        denominators,
                        covariance,
                        str(diagnostic["target"]),
                    )
                )
        output[variation] = prediction_set
    return output


def _compass_component_samples(
    groups: Mapping[
        tuple[Any, ...],
        Mapping[str, Mapping[str, Sequence[Mapping[str, Any]]]],
    ],
    variation: tuple[Any, ...],
    campaign_dir: Path,
    analysis: str,
    object_name: str,
    helicities: Sequence[str],
    cache: dict[tuple[str, str, str], experimental.BinSeries],
    components: Sequence[str] = ("P", "N"),
) -> dict[str, experimental.BinSeries]:
    """Load configured target samples after normalized signed-NLO addition."""

    family, channel, polarized, unpolarized, scale, mpi = variation
    samples: dict[str, experimental.BinSeries] = {}
    for component in components:
        key = (
            family,
            channel,
            component,
            polarized,
            unpolarized,
            scale,
            mpi,
        )
        if key not in groups:
            raise CampaignError(
                f"Missing SIDIS target component {component} for "
                f"{_variation_id(variation)}"
            )
        loaded = _load_sidis_object(
            groups[key],
            campaign_dir,
            analysis,
            object_name,
            cache,
            helicities,
        )
        for helicity, series in loaded.items():
            label = component if helicities == ("00",) else f"{component}:{helicity}"
            samples[label] = series
    return samples


def _compass_sidis_prediction_sets(
    groups: Mapping[
        tuple[Any, ...],
        Mapping[str, Mapping[str, Sequence[Mapping[str, Any]]]],
    ],
    campaign_dir: Path,
    measurement: Mapping[str, Any],
    snapshot: Mapping[str, Any],
) -> dict[tuple[Any, ...], dict[str, Any]]:
    """Construct all COMPASS A1 or multiplicity variations."""

    analysis = str(measurement["analysis"]["name"])
    variations = sorted(
        {_sidis_variation_key(key) for key in groups},
        key=lambda item: tuple(str(value) for value in item),
    )
    cache: dict[tuple[str, str, str], experimental.BinSeries] = {}
    output: dict[tuple[Any, ...], dict[str, Any]] = {}
    is_a1 = measurement["postprocessor"] == "compass_sidis_a1"
    is_ratio = measurement["postprocessor"] == "compass_sidis_charge_ratio"
    helicities = compass_sidis.HELICITIES if is_a1 else ("00",)
    schema6 = int(measurement.get("schema_version", 5)) >= 6
    postprocess_config = measurement.get("postprocess_config", {})
    target_outputs = postprocess_config.get(
        "target_outputs", {"D": {"P": .5, "N": .5}}
    )
    for variation in variations:
        prediction_set: dict[str, Any] = {}
        for species, dataset in snapshot["datasets"].items():
            raw = dataset["raw_objects"]
            published_target = dataset.get("published_target")
            if published_target is None:
                if len(target_outputs) != 1:
                    raise CampaignError(
                        f"{measurement['id']}/{species} must name its published_target"
                    )
                published_target = next(iter(target_outputs))
            if published_target not in target_outputs:
                raise CampaignError(
                    f"{measurement['id']}/{species} has unknown target output {published_target}"
                )
            target_weights = {
                str(component): float(weight)
                for component, weight in target_outputs[published_target].items()
            }
            components = tuple(target_weights)
            if is_a1:
                ordinary = _compass_component_samples(
                    groups, variation, campaign_dir, analysis,
                    str(raw["ordinary"]), helicities, cache, components,
                )
                inverse_d = _compass_component_samples(
                    groups, variation, campaign_dir, analysis,
                    str(raw["inverse_depolarization"]), helicities, cache,
                    components,
                )
                covariance = _compass_component_samples(
                    groups, variation, campaign_dir, analysis,
                    str(raw["covariance"]), helicities, cache, components,
                )
                if schema6:
                    result = compass_sidis.a1_target_combination(
                        ordinary, inverse_d, covariance, target_weights,
                        float(postprocess_config["longitudinal_target_scale"]),
                    )
                else:
                    result = compass_sidis.a1_deuteron(
                        ordinary, inverse_d, covariance
                    )
            else:
                numerator = _compass_component_samples(
                    groups, variation, campaign_dir, analysis,
                    str(raw["numerator"]), helicities, cache, components,
                )
                denominator = _compass_component_samples(
                    groups, variation, campaign_dir, analysis,
                    str(raw["denominator"]), helicities, cache, components,
                )
                covariance = _compass_component_samples(
                    groups, variation, campaign_dir, analysis,
                    str(raw["covariance"]), helicities, cache, components,
                )
                widths = [
                    math.prod(
                        float(point[f"{axis}_high"])
                        - float(point[f"{axis}_low"])
                        for axis in dataset.get("density_widths", ["z"])
                    )
                    for point in dataset["points"]
                ]
                if schema6 and is_ratio:
                    result = compass_sidis.charge_ratio_target_combination(
                        numerator, denominator, covariance, target_weights,
                    )
                elif schema6:
                    result = compass_sidis.multiplicity_target_combination(
                        numerator, denominator, covariance,
                        target_weights, widths,
                    )
                else:
                    result = compass_sidis.multiplicity_isoscalar(
                        numerator, denominator, covariance, widths,
                    )
            from compass_multiplicity_fiducial import unsupported_cells
            unsupported = unsupported_cells(snapshot, dataset, str(species), target_weights)
            bin_status = ["ok" if value is not None else "insufficient_mc_support"
                          for value in result["values"]]
            for index in unsupported:
                result["values"][index] = None
                result["errors"][index] = None
                bin_status[index] = "outside_nominal_beam_support"
            prediction_set[str(species)] = {
                "bin_status": bin_status,
                "unsupported_nominal_beam_bins": unsupported,
                "edges": result["edges"],
                "values": result["values"],
                "errors": result["errors"],
                "published_target": str(published_target),
            }
            if not is_a1:
                for slice_spec in dataset["slices"]:
                    observable = Path(str(slice_spec["rivet_path"])).name
                    indices = [int(index) for index in slice_spec["flat_bins"]]
                    prediction_set[observable] = {
                        "edges": [
                            float(value) for value in
                            slice_spec.get("edges", slice_spec.get("z_edges", []))
                        ],
                        "values": [result["values"][index] for index in indices],
                        "errors": [result["errors"][index] for index in indices],
                        "parent_dataset": str(species),
                        "flat_bins": indices,
                        "bin_status": [bin_status[index] for index in indices],
                        "unsupported_nominal_beam_bins": [
                            i for i, index in enumerate(indices) if index in unsupported],
                    }
        # HERMES publishes one-dimensional projections of the five independent
        # 3D binnings.  Load the dedicated event-aggregated raw objects so these
        # curves are ratios of integrated yields, never averages/projections of
        # already formed cell multiplicities.
        if not is_a1:
            for observable, projection in snapshot.get(
                "readable_projections", {}
            ).items():
                published_target = str(projection["published_target"])
                if published_target not in target_outputs:
                    raise CampaignError(
                        f"{measurement['id']}/{observable} has unknown target "
                        f"output {published_target}"
                    )
                target_weights = {
                    str(component): float(weight)
                    for component, weight in
                    target_outputs[published_target].items()
                }
                components = tuple(target_weights)
                raw = projection["raw_objects"]
                numerator = _compass_component_samples(
                    groups, variation, campaign_dir, analysis,
                    str(raw["numerator"]), helicities, cache, components,
                )
                denominator = _compass_component_samples(
                    groups, variation, campaign_dir, analysis,
                    str(raw["denominator"]), helicities, cache, components,
                )
                covariance = _compass_component_samples(
                    groups, variation, campaign_dir, analysis,
                    str(raw["covariance"]), helicities, cache, components,
                )
                widths = [
                    math.prod(
                        float(point[f"{axis}_high"])
                        - float(point[f"{axis}_low"])
                        for axis in projection["density_widths"]
                    )
                    for point in projection["points"]
                ]
                result = compass_sidis.multiplicity_target_combination(
                    numerator, denominator, covariance,
                    target_weights, widths,
                )
                prediction_set[str(observable)] = {
                    "edges": result["edges"],
                    "values": result["values"],
                    "errors": result["errors"],
                    "published_target": published_target,
                    "integrated_projection": True,
                }
        output[variation] = prediction_set
    return output


def _sidis_diagnostic_prediction_sets(
    groups: Mapping[
        tuple[Any, ...],
        Mapping[str, Mapping[str, Sequence[Mapping[str, Any]]]],
    ],
    campaign_dir: Path,
    measurement: Mapping[str, Any],
    snapshot: Mapping[str, Any],
) -> dict[tuple[Any, ...], dict[str, Any]]:
    """Construct target-combined pT-slope or azimuthal predictions."""

    analysis = str(measurement["analysis"]["name"])
    variations = sorted(
        {_sidis_variation_key(key) for key in groups},
        key=lambda item: tuple(str(value) for value in item),
    )
    cache: dict[tuple[str, str, str], experimental.BinSeries] = {}
    target_outputs = measurement["postprocess_config"]["target_outputs"]
    output: dict[tuple[Any, ...], dict[str, Any]] = {}
    for variation in variations:
        prediction_set: dict[str, Any] = {}
        for dataset_id, dataset in snapshot["datasets"].items():
            published_target = str(dataset["published_target"])
            if published_target not in target_outputs:
                raise CampaignError(
                    f"{measurement['id']}/{dataset_id} has unknown target "
                    f"output {published_target}"
                )
            target_weights = {
                str(component): float(weight)
                for component, weight in target_outputs[published_target].items()
            }
            components = tuple(target_weights)
            raw = dataset["raw_objects"]
            if measurement["postprocessor"] == "compass_sidis_pt_slope":
                spectra = _compass_component_samples(
                    groups, variation, campaign_dir, analysis,
                    str(raw["spectrum"]), ("00",), cache, components,
                )
                spectrum_binning = dataset["spectrum_binning"]
                result = compass_sidis.pt2_slope_target_combination(
                    spectra,
                    target_weights,
                    [float(value) for value in spectrum_binning["pt2_edges"]],
                    len(dataset["points"]),
                    _compass_component_samples(
                        groups, variation, campaign_dir, analysis,
                        str(raw["covariance"]), ("00",), cache, components),
                    float(snapshot["fit_policy"]["minimum_region_significance"]),
                )
            elif "inputs" in raw:
                inputs = _compass_component_samples(
                    groups, variation, campaign_dir, analysis,
                    str(raw["inputs"]), ("00",), cache, components)
                covariance = _compass_component_samples(
                    groups, variation, campaign_dir, analysis,
                    str(raw["covariance"]), ("00",), cache, components)
                result = compass_sidis.azimuthal_target_fit(
                    inputs, covariance, target_weights, int(dataset["harmonic"]),
                    len(dataset["points"]), snapshot["estimator"])
            else:
                numerator = _compass_component_samples(
                    groups, variation, campaign_dir, analysis,
                    str(raw["numerator"]), ("00",), cache, components,
                )
                denominator = _compass_component_samples(
                    groups, variation, campaign_dir, analysis,
                    str(raw["denominator"]), ("00",), cache, components,
                )
                covariance_positive = _compass_component_samples(
                    groups, variation, campaign_dir, analysis,
                    str(raw["covariance_positive"]), ("00",), cache, components,
                )
                covariance_negative = _compass_component_samples(
                    groups, variation, campaign_dir, analysis,
                    str(raw["covariance_negative"]), ("00",), cache, components,
                )
                result = compass_sidis.azimuthal_target_combination(
                    numerator,
                    denominator,
                    covariance_positive,
                    covariance_negative,
                    target_weights,
                )
            prediction_set[str(dataset_id)] = {
                **result,
                "published_target": published_target,
            }
            for slice_spec in dataset.get("slices", []):
                observable = Path(str(slice_spec["rivet_path"])).name
                indices = [int(index) for index in slice_spec["flat_bins"]]
                prediction_set[observable] = {
                    "edges": [float(value) for value in slice_spec["edges"]],
                    "values": [result["values"][index] for index in indices],
                    "errors": [result["errors"][index] for index in indices],
                    "published_target": published_target,
                    "parent_dataset": str(dataset_id),
                    "flat_bins": indices,
                }
        output[variation] = prediction_set
    return output


def _sidis_pull_bins(
    prediction: Mapping[str, Any], dataset: Mapping[str, Any]
) -> list[float | None]:
    """Return pulls in the full Rivet flat-bin layout."""

    pulls: list[float | None] = [None] * len(prediction["values"])
    for point_index, point in enumerate(dataset["points"]):
        flat_bin = int(point.get("flat_bin", point_index))
        theory = prediction["values"][flat_bin]
        theory_error = prediction["errors"][flat_bin]
        if theory is None or theory_error is None:
            continue
        total = math.sqrt(
            float(theory_error) ** 2
            + float(point["aparallel_stat"]) ** 2
            + float(point["aparallel_systematic"]) ** 2
        )
        if total > 0.0:
            pulls[flat_bin] = (
                float(theory) - float(point["aparallel"])
            ) / total
    return pulls


def _primary_reference_observable(
    measurement: Mapping[str, Any],
    snapshot: Mapping[str, Any],
    observable: str,
) -> bool:
    reference_path = _reference_path(measurement, observable, snapshot)
    if reference_path is None:
        return False
    if measurement["postprocessor"] == "star_jet_all":
        dataset = snapshot["datasets"].get(observable)
        return bool(dataset) and not bool(dataset.get("alternate_projection"))
    if measurement["postprocessor"] == "hermes_sidis":
        dataset = snapshot["datasets"].get(observable)
        return bool(dataset) and dataset.get("observable") == "A_parallel"
    if measurement["postprocessor"] in {
        "compass_sidis_a1", "compass_sidis_multiplicity",
        "hermes_sidis_multiplicity", "compass_sidis_charge_ratio",
        "compass_sidis_pt_slope", "sidis_azimuthal_moment",
    }:
        entry = _compass_reference_entry(snapshot, observable)
        return entry is not None and bool(entry[0].get("data_available", True))
    return True


def postprocess_sidis_diagnostic(
    args: argparse.Namespace, measurement: Mapping[str, Any]
) -> Path:
    """Postprocess the lower-dimensional pT and azimuthal SIDIS diagnostics."""

    campaign_dir = _campaign_dir(measurement["id"], args.tag)
    manifest_path = campaign_dir / experimental.MANIFEST_NAME
    if not manifest_path.exists():
        raise CampaignError(f"No prepared campaign at {campaign_dir}")
    manifest = _load_json(manifest_path)
    _assert_manifest_signature_current(manifest, measurement)
    groups = _logical_sidis_groups(manifest, campaign_dir, measurement)
    prediction_path = campaign_dir / "postprocess" / "prediction.yoda"
    if args.dry_run:
        print(json.dumps({
            "measurement": measurement["id"],
            "target_variation_groups": len(groups),
            "postprocessor": measurement["postprocessor"],
            "output": str(prediction_path),
        }, indent=2))
        return prediction_path

    snapshot = _measurement_snapshot(measurement)
    predictions = _sidis_diagnostic_prediction_sets(
        groups, campaign_dir, measurement, snapshot
    )
    nominal_mpi = measurement["families"]["nominal"]["mpi"]
    central_key = ("nominal", "sidis", 0, 0, 1.0, nominal_mpi)
    if central_key not in predictions:
        raise CampaignError(
            f"Missing central SIDIS diagnostic prediction {central_key}"
        )
    central = predictions[central_key]
    bands = aggregate_uncertainties(predictions, measurement)
    yoda = experimental._import_yoda()
    objects: list[Any] = []
    for observable, prediction in central.items():
        reference_path = _reference_path(measurement, observable, snapshot)
        if reference_path is None:
            continue
        target_output = str(prediction.get("published_target", "D"))
        target_weights = measurement["postprocess_config"]["target_outputs"].get(
            target_output
        )
        observable_definition = (
            "target-combined low-pT2 yield fitted to A exp(-pT2/<pT2>)"
            if measurement["postprocessor"] == "compass_sidis_pt_slope"
            else (
                "16-bin harmonic GLS fit excluding phi bins 0/15, divided by full-phi mean epsilon_n"
                if measurement["id"].startswith("COMPASS_")
                else "target-combined sum[cos(n phi)] / identified-hadron yield"
            )
        )
        objects.append(_estimate_with_bands(
            yoda,
            prediction,
            reference_path,
            {
                "Generator": "HerwigPol POWHEG NLO+PS",
                "HardProcessAccuracy": "NLO",
                "NLOCombination": "normalized POSNLO+NEGNLO bins",
                "TargetCombination": json.dumps(target_weights, sort_keys=True),
                "ObservableDefinition": observable_definition,
                "DiagnosticOnly": 1,
                "ExperimentalCorrectionsAppliedToHerwig": "none",
            },
            bands.get(observable),
        ))

    flat_predictions = {
        dataset_id: central[dataset_id]
        for dataset_id in snapshot["datasets"]
    }
    if any(
        dataset.get("data_available", True)
        for dataset in snapshot["datasets"].values()
    ):
        goodness: Mapping[str, Any] = (
            compass_sidis.diagonal_multiplicity_goodness_of_fit(
                flat_predictions, snapshot
            )
        )
    else:
        goodness = {
            "status": "not evaluated: official numerical release is not vendored",
            "policy": "no pseudo-data or covariance entries are fabricated",
        }

    rows: list[dict[str, Any]] = []
    for dataset_id, dataset in snapshot["datasets"].items():
        prediction = central[dataset_id]
        points = dataset.get("points")
        if not isinstance(points, list):
            points = [{} for _ in prediction["values"]]
        for index, point in enumerate(points):
            rows.append({
                "measurement": measurement["id"],
                "dataset": dataset_id,
                "published_target": dataset.get("published_target"),
                "harmonic": dataset.get("harmonic"),
                "projection": dataset.get("projection"),
                "bin": index + 1,
                "flat_bin": int(point.get("flat_bin", index)),
                "x_low": point.get("x_low"),
                "x_high": point.get("x_high"),
                "q2_low": point.get("q2_low"),
                "q2_high": point.get("q2_high"),
                "y_low": point.get("y_low"),
                "y_high": point.get("y_high"),
                "z_low": point.get("z_low"),
                "z_high": point.get("z_high"),
                "pt_low": point.get("pt_low"),
                "pt_high": point.get("pt_high"),
                "theory": prediction["values"][index],
                "mc_stat": prediction["errors"][index],
                "data": point.get("value"),
                "data_stat": point.get("stat"),
                "data_systematic": point.get("systematic"),
                "fit_status": prediction.get("fit_diagnostics", [{}]*len(points))[index].get("status"),
            })

    output_dir = campaign_dir / "postprocess"
    output_dir.mkdir(parents=True, exist_ok=True)
    experimental._write_yoda_objects(yoda, objects, prediction_path)
    if not experimental._nonempty(prediction_path):
        raise CampaignError(f"Postprocessing produced empty output {prediction_path}")
    summary = {
        "measurement": measurement["id"],
        "tag": args.tag,
        "diagnostic_only": True,
        "hard_process_accuracy": "POWHEG NLO+PS",
        "central_sample_label": measurement["families"]["nominal"]["label"],
        "nlo_combination": (
            "shards combined within each contribution, then normalized "
            "POSNLO and NEGNLO bins added"
        ),
        "raw_object_inventory": {
            dataset_id: dict(dataset["raw_objects"])
            for dataset_id, dataset in snapshot["datasets"].items()
        },
        "prediction_object_count": len(objects),
        "uncertainties": bands,
        "goodness_of_fit": goodness,
        "masked_bins": {
            dataset_id: [
                index + 1
                for index, value in enumerate(central[dataset_id]["values"])
                if value is None
            ]
            for dataset_id in snapshot["datasets"]
        },
        "systematic_model": snapshot.get(
            "systematics", "not available because numerical data are not vendored"
        ),
        "reference_provenance": snapshot["provenance"],
        "variations": {
            _variation_id(key): {
                dataset_id: prediction_set[dataset_id]
                for dataset_id in snapshot["datasets"]
            }
            for key, prediction_set in predictions.items()
        },
    }
    experimental.atomic_write_json(output_dir / "summary.json", summary)
    with (output_dir / "central.csv").open(
        "w", encoding="utf-8", newline=""
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    manifest["postprocess"] = {
        "created_at": experimental.utc_now(),
        "prediction": str(prediction_path.relative_to(campaign_dir)),
        "predictions": [{
            "family": "nominal",
            "label": measurement["families"]["nominal"]["label"],
            "path": str(prediction_path.relative_to(campaign_dir)),
        }],
        "summary": "postprocess/summary.json",
        "central_csv": "postprocess/central.csv",
    }
    manifest["updated_at"] = experimental.utc_now()
    manifest["history"].append(
        {"at": experimental.utc_now(), "action": "postprocess"}
    )
    experimental.atomic_write_json(manifest_path, manifest)
    print(f"Wrote normalized SIDIS diagnostic predictions to {prediction_path}")
    return prediction_path


def postprocess_compass_sidis(
    args: argparse.Namespace, measurement: Mapping[str, Any]
) -> Path:
    """Postprocess one schema-5/6 SIDIS campaign at normalized-bin level."""

    campaign_dir = _campaign_dir(measurement["id"], args.tag)
    manifest_path = campaign_dir / experimental.MANIFEST_NAME
    if not manifest_path.exists():
        raise CampaignError(f"No prepared campaign at {campaign_dir}")
    manifest = _load_json(manifest_path)
    _assert_manifest_signature_current(manifest, measurement)
    groups = _logical_sidis_groups(manifest, campaign_dir, measurement)
    prediction_path = campaign_dir / "postprocess" / "prediction.yoda"
    if args.dry_run:
        print(
            json.dumps(
                {
                    "measurement": measurement["id"],
                    "target_variation_groups": len(groups),
                    "postprocessor": measurement["postprocessor"],
                    "output": str(prediction_path),
                },
                indent=2,
            )
        )
        return prediction_path

    snapshot = _measurement_snapshot(measurement)
    predictions = _compass_sidis_prediction_sets(
        groups, campaign_dir, measurement, snapshot
    )
    nominal_mpi = measurement["families"]["nominal"]["mpi"]
    central_key = ("nominal", "sidis", 0, 0, 1.0, nominal_mpi)
    if central_key not in predictions:
        raise CampaignError(f"Missing central COMPASS SIDIS prediction {central_key}")
    central = predictions[central_key]
    bands = aggregate_uncertainties(predictions, measurement)
    yoda = experimental._import_yoda()
    objects: list[Any] = []
    is_a1 = measurement["postprocessor"] == "compass_sidis_a1"
    is_ratio = measurement["postprocessor"] == "compass_sidis_charge_ratio"

    for observable, prediction in central.items():
        reference_path = _reference_path(measurement, observable, snapshot)
        if reference_path is None:
            continue
        target_output = prediction.get("published_target", "D")
        target_weights = measurement.get("postprocess_config", {}).get(
            "target_outputs", {}
        ).get(target_output)
        annotations = {
            "Generator": "HerwigPol POWHEG NLO+PS",
            "HardProcessAccuracy": "NLO",
            "NLOCombination": "normalized POSNLO+NEGNLO bins",
            "TargetCombination": json.dumps(
                target_weights,
                sort_keys=True,
            ) if target_weights is not None else (
                "sigma_UU=(p+n)/2; sigma_LL=0.925*(p+n)/2"
                if is_a1 else "P/N isoscalar sum before numerator/DIS ratio"
            ),
            "ObservableDefinition": (
                "helicity-signed inverse-D yield divided by ordinary yield"
                if is_a1
                else (
                    "negative identified-hadron yield divided by positive "
                    "identified-hadron yield"
                    if is_ratio else
                    "hadron yield per inclusive-DIS event and published density widths"
                )
            ),
            "UnsupportedNominalBeamBins": json.dumps(
                [i+1 for i in prediction.get("unsupported_nominal_beam_bins", [])]),
            "ExperimentalCorrectionsAppliedToHerwig": "none",
        }
        objects.append(
            _estimate_with_bands(
                yoda,
                prediction,
                reference_path,
                annotations,
                bands.get(observable),
            )
        )

    flat_predictions = {
        species: central[species] for species in snapshot["datasets"]
    }
    if is_a1:
        goodness = compass_sidis.a1_goodness_of_fit(
            flat_predictions, snapshot
        )
    elif is_ratio:
        goodness = compass_sidis.charge_ratio_goodness_of_fit(
            flat_predictions, snapshot
        )
    elif measurement["postprocessor"] == "hermes_sidis_multiplicity":
        goodness = compass_sidis.hermes_multiplicity_goodness_of_fit(
            flat_predictions, snapshot
        )
    elif any(
        "systematic_correlated_80pct" in dataset["points"][0]
        for dataset in snapshot["datasets"].values()
        if dataset["points"]
    ):
        goodness = compass_sidis.multiplicity_goodness_of_fit(
            flat_predictions, snapshot
        )
    else:
        goodness = compass_sidis.diagonal_multiplicity_goodness_of_fit(
            flat_predictions, snapshot
        )

    rows: list[dict[str, Any]] = []
    for species, dataset in snapshot["datasets"].items():
        prediction = central[species]
        for index, point in enumerate(dataset["points"]):
            correction = point.get("corrections", {})
            rows.append(
                {
                    "measurement": measurement["id"],
                    "species": species,
                    "published_target": dataset.get("published_target"),
                    "binning": dataset.get("binning"),
                    "bin": index + 1,
                    "flat_bin": int(point.get("flat_bin", index)),
                    "slice": point.get("slice"),
                    "slice_bin": point.get("slice_bin"),
                    "x_low": point.get("x_low"),
                    "x_high": point.get("x_high"),
                    "x_mean": point.get("x_mean"),
                    "y_low": point.get("y_low"),
                    "y_high": point.get("y_high"),
                    "y_mean": point.get("y_mean"),
                    "q2_low": point.get("q2_low"),
                    "q2_high": point.get("q2_high"),
                    "q2_mean": point.get("q2_mean"),
                    "z_low": point.get("z_low"),
                    "z_high": point.get("z_high"),
                    "z_mean": point.get("z_mean"),
                    "pt2_low": point.get("pt2_low"),
                    "pt2_high": point.get("pt2_high"),
                    "pt2_mean": point.get("pt2_mean"),
                    "phperp_low": point.get("phperp_low"),
                    "phperp_high": point.get("phperp_high"),
                    "phperp_mean": point.get("phperp_mean"),
                    "momentum_low": point.get("momentum_low"),
                    "momentum_high": point.get("momentum_high"),
                    "momentum_mean": point.get("momentum_mean"),
                    "theory": prediction["values"][index],
                    "mc_stat": prediction["errors"][index],
                    "theory_status": prediction["bin_status"][index],
                    "data": point.get("a1", point.get("value")),
                    "data_stat": point.get("stat"),
                    "data_systematic": point.get("systematic"),
                    "radiative_hadron": correction.get("radiative_hadron"),
                    "radiative_dis": correction.get("radiative_dis"),
                    "dvm_hadron": correction.get("dvm_hadron"),
                    "dvm_dis": correction.get("dvm_dis"),
                }
            )

    output_dir = campaign_dir / "postprocess"
    output_dir.mkdir(parents=True, exist_ok=True)
    experimental._write_yoda_objects(yoda, objects, prediction_path)
    if not experimental._nonempty(prediction_path):
        raise CampaignError(f"Postprocessing produced empty output {prediction_path}")
    summary = {
        "measurement": measurement["id"],
        "tag": args.tag,
        "hard_process_accuracy": "POWHEG NLO+PS",
        "central_sample_label": measurement["families"]["nominal"]["label"],
        "nlo_combination": (
            "shards combined within each contribution, then normalized "
            "POSNLO and NEGNLO bins added"
        ),
        "raw_object_inventory": {
            species: dict(dataset["raw_objects"])
            for species, dataset in snapshot["datasets"].items()
        },
        "prediction_object_count": len(objects),
        "uncertainties": bands,
        "correlated_goodness_of_fit": goodness,
        "masked_bins": {
            species: [
                index + 1
                for index, value in enumerate(central[species]["values"])
                if value is None
            ]
            for species in snapshot["datasets"]
        },
        "unsupported_nominal_beam_bins": {
            species: [i + 1 for i in central[species]["unsupported_nominal_beam_bins"]]
            for species in snapshot["datasets"]},
        "systematic_model": snapshot["systematics"],
        "correction_policy": snapshot.get("corrections"),
        "reference_provenance": snapshot["provenance"],
        "variations": {
            _variation_id(key): {
                species: prediction_set[species]
                for species in snapshot["datasets"]
            }
            for key, prediction_set in predictions.items()
        },
    }
    experimental.atomic_write_json(output_dir / "summary.json", summary)
    with (output_dir / "central.csv").open(
        "w", encoding="utf-8", newline=""
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    manifest["postprocess"] = {
        "created_at": experimental.utc_now(),
        "prediction": str(prediction_path.relative_to(campaign_dir)),
        "predictions": [
            {
                "family": "nominal",
                "label": measurement["families"]["nominal"]["label"],
                "path": str(prediction_path.relative_to(campaign_dir)),
            }
        ],
        "summary": "postprocess/summary.json",
        "central_csv": "postprocess/central.csv",
    }
    manifest["updated_at"] = experimental.utc_now()
    manifest["history"].append(
        {"at": experimental.utc_now(), "action": "postprocess"}
    )
    experimental.atomic_write_json(manifest_path, manifest)
    print(f"Wrote normalized SIDIS predictions to {prediction_path}")
    return prediction_path


def postprocess_sidis(
    args: argparse.Namespace, measurement: Mapping[str, Any]
) -> Path:
    campaign_dir = _campaign_dir(measurement["id"], args.tag)
    manifest_path = campaign_dir / experimental.MANIFEST_NAME
    if not manifest_path.exists():
        raise CampaignError(f"No prepared campaign at {campaign_dir}")
    manifest = _load_json(manifest_path)
    _assert_manifest_signature_current(manifest, measurement)
    groups = _logical_sidis_groups(manifest, campaign_dir, measurement)
    prediction_path = campaign_dir / "postprocess" / "prediction.yoda"
    if args.dry_run:
        print(
            json.dumps(
                {
                    "measurement": measurement["id"],
                    "target_variation_groups": len(groups),
                    "output": str(prediction_path),
                },
                indent=2,
            )
        )
        return prediction_path

    snapshot = _measurement_snapshot(measurement)
    include_diagnostics = bool(
        getattr(args, "include_diagnostics", False)
    )
    predictions = _sidis_prediction_sets(
        groups,
        campaign_dir,
        measurement,
        snapshot,
        include_diagnostics,
    )
    nominal_mpi = measurement["families"]["nominal"]["mpi"]
    central_key = ("nominal", "sidis", 0, 0, 1.0, nominal_mpi)
    if central_key not in predictions:
        raise CampaignError(f"Missing central SIDIS prediction {central_key}")
    central = predictions[central_key]
    bands = aggregate_uncertainties(predictions, measurement)
    yoda = experimental._import_yoda()
    objects: list[Any] = []
    rows: list[dict[str, Any]] = []
    goodness: dict[str, Any] = {}
    pulls: dict[str, Any] = {}

    for observable, prediction in central.items():
        reference_path = _reference_path(
            measurement, observable, snapshot
        )
        dataset = snapshot["datasets"].get(observable)
        if dataset is None and observable.startswith(
            "AParallelChargeDifference_"
        ):
            dataset = snapshot["datasets"].get(
                observable.removeprefix("AParallelChargeDifference_")
            )
        if dataset is None and observable.startswith("UnpolarizedCosPhi_"):
            dataset = snapshot["diagnostics"].get(
                observable.removeprefix("UnpolarizedCosPhi_")
            )
        target = str(dataset.get("target", "unknown")) if dataset else "unknown"
        path = reference_path or f"/{measurement['analysis']['name']}/DIAGNOSTICS/{observable}"
        annotations = {
            "Generator": "HerwigPol POWHEG NLO+PS",
            "HardProcessAccuracy": "NLO",
            "PDFInputs": "NNPDF40 NLO + NNPDFpol2.0 NLO",
            "SampleLabel": measurement["families"]["nominal"]["label"],
            "HelicityCombination": "independent PP,PM,MP,MM samples",
            "NLOCombination": "normalized POSNLO+NEGNLO bins",
            "TargetModel": (
                "proton"
                if target == "proton"
                else "deuteron impulse approximation"
            ),
            "ChargeDifferenceData": (
                "published as A1; prediction retained as auxiliary A_parallel"
                if observable.startswith("AParallelChargeDifference_")
                else "not applicable"
            ),
            "ObservableDefinition": (
                "unpolarized helicity-summed yield 2<cos(phi)>; not the "
                "published HERMES A_parallel(phi) fit amplitude"
                if observable.startswith("UnpolarizedCosPhi_")
                else "see measurement definition"
            ),
        }
        objects.append(
            _estimate_with_bands(
                yoda,
                prediction,
                path,
                annotations,
                bands.get(observable) if reference_path else None,
            )
        )
        dataset = snapshot["datasets"].get(observable)
        if dataset and dataset["projection"] != "difference":
            goodness[observable] = sidis.projection_goodness_of_fit(
                prediction, dataset
            )
            pull_values = _sidis_pull_bins(prediction, dataset)
            pulls[observable] = pull_values
            if any(value is not None for value in pull_values):
                objects.append(
                    experimental._estimate_from_values(
                        yoda,
                        prediction["edges"],
                        f"/{measurement['analysis']['name']}/Pull_{observable}",
                        pull_values,
                        [
                            0.0 if value is not None else None
                            for value in pull_values
                        ],
                        {
                            "Observable": "(theory-data)/combined uncertainty",
                            "CovarianceNote": (
                                "Displayed pointwise pulls; covariance-aware "
                                "decorrelated pulls are in summary.json"
                            ),
                        },
                    )
                )
        for index, (value, error) in enumerate(
            zip(prediction["values"], prediction["errors"]), start=1
        ):
            rows.append(
                {
                    "observable": observable,
                    "bin": index,
                    "low": prediction["edges"][index - 1],
                    "high": prediction["edges"][index],
                    "value": value,
                    "mc_stat": error,
                }
            )

    output_dir = campaign_dir / "postprocess"
    output_dir.mkdir(parents=True, exist_ok=True)
    experimental._write_yoda_objects(yoda, objects, prediction_path)
    summary = {
        "measurement": measurement["id"],
        "tag": args.tag,
        "hard_process_accuracy": "POWHEG NLO+PS",
        "central_sample_label": measurement["families"]["nominal"]["label"],
        "uncertainties": bands,
        "goodness_of_fit_by_projection": goodness,
        "combined_goodness_of_fit": None,
        "combined_goodness_of_fit_note": (
            "No cross-projection covariance is published; 1D, 2D, 3D, "
            "and charge-difference projections are never combined."
        ),
        "pointwise_pulls": pulls,
        "charge_difference_note": (
            "The APS supplement publishes charge-difference A1 only. "
            "Generator A_parallel charge differences are diagnostic and are "
            "not used in goodness-of-fit calculations."
        ),
        "azimuthal_note": (
            "The optional UnpolarizedCosPhi outputs are helicity-summed yield "
            "moments. They are not overlaid on the published HERMES "
            "A_parallel(phi) cosine-fit amplitudes."
        ),
        "include_diagnostics": include_diagnostics,
        "variations": {
            _variation_id(key): value for key, value in predictions.items()
        },
    }
    experimental.atomic_write_json(output_dir / "summary.json", summary)
    if rows:
        with (output_dir / "central.csv").open(
            "w", encoding="utf-8", newline=""
        ) as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    manifest["postprocess"] = {
        "created_at": experimental.utc_now(),
        "prediction": str(prediction_path.relative_to(campaign_dir)),
        "predictions": [
            {
                "family": "nominal",
                "label": measurement["families"]["nominal"]["label"],
                "path": str(prediction_path.relative_to(campaign_dir)),
            }
        ],
        "summary": "postprocess/summary.json",
        "include_diagnostics": include_diagnostics,
    }
    manifest["updated_at"] = experimental.utc_now()
    manifest["history"].append(
        {"at": experimental.utc_now(), "action": "postprocess"}
    )
    experimental.atomic_write_json(manifest_path, manifest)
    print(f"Wrote normalized SIDIS predictions to {prediction_path}")
    return prediction_path


def postprocess_pp(args: argparse.Namespace, measurement: Mapping[str, Any]) -> Path:
    campaign_dir = _campaign_dir(measurement["id"], args.tag)
    manifest_path = campaign_dir/experimental.MANIFEST_NAME
    if not manifest_path.exists():
        raise CampaignError(f"No prepared campaign at {campaign_dir}")
    manifest = _load_json(manifest_path)
    _assert_manifest_signature_current(manifest, measurement)
    groups = _logical_groups(manifest, campaign_dir)
    if args.dry_run:
        print(json.dumps({"measurement": measurement["id"], "groups": len(groups),
                          "output": str(campaign_dir/"postprocess"/"prediction.yoda")},
                         indent=2))
        return campaign_dir/"postprocess"/"prediction.yoda"
    snapshot = _measurement_snapshot(measurement)
    analysis = str(measurement["analysis"]["name"])
    predictions: dict[tuple[Any, ...], dict[str, Any]] = {}
    raw_samples: dict[tuple[Any, ...], dict[str, Mapping[str, experimental.BinSeries]]] = {}
    raw_statistics: dict[
        tuple[Any, ...], dict[str, Mapping[str, experimental.BinSeries]]
    ] = {}
    for key, helicity_jobs in groups.items():
        family, channel, polarized, unpolarized, scale, mpi = key
        required = set(measurement["families"][family]["helicities"])
        if set(helicity_jobs) != required:
            raise CampaignError(f"Incomplete helicity matrix for {_variation_id(key)}")
        channel_spec = measurement["channels"][channel]
        object_names = (
            channel_spec.get("raw_objects")
            or {"yield": channel_spec["raw_object"]}
        )
        statistic_names = channel_spec.get("statistics_objects", {})
        requested = {
            **{
                f"raw::{observable}": _raw_object_source(
                    measurement, object_spec
                )
                for observable, object_spec in object_names.items()
            },
            **{
                f"stat::{statistic}": _raw_object_source(
                    measurement, object_spec
                )
                for statistic, object_spec in statistic_names.items()
            },
        }
        loaded_by_helicity = {
            helicity: _load_series_many(
                jobs, campaign_dir, requested
            )
            for helicity, jobs in helicity_jobs.items()
        }
        objects: dict[str, Mapping[str, experimental.BinSeries]] = {
            observable: {
                helicity: loaded[f"raw::{observable}"]
                for helicity, loaded in loaded_by_helicity.items()
            }
            for observable in object_names
        }
        if (measurement["process_kind"] == "polarized_pp_jets" and
                measurement["families"][family].get("observable_level") == "hard_parton"):
            for helicity, jobs in helicity_jobs.items():
                provenance = _load_series(
                    jobs, campaign_dir, analysis, "ProvenanceStatus"
                )
                if len(provenance.values) < 2 or provenance.values[1] != 0.0:
                    raise CampaignError(
                        f"Hard-parton provenance failed for "
                        f"{_variation_id(key)}/{helicity}"
                    )
        raw_samples[key] = objects
        statistics: dict[str, Mapping[str, experimental.BinSeries]] = {
            statistic: {
                helicity: loaded[f"stat::{statistic}"]
                for helicity, loaded in loaded_by_helicity.items()
            }
            for statistic in statistic_names
        }
        raw_statistics[key] = statistics
        if required != set(DENOMINATOR):
            continue
        if measurement["postprocessor"] == "star_weak_bosons":
            predictions[key] = _star_prediction(channel, objects["yield"])
        elif measurement["postprocessor"] == "phenix_prompt_photon":
            predictions[key] = _phenix_prediction(objects)
        elif measurement["postprocessor"] == "star_jet_all":
            predictions[key] = _star_jet_prediction(objects)
            _apply_star_display_binning(predictions[key], snapshot)
        elif measurement["postprocessor"] == "mc_poldijets":
            predictions[key] = _mc_poldijets_prediction(
                objects,
                channel_spec.get("helicity_resolved_observables", []),
            )
        elif measurement["postprocessor"] == "mc_poljetshapes":
            predictions[key] = _mc_poljetshapes_prediction(
                objects,
                channel_spec.get("angular_observables", []),
                channel_spec.get("conditional_angular_observables", {}),
                channel_spec.get("support_observables", []),
            )
        else:
            raise CampaignError(f"Unknown postprocessor {measurement['postprocessor']}")

    configured_families = list(manifest["configuration"]["families"])
    physical_families = [
        family_id
        for family_id in configured_families
        if set(measurement["families"][family_id]["helicities"])
        == set(DENOMINATOR)
    ]
    if not physical_families:
        raise CampaignError("No four-helicity physics family was generated")
    central_family = (
        "nominal" if "nominal" in physical_families else physical_families[0]
    )
    nominal_available = central_family == "nominal"
    nominal_mpi = measurement["families"]["nominal"]["mpi"]
    central_mpi = measurement["families"][central_family]["mpi"]
    central_keys = [(central_family, channel, 0, 0, 1.0, central_mpi)
                    for channel in measurement["channels"]]
    missing_central = [key for key in central_keys if key not in predictions]
    if missing_central:
        raise CampaignError(f"Missing central predictions: {missing_central}")
    _mask_spin_asymmetry_support(measurement, predictions, raw_samples)
    bands = (
        aggregate_uncertainties(predictions, measurement)
        if nominal_available
        else {}
    )
    correlated_goodness_of_fit = None
    if measurement["postprocessor"] == "star_jet_all" and nominal_available:
        correlated_goodness_of_fit = _star_correlated_goodness_of_fit(
            predictions[central_keys[0]],
            snapshot,
            measurement,
            require_complete=not bool(
                manifest["configuration"].get("smoke", False)
            ),
        )
    shape_assessment: dict[str, Any] = {}
    comparison_predictions: dict[str, dict[str, Any]] = {}
    sensitivity_ranking: list[dict[str, Any]] = []
    if measurement["postprocessor"] == "mc_poljetshapes":
        (
            shape_assessment,
            comparison_predictions,
            sensitivity_ranking,
        ) = _mc_poljetshapes_assessment(
            measurement, predictions, raw_samples, manifest
        )
    shard_covariance_summary: dict[str, Any] = {}
    shard_covariance_rows: list[dict[str, Any]] = []
    if measurement["postprocessor"] == "mc_poljetshapes":
        (
            shard_covariance_summary,
            shard_covariance_rows,
        ) = _mc_poljetshapes_shard_covariance(
            measurement,
            groups,
            campaign_dir,
            raw_samples,
            predictions,
        )
    yoda = experimental._import_yoda()
    include_diagnostics = bool(
        getattr(args, "include_diagnostics", False)
    )
    objects_by_family: dict[str, list[Any]] = {
        family: [] for family in manifest["configuration"]["families"]
    }
    summary_variations: dict[str, Any] = {}
    central_rows: list[dict[str, Any]] = []
    pulls_summary: dict[str, Any] = {}

    for key, prediction_set in predictions.items():
        variation = _variation_id(key)
        is_nominal_central = (
            key[0] == "nominal" and key[2:5] == (0,0,1.0)
            and key[5] == nominal_mpi
        )
        is_family_central = key[2:5] == (0,0,1.0)
        summary_variations[variation] = prediction_set
        for observable, prediction in prediction_set.items():
            reference_path = _reference_path(measurement, observable, snapshot)
            primary = _primary_reference_observable(
                measurement, snapshot, observable
            )
            mask = (
                _comparison_mask(
                    measurement, snapshot, observable, len(prediction["values"])
                )
                if primary
                else [True] * len(prediction["values"])
            )
            displayed_prediction = star_policy.masked_copy(prediction, mask)
            if not is_family_central:
                # PDF replicas and scale points contribute to named bands,
                # not a forest of individual curves.
                continue
            if not primary and not include_diagnostics:
                continue
            path = (
                reference_path
                if reference_path
                else f"/{analysis}/DIAGNOSTICS/{observable}"
            )
            annotation = {
                "Generator": "HerwigPol full event simulation",
                "HardProcessAccuracy": "LO",
                "PDFInputs": "NNPDF40 NLO + NNPDFpol2.0 NLO",
                "SampleLabel": measurement["families"][key[0]]["label"],
                "Channel": key[1], "PolarizedPDFMember": key[2],
                "UnpolarizedPDFMember": key[3], "HardScaleFactor": key[4],
                "MPI": key[5], "HelicityCombination": "independent PP,PM,MP,MM samples",
                "ShowerSpinCorrelations": measurement["families"][key[0]].get(
                    "shower_spin_correlations", "on"
                ),
                "HardProcessSpin": measurement["families"][key[0]].get(
                    "hard_process_spin", "on"
                ),
                "JetKtMinGeV": manifest["configuration"].get(
                    "jet_kt_min_gev"
                ),
            }
            helicity_match = re.match(
                r"^Sigma(PP|PM|MP|MM)_", observable
            )
            if helicity_match:
                annotation["HelicityCombination"] = helicity_match.group(1)
            annotation["ComparisonRole"] = (
                "primary with excluded bins stored under DIAGNOSTICS"
                if not all(mask) else "primary"
            )
            objects_by_family.setdefault(str(key[0]), []).append(_estimate_with_bands(
                yoda, displayed_prediction, path, annotation,
                _masked_bands(bands[observable], mask)
                if (
                    is_nominal_central
                    and primary
                    and observable in bands
                )
                else None))
            if is_nominal_central and primary and not all(mask):
                objects_by_family["nominal"].append(
                    _estimate_with_bands(
                        yoda,
                        prediction,
                        f"/{analysis}/DIAGNOSTICS/{observable}_full",
                        {
                            **annotation,
                            "ComparisonRole":
                            "diagnostic_only below the configured primary threshold",
                        },
                    )
                )
            if is_nominal_central and primary:
                points = _reference_points(measurement, snapshot, observable)
                if points:
                    pull_values, chi2, count = _pulls(prediction, points, mask)
                    pulls_summary[observable] = {"chi2": chi2, "points": count,
                                                 "values": pull_values,
                                                 "comparison_roles":
                                                 star_policy.comparison_roles(mask),
                                                 "note": "PDF fit overlap prevents interpreting this as an independent PDF validation"}
                    if count:
                        objects_by_family["nominal"].append(
                            experimental._estimate_from_values(
                                yoda,
                                prediction["edges"],
                                f"/{analysis}/Pull_{observable}",
                                pull_values,
                                [
                                    0.0 if value is not None else None
                                    for value in pull_values
                                ],
                                {
                                    "Observable":
                                    "(theory-data)/combined uncertainty"
                                },
                            )
                        )
                for index, (value, error) in enumerate(
                    zip(prediction["values"], prediction["errors"]), start=1
                ):
                    analysis_edges = prediction.get(
                        "analysis_edges", prediction["edges"]
                    )
                    central_rows.append(
                        {
                            "observable": observable,
                            "bin": index,
                            "analysis_low": analysis_edges[index - 1],
                            "analysis_high": analysis_edges[index],
                            "display_low": prediction["edges"][index - 1],
                            "display_high": prediction["edges"][index],
                            "value": value,
                            "mc_stat": error,
                            "comparison_role":
                            star_policy.comparison_roles(mask)[index - 1],
                        }
                    )

    # Add the explicitly requested unpolarized closure when the paper profile
    # has generated it.  It is kept separate from the nominal physics curves.
    for channel in measurement["channels"]:
        if not include_diagnostics:
            break
        physical_key = ("nominal", channel, 0, 0, 1.0, nominal_mpi)
        closure_mpi = measurement["families"].get(
            "unpolarized_closure", {}
        ).get("mpi", nominal_mpi)
        closure_key = (
            "unpolarized_closure", channel, 0, 0, 1.0, closure_mpi
        )
        if closure_key not in raw_samples:
            continue
        physical_objects = raw_samples[physical_key]
        closure_objects = raw_samples[closure_key]
        for object_label in physical_objects:
            unpolarized_series = closure_objects[object_label].get("00")
            if unpolarized_series is None:
                raise CampaignError(f"Missing 00 closure sample for {channel}/{object_label}")
            closure = _unpolarized_closure(physical_objects[object_label],
                                            unpolarized_series)
            objects_by_family.setdefault("unpolarized_closure", []).append(
                _estimate_with_bands(
                yoda, closure, f"/{analysis}/UnpolarizedClosure_{channel}_{object_label}",
                {"Observable": "(sigma_UU-sigma_00)/(sigma_UU+sigma_00)",
                 "SampleLabel": measurement["families"]["unpolarized_closure"]["label"]})
            )

    output_dir = campaign_dir/"postprocess"
    output_dir.mkdir(parents=True, exist_ok=True)
    if shard_covariance_summary:
        experimental.atomic_write_json(
            output_dir / "shard-block-covariance.json",
            {
                "measurement": measurement["id"],
                "tag": args.tag,
                **shard_covariance_summary,
            },
        )
        if shard_covariance_rows:
            with (output_dir / "shard-block-covariance.csv").open(
                "w", encoding="utf-8", newline=""
            ) as stream:
                writer = csv.DictWriter(
                    stream,
                    fieldnames=(
                        "sample", "observable", "bin_i", "bin_j",
                        "covariance", "correlation",
                    ),
                )
                writer.writeheader()
                writer.writerows(shard_covariance_rows)
    prediction_path = output_dir/"prediction.yoda"
    prediction_entries: list[dict[str, str]] = []
    primary_prediction_path: Path | None = None
    for family_id, objects in objects_by_family.items():
        if not objects:
            continue
        destination = (
            prediction_path
            if family_id == "nominal"
            else output_dir / f"prediction-{family_id}.yoda"
        )
        experimental._write_yoda_objects(yoda, objects, destination)
        prediction_entries.append(
            {
                "family": family_id,
                "label": measurement["families"][family_id]["label"],
                "path": str(destination.relative_to(campaign_dir)),
            }
        )
        if family_id == central_family:
            primary_prediction_path = destination
    if primary_prediction_path is None or not experimental._nonempty(
        primary_prediction_path
    ):
        raise CampaignError("Central postprocessing produced no prediction objects")
    comparison_path: Path | None = None
    if comparison_predictions:
        comparison_path = output_dir / "comparison-on-minus-off.yoda"
        comparison_objects = [
            experimental._estimate_from_values(
                yoda,
                prediction["edges"],
                f"/{analysis}/COMPARISON/OnMinusOff_{observable}",
                prediction["values"],
                prediction["errors"],
                {
                    "Observable": _comparison_label(measurement, "difference"),
                    "Uncertainty": "independent samples added in quadrature",
                },
            )
            for observable, prediction in sorted(
                comparison_predictions.items()
            )
        ]
        experimental._write_yoda_objects(
            yoda, comparison_objects, comparison_path
        )
        experimental.atomic_write_json(
            output_dir / "comparison.json",
            {
                "measurement": measurement["id"],
                "tag": args.tag,
                "differences": comparison_predictions,
                "sensitivity_ranking": sensitivity_ranking,
                "uncertainty": ("independent full-spin/LHE-like samples"
                                if "comparison_pair" in measurement else
                                "independent spin-on/off samples"),
            },
        )
        if sensitivity_ranking:
            with (output_dir / "sensitivity-ranking.csv").open(
                "w", encoding="utf-8", newline=""
            ) as stream:
                writer = csv.DictWriter(
                    stream, fieldnames=list(sensitivity_ranking[0])
                )
                writer.writeheader()
                writer.writerows(sensitivity_ranking)
        moment_rows: list[dict[str, Any]] = []
        for observable, prediction in sorted(comparison_predictions.items()):
            if not observable.startswith(
                (
                    "A2UU_", "A2LL_", "B2UU_", "B2LL_",
                    "C2UU_", "C2LL_", "S2UU_", "S2LL_",
                )
            ):
                continue
            for index, (value, error) in enumerate(
                zip(prediction["values"], prediction["errors"])
            ):
                moment_rows.append(
                    {
                        "observable": observable,
                        "bin": index + 1,
                        "low": prediction["edges"][index],
                        "high": prediction["edges"][index + 1],
                        "on_minus_off": value,
                        "independent_sample_error": error,
                        "significance": (
                            None if value is None or error in (None, 0.0)
                            else float(value) / float(error)
                        ),
                    }
                )
        if moment_rows:
            with (output_dir / "moment-differences.csv").open(
                "w", encoding="utf-8", newline=""
            ) as stream:
                writer = csv.DictWriter(
                    stream, fieldnames=list(moment_rows[0])
                )
                writer.writeheader()
                writer.writerows(moment_rows)

    statistics_summary: dict[str, Any] = {}
    for key, statistic_set in raw_statistics.items():
        if key[2:5] != (0, 0, 1.0):
            continue
        variation = _variation_id(key)
        statistics_summary[variation] = {}
        for statistic, helicity_series in statistic_set.items():
            statistics_summary[variation][statistic] = {}
            for helicity, series in helicity_series.items():
                statistics_summary[variation][statistic][helicity] = {
                    "edges": list(series.edges),
                    "values_per_generated_event": list(series.values),
                    "errors_per_generated_event": [
                        math.sqrt(max(0.0, variance))
                        for variance in series.variances
                    ],
                    "projected_entries": [
                        value * int(manifest["configuration"]["lo_events"])
                        for value in series.values
                    ],
                }
    if statistics_summary:
        experimental.atomic_write_json(
            output_dir / "cutflow-accepted-counts.json",
            {
                "measurement": measurement["id"],
                "tag": args.tag,
                "statistics": statistics_summary,
                "definitions": snapshot.get("statistics", {}),
            },
        )
    if shape_assessment:
        experimental.atomic_write_json(
            output_dir / "statistics-projection.json", shape_assessment
        )
    summary = {
        "measurement": measurement["id"], "tag": args.tag,
        "hard_process_accuracy": "LO",
        "central_sample_label": measurement["families"][central_family]["label"],
        "uncertainties": bands, "pulls": pulls_summary,
        "correlated_goodness_of_fit": correlated_goodness_of_fit,
        "global_experimental_uncertainties": measurement.get("global_uncertainties", {}),
        "jet_kt_min_gev": manifest["configuration"].get("jet_kt_min_gev"),
        "include_diagnostics": include_diagnostics,
        "comparison_policy": star_policy.policy_payload(measurement),
        "comparison_policy_sha256": star_policy.policy_sha256(measurement),
        "spin_density_policy": copy.deepcopy(
            measurement.get("spin_density_policy", {})
        ),
        "primary_covariance_points": 0 if correlated_goodness_of_fit is None else correlated_goodness_of_fit.get("points"),
        "variations": summary_variations,
        "cutflow_and_accepted_counts": statistics_summary,
        "shower_spin_comparison": {
            "sensitivity_ranking": sensitivity_ranking,
            "differences": comparison_predictions,
        },
        "statistics_projection": shape_assessment,
        "shard_block_covariance": shard_covariance_summary,
    }
    if "comparison_pair" in measurement:
        summary["comparison_pair"] = list(_comparison_pair(measurement))
        summary["shower_spin_policy"] = manifest["configuration"]["shower_spin_policy"]
        summary["analysis_instances"] = manifest["configuration"].get("analysis_instances", {})
        summary["runtime_provenance"] = manifest.get("runtime", {}).get("provenance", {})
        summary["loaded_library_fingerprints"] = manifest.get("runtime", {}).get("loaded_library_fingerprints", {})
    experimental.atomic_write_json(output_dir/"summary.json", summary)
    if central_rows:
        with (output_dir/"central.csv").open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(central_rows[0]))
            writer.writeheader(); writer.writerows(central_rows)
    manifest["postprocess"] = {
        "created_at": experimental.utc_now(),
        "prediction": str(primary_prediction_path.relative_to(campaign_dir)),
        "predictions": prediction_entries,
        "summary": "postprocess/summary.json",
        "include_diagnostics": include_diagnostics,
    }
    if comparison_path is not None:
        manifest["postprocess"]["comparison"] = str(
            comparison_path.relative_to(campaign_dir)
        )
        manifest["postprocess"]["comparison_summary"] = (
            "postprocess/comparison.json"
        )
        manifest["postprocess"]["sensitivity_ranking"] = (
            "postprocess/sensitivity-ranking.csv"
        )
        manifest["postprocess"]["moment_differences"] = (
            "postprocess/moment-differences.csv"
        )
    if statistics_summary:
        manifest["postprocess"]["cutflow_accepted_counts"] = (
            "postprocess/cutflow-accepted-counts.json"
        )
    if shape_assessment:
        manifest["postprocess"]["statistics_projection"] = (
            "postprocess/statistics-projection.json"
        )
    if shard_covariance_summary:
        manifest["postprocess"]["shard_block_covariance"] = (
            "postprocess/shard-block-covariance.json"
        )
        if shard_covariance_rows:
            manifest["postprocess"]["shard_block_covariance_csv"] = (
                "postprocess/shard-block-covariance.csv"
            )
    manifest["updated_at"] = experimental.utc_now()
    manifest["history"].append({"at": experimental.utc_now(), "action": "postprocess"})
    experimental.atomic_write_json(manifest_path, manifest)
    print(f"Wrote normalized-helicity predictions to {primary_prediction_path}")
    return primary_prediction_path


def _prediction_plot_argument(
    prediction: Path, label: str, family: Mapping[str, Any]
) -> str:
    """Return one rivet-mkhtml input with deterministic family styling."""

    argument = f"{prediction}:Title={label}"
    for key, value in sorted(family.get("plot_options", {}).items()):
        argument += f":{key}={value}"
    return argument


def _configure_plot_text_rendering(
    output: Path, environment: Mapping[str, str]
) -> dict[str, Any]:
    """Use Matplotlib MathText when Rivet's TeX toolchain is incomplete."""

    search_path = str(environment.get("PATH", os.environ.get("PATH", "")))
    required_tools = ("latex", "dvipng")
    missing_tools = [
        tool for tool in required_tools
        if shutil.which(tool, path=search_path) is None
    ]
    if not missing_tools:
        return {"mode": "latex", "missing_tools": []}

    style_path = output / "default.mplstyle"
    if not style_path.is_file():
        raise CampaignError(
            f"Cannot configure MathText fallback; missing {style_path}"
        )
    style = style_path.read_text(encoding="utf-8")
    configured, replacements = re.subn(
        r"(?m)^(\s*text\.usetex\s*:\s*)True(\s*(?:#.*)?)$",
        r"\g<1>False\g<2>",
        style,
    )
    if replacements != 1:
        if not re.search(
            r"(?m)^\s*text\.usetex\s*:\s*False\s*(?:#.*)?$", style
        ):
            raise CampaignError(
                f"Could not disable text.usetex in {style_path}"
            )
        configured = style
    if configured != style:
        experimental.atomic_write_text(style_path, configured)
    return {"mode": "mathtext", "missing_tools": missing_tools}


def _label_companion_jet_spectrum(script: Path) -> None:
    """Companion spectra live in the primary output namespace; label their axes."""
    match = re.search(r"_jet([1-4])_(pt|eta)\.py$", script.name)
    if not match:
        return
    variable = r"$p_T$ [GeV]" if match.group(2) == "pt" else r"$\eta$"
    label = f"Jet {match.group(1)} {variable}"
    source = script.read_text(encoding="utf-8")
    replacement = f"ax.set_xlabel({label!r})"
    source = source.replace("ax.set_xlabel(ax_xLabel)", replacement)
    experimental.atomic_write_text(script, source)


def _write_mc_poljetshapes_ratio_yoda(
    measurement: Mapping[str, Any],
    ratios: Mapping[str, Mapping[str, Any]],
    destination: Path,
) -> Path:
    """Write plot-only spin-on/spin-off ratios under the original paths."""

    yoda = experimental._import_yoda()
    analysis = str(measurement["analysis"]["name"])
    objects = [
        experimental._estimate_from_values(
            yoda,
            prediction["edges"],
            f"/{analysis}/{observable}",
            prediction["values"],
            prediction["errors"],
            {
                "Observable": _comparison_label(measurement, "ratio"),
                "Uncertainty": (
                    "independent numerator and denominator MC errors "
                    "propagated in quadrature"
                ),
                "Numerator": measurement["families"]["nominal"]["label"],
                "Denominator": measurement["families"][_comparison_pair(measurement)[1]]["label"],
            },
        )
        for observable, prediction in sorted(ratios.items())
    ]
    if not objects:
        raise CampaignError("No MC_POLJETSHAPES ratios were available to write")
    experimental._write_yoda_objects(yoda, objects, destination)
    experimental.atomic_write_json(
        destination.with_suffix(".json"),
        {
            "measurement": measurement["id"],
            "definition": _comparison_label(measurement, "ratio"),
            "uncertainty": (
                "independent numerator and denominator MC errors propagated "
                "in quadrature"
            ),
            "excluded_signed_or_zero_crossing_prefixes": [
                "ALL_", "DeltaSigmaLL_", "A2UU_", "A2LL_", "B2UU_", "B2LL_",
                "C2UU_", "C2LL_", "S2UU_", "S2LL_",
            ],
            "ratios": ratios,
        },
    )
    return destination


def _configure_mc_poljetshapes_ratio_script(
    script: Path, label: str = "Shower spin on / off"
) -> None:
    """Give a ratio-only Rivet script a linear scale and unity reference."""

    source = script.read_text(encoding="utf-8")
    if "# MC_POLJETSHAPES spin-on/spin-off ratio presentation" in source:
        return
    source, replacements = re.subn(
        r"(?m)^ax_yScale\s*=\s*['\"](?:linear|log)['\"]$",
        "ax_yScale = 'linear'",
        source,
        count=1,
    )
    if replacements != 1:
        raise CampaignError(f"Could not set linear ratio scale in {script}")
    marker = "\n\nlegend_handles = dict() # keep track of handles for the legend"
    if marker not in source:
        raise CampaignError(
            f"Could not locate generated-data marker in ratio plot {script}"
        )
    presentation = r'''

# MC_POLJETSHAPES spin-on/spin-off ratio presentation
ax_yLabel = 'Shower spin on / off'
_ratio_extent = [1.0]
for _ratio_label, _ratio_values in dataf.get('yvals', {}).items():
    _ratio_errors = dataf.get('yerrs', {}).get(_ratio_label, ([], []))
    _ratio_down = _ratio_errors[0] if len(_ratio_errors) > 0 else []
    _ratio_up = _ratio_errors[1] if len(_ratio_errors) > 1 else []
    for _ratio_index, _ratio_value in enumerate(_ratio_values):
        if not np.isfinite(_ratio_value):
            continue
        _ratio_extent.append(float(_ratio_value))
        if _ratio_index < len(_ratio_down) and np.isfinite(_ratio_down[_ratio_index]):
            _ratio_extent.append(float(_ratio_value - _ratio_down[_ratio_index]))
        if _ratio_index < len(_ratio_up) and np.isfinite(_ratio_up[_ratio_index]):
            _ratio_extent.append(float(_ratio_value + _ratio_up[_ratio_index]))
_ratio_low = min(_ratio_extent)
_ratio_high = max(_ratio_extent)
_ratio_span = max(_ratio_high - _ratio_low, 0.08)
yLims = (_ratio_low - 0.12*_ratio_span, _ratio_high + 0.12*_ratio_span)
ax.axhline(1.0, color='#666666', linestyle='--', linewidth=1.0, zorder=1)
'''
    presentation = presentation.replace("'Shower spin on / off'", repr(label))
    source = source.replace(marker, presentation + marker, 1)
    experimental.atomic_write_text(script, source)


def _write_mc_poljetshapes_focus_index(
    output: Path, measurement: Mapping[str, Any]
) -> Path:
    """Write and link a curated physics-first view of the complete gallery."""

    analysis = str(measurement["analysis"]["name"])
    focus_dir = output / "focus"
    focus_dir.mkdir(parents=True, exist_ok=True)

    def asset_panel(base: Path, panel_title: str) -> str:
        png = output / base.with_suffix(".png")
        pdf = output / base.with_suffix(".pdf")
        if not experimental._nonempty(png):
            return (
                '<div class="panel missing"><h4>'
                + html.escape(panel_title)
                + "</h4><p>Not rendered: the result was empty or fully masked.</p></div>"
            )
        png_link = "../" + base.with_suffix(".png").as_posix()
        links = [f'<a href="{html.escape(png_link)}">PNG</a>']
        if experimental._nonempty(pdf):
            pdf_link = "../" + base.with_suffix(".pdf").as_posix()
            links.append(f'<a href="{html.escape(pdf_link)}">PDF</a>')
        return (
            '<div class="panel"><h4>'
            + html.escape(panel_title)
            + "</h4>"
            + f'<a href="{html.escape(png_link)}"><img '
            + f'src="{html.escape(png_link)}" alt="{html.escape(panel_title)}"></a>'
            + f"<p>{' &middot; '.join(links)}</p></div>"
        )

    sections: list[str] = []
    rendered_cards = 0
    for section in MC_POLJETSHAPES_FOCUS_SECTIONS:
        cards: list[str] = []
        for stem, label, companion in section["plots"]:
            main_base = Path(analysis) / str(stem)
            if not experimental._nonempty(output / main_base.with_suffix(".png")):
                continue
            if companion == "ratio":
                auxiliary_base = Path("ratios") / analysis / str(stem)
                auxiliary_title = (_comparison_label(measurement, "ratio")
                                   if "comparison_pair" in measurement else "Spin on / spin off")
            elif companion == "difference":
                auxiliary_base = (
                    Path(analysis) / "COMPARISON" / f"OnMinusOff_{stem}"
                )
                auxiliary_title = (_comparison_label(measurement, "difference")
                                   if "comparison_pair" in measurement else "Spin on minus spin off")
            else:
                raise CampaignError(
                    f"Unknown focused-gallery companion {companion!r}"
                )
            cards.append(
                '<article class="card"><h3>'
                + html.escape(str(label))
                + '</h3><div class="panels">'
                + asset_panel(main_base, "Red/blue family overlay")
                + asset_panel(auxiliary_base, auxiliary_title)
                + "</div></article>"
            )
        if not cards:
            continue
        rendered_cards += len(cards)
        sections.append(
            "<section><h2>"
            + html.escape(str(section["title"]))
            + "</h2><p>"
            + html.escape(str(section["why"]))
            + "</p>"
            + "\n".join(cards)
            + "</section>"
        )
    if not sections:
        raise CampaignError(
            "No recommended MC_POLJETSHAPES plots were rendered for the focus page"
        )

    covariance_panel = ""
    covariance_path = (
        output.parent.parent / "postprocess" / "shard-block-covariance.json"
    )
    if covariance_path.is_file():
        covariance = _load_json(covariance_path)
        rows: list[str] = []
        for observable, result in sorted(
            covariance.get("spin_on_minus_off", {}).items()
        ):
            chi2 = result.get("chi2")
            rank = int(result.get("rank", 0))
            rows.append(
                "<tr><td>"
                + html.escape(str(observable))
                + "</td><td>"
                + ("&mdash;" if chi2 is None else f"{float(chi2):.3g}")
                + "</td><td>"
                + str(rank)
                + "</td><td>"
                + (
                    "&mdash;"
                    if chi2 is None or rank == 0
                    else f"{float(chi2)/rank:.3g}"
                )
                + "</td></tr>"
            )
        family_status = ", ".join(
            f"{name}: {details.get('status', 'unknown')} "
            f"({details.get('blocks', 0)} blocks)"
            for name, details in sorted(covariance.get("families", {}).items())
        )
        covariance_panel = (
            "<section><h2>Shard-block covariance audit</h2><p>"
            "Delete-one common-shard jackknife covariance is used for the "
            "resolved-radiation spectra, rate scans, tails, and conditional "
            "moments. Spin-on and spin-off matrices are added as independent "
            "samples. "
            + html.escape(family_status)
            + ".</p>"
            + (
                '<table><thead><tr><th>Observable</th><th>chi2</th>'
                '<th>rank</th><th>chi2/rank</th></tr></thead><tbody>'
                + "".join(rows)
                + "</tbody></table>"
                if rows
                else "<p>The covariance product is present, but a complete "
                     "two-family difference matrix was not available.</p>"
            )
            + '<p><a href="../../../postprocess/shard-block-covariance.json">'
              "Full covariance JSON</a> &middot; "
              '<a href="../../../postprocess/shard-block-covariance.csv">'
              "flat covariance CSV</a></p></section>"
        )

    title = f"{measurement['title']} — recommended comparisons"
    document = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{html.escape(title)}</title>
  <style>
    body {{ font-family: sans-serif; margin: 2rem auto; max-width: 1450px; padding: 0 1rem; }}
    .notice {{ background: #f3f6f8; border-left: 5px solid #0077BB; padding: 0.8rem 1rem; }}
    section {{ border-top: 2px solid #bbb; margin-top: 2.5rem; padding-top: 1rem; }}
    .card {{ border-top: 1px solid #ddd; margin-top: 1.5rem; padding-top: 0.5rem; }}
    .panels {{ display: grid; gap: 1rem; grid-template-columns: repeat(2, minmax(0, 1fr)); }}
    .panel {{ min-width: 0; }}
    .panel h4 {{ margin-bottom: 0.4rem; }}
    .panel img {{ height: auto; max-width: 100%; }}
    .missing {{ background: #fafafa; border: 1px dashed #bbb; padding: 1rem; }}
    table {{ border-collapse: collapse; width: 100%; }}
    th, td {{ border: 1px solid #ccc; padding: 0.45rem; text-align: left; }}
    th {{ background: #f3f6f8; }}
    @media (max-width: 850px) {{ .panels {{ grid-template-columns: 1fr; }} }}
  </style>
</head>
<body>
  <p><a href="../index.html">Back to the complete gallery</a></p>
  <h1>{html.escape(title)}</h1>
  <div class="notice">
    <p><strong>Reading convention:</strong> red is shower spin on and blue is
    shower spin off. Both samples retain polarized beams and the polarized hard
    process. Ratio panels are shower spin on divided by shower spin off and
    propagate both independent Monte Carlo errors.</p>
    <p>Ratios are used only for positive cross sections, normalized shapes, and
    rate observables. Signed Delta-sigma LL, A_LL, and angular moments use
    on-minus-off differences instead. Isolated bins should not be interpreted
    without a coherent angular pattern and the null tests.</p>
  </div>
  <p>{rendered_cards} recommended plot pairs are shown below.</p>
  {covariance_panel}
  {''.join(sections)}
</body>
</html>
"""
    if "comparison_pair" in measurement:
        document = document.replace("red is shower spin on and blue is\n    shower spin off", "red is full-spin showering and blue is LHE-like showering")
        document = document.replace("shower spin on divided by shower spin off",
                                    html.escape(_comparison_label(measurement, "ratio")))
        document = document.replace("spin on minus spin off",
                                    html.escape(_comparison_label(measurement, "difference")))
        document = document.replace("Both samples retain polarized beams and the polarized hard\n    process.",
            "Both samples retain the same polarized hard process and ordinary shower-generated spin correlations. "
            "The blue sample discards hard spin input and polarized backward-ISR conditioning; "
            "LHE equivalence requires separately documented closure validation.")
    index = focus_dir / "index.html"
    experimental.atomic_write_text(index, document)

    root_index = output / "index.html"
    root = root_index.read_text(encoding="utf-8")
    if 'href="focus/index.html"' not in root:
        marker = "</h1>"
        link = (
            '</h1><p class="focus-link"><strong><a href="focus/index.html">'
            "Open the recommended shower-spin comparison panel"
            "</a></strong> — curated moments, shapes, null tests, rates, and "
            "spin-on/spin-off ratios.</p>"
        )
        if marker not in root:
            raise CampaignError(f"Could not link focused gallery from {root_index}")
        root = root.replace(marker, link, 1)
        experimental.atomic_write_text(root_index, root)
    return index


def plot_pp(args: argparse.Namespace, measurement: Mapping[str, Any]) -> Path:
    campaign_dir = _campaign_dir(measurement["id"], args.tag)
    manifest_path = campaign_dir/experimental.MANIFEST_NAME
    if not manifest_path.exists():
        raise CampaignError(f"No prepared campaign at {campaign_dir}")
    manifest = _load_json(manifest_path)
    plot_metadata_refresh = _assert_manifest_signature_current(
        manifest,
        measurement,
        allow_plot_metadata_refresh=bool(
            getattr(args, "allow_plot_metadata_refresh", False)
        ),
    )
    postprocess = manifest.get("postprocess", {})
    entries = postprocess.get("predictions")
    if not isinstance(entries, list) or not entries:
        entries = [
            {
                "family": "nominal",
                "label": measurement["families"]["nominal"]["label"],
                "path": postprocess.get(
                    "prediction", "postprocess/prediction.yoda"
                ),
            }
        ]
    external_nominal_prediction: Path | None = None
    nominal_prediction_option = getattr(args, "nominal_prediction", None)
    if nominal_prediction_option is not None:
        if not bool(getattr(args, "plot_comparisons", False)):
            raise CampaignError(
                "--nominal-prediction requires --plot-comparisons"
            )
        if any(str(entry.get("family")) == "nominal" for entry in entries):
            raise CampaignError(
                "--nominal-prediction is only valid for a campaign without "
                "its own nominal family"
            )
        external_nominal_prediction = Path(nominal_prediction_option).expanduser()
        if not external_nominal_prediction.is_absolute():
            external_nominal_prediction = (
                Path.cwd() / external_nominal_prediction
            ).resolve()
        else:
            external_nominal_prediction = external_nominal_prediction.resolve()
        entries = [
            {
                "family": "nominal",
                "label": measurement["families"]["nominal"]["label"],
                "path": str(external_nominal_prediction),
            },
            *entries,
        ]
    if not bool(getattr(args, "plot_comparisons", False)):
        entries = [
            entry for entry in entries
            if str(entry.get("family")) == "nominal"
        ] or entries[:1]
    predictions: list[tuple[Path, str, str]] = []
    nominal_prediction: Path | None = None
    for entry in entries:
        prediction = campaign_dir / str(entry["path"])
        if not experimental._nonempty(prediction):
            raise CampaignError(f"Missing postprocessed YODA {prediction}")
        family_id = str(entry.get("family", "nominal"))
        predictions.append(
            (
                prediction,
                str(entry.get("label", family_id)),
                family_id,
            )
        )
        if family_id == "nominal":
            nominal_prediction = prediction
    if nominal_prediction is None:
        nominal_prediction = predictions[0][0]
    campaign_summary_path = campaign_dir/"postprocess"/"summary.json"
    campaign_summary = (
        _load_json(campaign_summary_path)
        if campaign_summary_path.is_file() else {}
    )
    ratio_predictions: dict[str, dict[str, Any]] = {}
    if (
        measurement["postprocessor"] == "mc_poljetshapes"
        and {family_id for _, _, family_id in predictions}
        >= {"nominal", _comparison_pair(measurement)[1]}
    ):
        ratio_predictions = _mc_poljetshapes_plot_ratios(
            measurement, campaign_summary
        )
        if not ratio_predictions:
            raise CampaignError(
                "MC_POLJETSHAPES comparison plotting requires both central "
                "families in postprocess/summary.json"
            )
    runtime = manifest.get("runtime") or _runtime(measurement)
    output = campaign_dir/"plots"/"html"
    ratio_output = output/"ratios"
    ratio_yoda_path = campaign_dir/"plots"/"spin-on-over-off.yoda"
    # Generate Rivet's plotting scripts without starting its multiprocessing
    # manager.  The latter requires a local IPC socket and is unavailable in
    # some batch/sandbox environments; executing the generated scripts
    # sequentially gives identical plots and a deterministic failure point.
    safe_wrapper = DISPOL_ROOT/"scripts"/"rivet_mkhtml_safe.py"
    command = [sys.executable, str(safe_wrapper),
               runtime["tools"]["rivet-mkhtml"], "--dry-run",
               "--offline", "--deviation", "--pwd", "-o", str(output)]
    include_diagnostics = bool(
        getattr(args, "include_diagnostics", False)
    )
    if not include_diagnostics:
        command.extend(["-M", r".*/DIAGNOSTICS/.*"])
        if measurement["postprocessor"] == "star_jet_all":
            snapshot = _measurement_snapshot(measurement)
            for dataset in snapshot["datasets"].values():
                if dataset.get("alternate_projection"):
                    command.extend(
                        ["-M", re.escape(str(dataset["rivet_path"]))]
                    )
    for prediction, label, family_id in predictions:
        command.append(
            _prediction_plot_argument(
                prediction, label, measurement["families"][family_id]
            )
        )
    comparison_plot_path: Path | None = None
    if measurement["postprocessor"] == "mc_poljetshapes":
        comparison_entry = postprocess.get("comparison")
        if comparison_entry:
            comparison_plot_path = campaign_dir / str(comparison_entry)
            if not experimental._nonempty(comparison_plot_path):
                raise CampaignError(
                    f"Missing shower-spin comparison YODA {comparison_plot_path}"
                )
            command.append(
                f"{comparison_plot_path}:Title={_comparison_label(measurement, 'difference')}:"
                "LineColor=#882255"
            )
    ratio_command: list[str] | None = None
    if ratio_predictions:
        ratio_command = [
            sys.executable, str(safe_wrapper),
            runtime["tools"]["rivet-mkhtml"], "--dry-run", "--offline",
            "--no-ratio", "--pwd", "-o", str(ratio_output),
            f"{ratio_yoda_path}:Title={_comparison_label(measurement, 'ratio')}:LineColor=#CC3311",
        ]
    if args.dry_run:
        print(" ".join(command))
        if ratio_command is not None:
            print(" ".join(ratio_command))
        return output
    # rivet-mkhtml does not remove scripts for objects that disappeared after
    # a new postprocessing pass (for example an all-empty smoke-test pull).
    # Recreate only this generated plot subtree so stale scripts cannot make a
    # later, otherwise valid plot command fail.
    if output.exists():
        shutil.rmtree(output)
    output.mkdir(parents=True, exist_ok=True)
    environment = experimental.analysis_environment(measurement, campaign_dir, runtime)
    mpl_cache = campaign_dir/"work"/"matplotlib-cache"
    mpl_cache.mkdir(parents=True, exist_ok=True)
    environment["MPLCONFIGDIR"] = str(mpl_cache)
    environment["DISPOL_FORCE_NO_ERROR_BANDS"] = "1"
    if ratio_predictions:
        _write_mc_poljetshapes_ratio_yoda(
            measurement, ratio_predictions, ratio_yoda_path
        )
    experimental._run_logged(command, campaign_dir, environment,
                             campaign_dir/"logs"/"rivet-mkhtml-generate.log")
    if ratio_command is not None:
        experimental._run_logged(
            ratio_command, campaign_dir, environment,
            campaign_dir/"logs"/"rivet-mkhtml-ratios-generate.log",
        )
    plot_scripts = sorted(
        path for path in output.rglob("*.py")
        if not path.name.endswith("__data.py")
    )
    if not plot_scripts:
        raise CampaignError(f"rivet-mkhtml generated no plot scripts below {output}")
    text_rendering = _configure_plot_text_rendering(output, environment)
    ratio_text_rendering: dict[str, Any] | None = None
    if ratio_command is not None:
        ratio_text_rendering = _configure_plot_text_rendering(
            ratio_output, environment
        )
    script_log = campaign_dir/"logs"/"rivet-plot-scripts.log"
    rendered_scripts: list[Path] = []
    snapshot = _measurement_snapshot(measurement)
    summary = campaign_summary
    external_nominal_summary: Path | None = None
    if external_nominal_prediction is not None:
        candidate = external_nominal_prediction.parent / "summary.json"
        if candidate.is_file():
            external_nominal_summary = candidate
            summary = _load_json(candidate)
    with script_log.open("w", encoding="utf-8") as log:
        log.write(
            "text_rendering: "
            + json.dumps(text_rendering, sort_keys=True)
            + "\n"
        )
        for script in plot_scripts:
            log.write(f"script: {script}\n")
            log.flush()
            if "comparison_pair" in measurement:
                _label_companion_jet_spectrum(script)
            is_spin_ratio = ratio_output in script.parents
            if is_spin_ratio:
                _configure_mc_poljetshapes_ratio_script(
                    script, _comparison_label(measurement, "ratio")
                )
            else:
                experimental.add_theory_uncertainty_overlay(
                    script,
                    _pp_theory_uncertainty_bands(
                        measurement, snapshot, summary, script.stem
                    ),
                    nominal_prediction.name,
                )
                experimental.add_experimental_error_overlay(
                    script,
                    _pp_reference_overlay_points(
                        measurement, snapshot, script.stem
                    ),
                    show_statistical=bool(
                        getattr(args, "plot_data_components", False)
                    ),
                )
            if not experimental.plot_script_has_finite_y(script):
                log.write(
                    "skipped: non-renderable finite-data/axis-limit state\n"
                )
                continue
            completed = subprocess.run(
                [sys.executable, str(script)], cwd=script.parent,
                env=environment, stdout=log, stderr=subprocess.STDOUT,
            )
            if completed.returncode != 0:
                raise CampaignError(
                    f"Generated Rivet plot script failed with status "
                    f"{completed.returncode}: {script}; see {script_log}"
                )
            rendered_scripts.append(script)
    index = experimental.write_plot_indexes(output, measurement, rendered_scripts)
    if not index.is_file() or not any(output.rglob("*.png")):
        raise CampaignError(f"Rivet plotting produced no complete HTML below {output}")
    focus_index: Path | None = None
    if ratio_predictions:
        focus_index = _write_mc_poljetshapes_focus_index(output, measurement)
    manifest["plots"] = {
        "created_at": experimental.utc_now(),
        "index": str(index.relative_to(campaign_dir)),
        "plot_metadata_sha256": experimental.sha256_file(
            DISPOL_ROOT / str(measurement["analysis"]["plot"])
        ),
        "text_rendering": text_rendering,
    }
    if ratio_text_rendering is not None:
        manifest["plots"]["ratio_text_rendering"] = ratio_text_rendering
    if ratio_predictions:
        manifest["plots"]["spin_on_over_off"] = {
            "path": str(ratio_yoda_path.relative_to(campaign_dir)),
            "json": str(
                ratio_yoda_path.with_suffix(".json").relative_to(campaign_dir)
            ),
            "sha256": experimental.sha256_file(ratio_yoda_path),
            "observable_count": len(ratio_predictions),
            "uncertainty": "independent numerator and denominator MC errors",
        }
    if focus_index is not None:
        manifest["plots"]["focus_index"] = str(
            focus_index.relative_to(campaign_dir)
        )
    if comparison_plot_path is not None:
        manifest["plots"]["comparison_prediction"] = {
            "path": str(comparison_plot_path.relative_to(campaign_dir)),
            "sha256": experimental.sha256_file(comparison_plot_path),
        }
    if external_nominal_prediction is not None:
        manifest["plots"]["external_nominal_prediction"] = {
            "path": str(external_nominal_prediction),
            "sha256": experimental.sha256_file(external_nominal_prediction),
        }
        if external_nominal_summary is not None:
            manifest["plots"]["external_nominal_summary"] = {
                "path": str(external_nominal_summary),
                "sha256": experimental.sha256_file(external_nominal_summary),
            }
    if plot_metadata_refresh is not None:
        manifest["plots"]["presentation_only_refresh"] = plot_metadata_refresh
    manifest["updated_at"] = experimental.utc_now()
    manifest["history"].append(
        {
            "at": experimental.utc_now(),
            "action": (
                "plot-metadata-refresh"
                if plot_metadata_refresh is not None else "plot"
            ),
        }
    )
    experimental.atomic_write_json(manifest_path, manifest)
    print(f"Wrote Rivet HTML to {index}")
    return output


def _legacy_arguments(args: argparse.Namespace) -> list[str]:
    command = [args.command]
    if args.command != "list":
        command.extend(["--measurement", args.measurement])
    if args.command in {"prepare", "campaign", "full"}:
        command.extend(["--tag", args.tag])
        for option, value in (("--jobs", args.jobs), ("--shards", args.shards),
                              ("--seed-base", args.seed_base),
                              ("--posnlo-events", args.posnlo_events),
                              ("--negnlo-events", args.negnlo_events),
                              ("--lo-events", args.lo_events),
                              ("--progress-interval", args.progress_interval),
                              ("--max-listed", args.max_listed)):
            if value is not None:
                command.extend([option, str(value)])
        if args.smoke: command.append("--smoke")
        if args.dry_run: command.append("--dry-run")
        if args.recover_failed: command.append("--recover-failed")
        command.extend(["--profile", args.profile])
        for option, value in (
            ("--polarized-pdf-members", args.polarized_pdf_members),
            ("--unpolarized-pdf-members", args.unpolarized_pdf_members),
            ("--scales", args.scales),
        ):
            if value is not None:
                command.extend([option, str(value)])
        if args.comparisons:
            command.append("--comparisons")
        if args.command == "full":
            if args.plot_comparisons:
                command.append("--plot-comparisons")
            if args.plot_data_components:
                command.append("--plot-data-components")
    elif args.command in {"postprocess", "plot"}:
        command.extend(["--tag", args.tag])
        if args.dry_run: command.append("--dry-run")
        if args.command == "plot":
            if args.allow_plot_metadata_refresh:
                command.append("--allow-plot-metadata-refresh")
            if args.plot_comparisons:
                command.append("--plot-comparisons")
            if args.plot_data_components:
                command.append("--plot-data-components")
    return command


def _add_measurement(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--measurement", required=True)


def _add_tag(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--tag", required=True)


def _add_campaign_options(parser: argparse.ArgumentParser) -> None:
    _add_measurement(parser); _add_tag(parser)
    parser.add_argument("--jobs", type=int)
    parser.add_argument("--shards", type=int)
    parser.add_argument("--seed-base", type=int)
    parser.add_argument("--posnlo-events", type=int,
                        help="DIS POSNLO events per helicity")
    parser.add_argument("--negnlo-events", type=int,
                        help="DIS NEGNLO events per helicity")
    parser.add_argument("--lo-events", type=int,
                        help="LO events per physical RHIC helicity or DIS comparison")
    parser.add_argument("--profile", choices=("central", "paper"), default="central")
    parser.add_argument(
        "--families",
        help=(
            "nominal, all, or a comma-separated list of physics-family IDs; "
            "the default is nominal even with --profile paper"
        ),
    )
    parser.add_argument("--comparisons", action="store_true",
                        help=(
                            "Backward-compatible alias for --profile paper "
                            "--families all"
                        ))
    parser.add_argument("--polarized-pdf-members",
                        help="central, all, or comma-separated member numbers")
    parser.add_argument("--unpolarized-pdf-members",
                        help="central, all, or comma-separated member numbers")
    parser.add_argument("--scales", help="central, all, or comma-separated 0.5,1,2")
    parser.add_argument(
        "--jet-kt-min-gev",
        type=float,
        help=(
            "polarized-pp jet generator cut in GeV; only the pinned 3, 4, "
            "and 5 GeV scan points are accepted"
        ),
    )
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--recover-failed", "--rerun-failed-random-seed",
                        dest="recover_failed", action="store_true")
    parser.add_argument("--progress-interval", type=float, default=5.0)
    parser.add_argument("--max-listed", type=int, default=12)


def _add_plot_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--plot-comparisons",
        action="store_true",
        help=(
            "Overlay separately postprocessed comparison-family predictions; "
            "the nominal full-spin prediction is the default"
        ),
    )
    parser.add_argument(
        "--plot-data-components",
        action="store_true",
        help=(
            "Overlay statistical-only experimental bars inside the total "
            "data uncertainties; the default shows one total error bar"
        ),
    )
    parser.add_argument(
        "--nominal-prediction",
        type=Path,
        help=(
            "Existing nominal prediction.yoda to prepend when plotting a "
            "comparison-family-only campaign; its path and SHA-256 are "
            "recorded in the plot manifest"
        ),
    )
    parser.add_argument(
        "--include-diagnostics",
        action="store_true",
        help=(
            "Include alternate projections, closure tests, azimuthal moments, "
            "and other diagnostic-only observables"
        ),
    )


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("list")
    fetch = commands.add_parser("fetch-data"); _add_measurement(fetch)
    fetch.add_argument(
        "--source-file",
        type=Path,
        help="Local copy of the pinned APS HERMES supplemental ZIP",
    )
    for name in ("prepare", "campaign"):
        _add_campaign_options(commands.add_parser(name))
    full = commands.add_parser("full")
    _add_campaign_options(full)
    _add_plot_options(full)
    for name in ("postprocess",):
        child = commands.add_parser(name); _add_measurement(child); _add_tag(child)
        child.add_argument("--dry-run", action="store_true")
        child.add_argument("--include-diagnostics", action="store_true")
    plot = commands.add_parser("plot")
    _add_measurement(plot); _add_tag(plot)
    plot.add_argument("--dry-run", action="store_true")
    plot.add_argument(
        "--allow-plot-metadata-refresh",
        action="store_true",
        help=(
            "Replot a complete campaign only when restoring its historical "
            "Rivet .plot file reproduces the immutable generation signature"
        ),
    )
    _add_plot_options(plot)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = make_parser().parse_args(argv)
    args._tracker_started_at = time.time()
    try:
        registry = discover_all()
        if args.command == "list":
            for identifier, measurement in registry.items():
                print(f"{identifier}\t{measurement['process_kind']}\t{measurement.get('title','')}")
            return 0
        if args.measurement not in registry:
            raise CampaignError(f"Unknown measurement {args.measurement!r}")
        measurement = registry[args.measurement]
        if measurement["process_kind"] == "fixed_target_dis":
            if getattr(args, "jet_kt_min_gev", None) is not None:
                raise CampaignError(
                    "--jet-kt-min-gev is only valid for polarized pp jet "
                    "measurements"
                )
            if getattr(args, "nominal_prediction", None) is not None:
                raise CampaignError(
                    "--nominal-prediction is only valid for polarized pp "
                    "campaigns"
                )
            return experimental.main(_legacy_arguments(args))
        if getattr(args, "comparisons", False):
            args.profile = "paper"
            if getattr(args, "families", None) is None:
                args.families = "all"
        if args.command == "fetch-data":
            if measurement["reference"].get("kind") == "internal_observable_definition":
                snapshot = _measurement_snapshot(measurement)
                print(
                    f"Validated internal definition with "
                    f"{len(snapshot['observables'])} observables; "
                    "no external data or reference YODA is required"
                )
            else:
                outputs = fetch_and_validate(
                    args.measurement,
                    source_file=getattr(args, "source_file", None),
                )
                print(
                    f"Checksum- and schema-validated {len(outputs)} source(s); "
                    "refreshed reference YODA"
                )
        elif args.command == "prepare":
            prepare_pp(args, measurement)
        elif args.command == "campaign":
            run_pp(args, measurement)
        elif args.command == "postprocess":
            if measurement["postprocessor"] in {
                "compass_sidis_pt_slope", "sidis_azimuthal_moment",
            }:
                postprocess_sidis_diagnostic(args, measurement)
            elif measurement["postprocessor"] in {
                "compass_sidis_a1", "compass_sidis_multiplicity",
                "hermes_sidis_multiplicity", "compass_sidis_charge_ratio",
            }:
                postprocess_compass_sidis(args, measurement)
            elif measurement["process_kind"] == "polarized_sidis":
                postprocess_sidis(args, measurement)
            else:
                postprocess_pp(args, measurement)
        elif args.command == "plot":
            plot_pp(args, measurement)
        elif args.command == "full":
            prepare_pp(args, measurement)
            if not args.dry_run:
                run_pp(args, measurement)
                if measurement["postprocessor"] in {
                    "compass_sidis_pt_slope", "sidis_azimuthal_moment",
                }:
                    postprocess_sidis_diagnostic(args, measurement)
                elif measurement["postprocessor"] in {
                    "compass_sidis_a1", "compass_sidis_multiplicity",
                    "hermes_sidis_multiplicity", "compass_sidis_charge_ratio",
                }:
                    postprocess_compass_sidis(args, measurement)
                elif measurement["process_kind"] == "polarized_sidis":
                    postprocess_sidis(args, measurement)
                else:
                    postprocess_pp(args, measurement)
                plot_pp(args, measurement)
        else:
            raise CampaignError(f"Unsupported command {args.command}")
    except (CampaignError, ReferenceDataError, experimental.CampaignError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
