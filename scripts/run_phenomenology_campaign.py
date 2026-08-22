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
import json
import math
import os
import re
import shutil
import subprocess
import sys
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import run_experimental_campaign as experimental
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


def _validate_descriptor(measurement: Mapping[str, Any], path: Path) -> None:
    required = {"id", "schema_version", "process_kind", "analysis", "reference",
                "cards", "channels", "families", "campaign", "postprocessor",
                "scales", "pdf_ensembles"}
    missing = sorted(required-set(measurement))
    if missing:
        raise CampaignError(f"{path} is missing: {', '.join(missing)}")
    process_kind = str(measurement["process_kind"])
    if int(measurement["schema_version"]) not in {3, 4} or process_kind not in {
        "polarized_pp", "polarized_pp_jets", "polarized_sidis"
    }:
        raise CampaignError(f"Unsupported schema/process kind in {path}")
    if path.stem != measurement["id"]:
        raise CampaignError(f"Registry filename must match id in {path}")
    if set(measurement["cards"]["helicities"]) != {"PP", "PM", "MP", "MM"}:
        raise CampaignError(f"{measurement['id']} must define four physical helicities")
    contributions = set(measurement["cards"]["contributions"])
    if process_kind in {"polarized_pp", "polarized_pp_jets"} and contributions != {"LO"}:
        raise CampaignError(f"{measurement['id']} RHIC hard processes must be labelled LO")
    if process_kind == "polarized_sidis":
        if contributions != {"POSNLO", "NEGNLO"}:
            raise CampaignError(f"{measurement['id']} SIDIS requires POSNLO and NEGNLO")
        if set(measurement["cards"].get("target_components", {})) != {"P", "N"}:
            raise CampaignError(f"{measurement['id']} SIDIS requires proton and neutron components")
    if "nominal" not in measurement["families"]:
        raise CampaignError(f"{measurement['id']} has no nominal family")
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


def _variation_points(args: argparse.Namespace) -> list[tuple[int, int, float]]:
    polarized = _member_selector(args.polarized_pdf_members, args.profile)
    unpolarized = _member_selector(args.unpolarized_pdf_members, args.profile)
    scales = _scale_selector(args.scales, args.profile)
    # The two PDF ensembles are varied independently; a Cartesian product
    # would not represent the prescribed uncertainty construction.
    points = {(0, 0, 1.0)}
    points.update((member, 0, 1.0) for member in polarized)
    points.update((0, member, 1.0) for member in unpolarized)
    points.update((0, 0, scale) for scale in scales)
    return sorted(points, key=lambda item: (item[2] != 1.0, item[0] != 0,
                                             item[1] != 0, item))


def _scale_token(value: float) -> str:
    return {0.5: "0p5", 1.0: "1", 2.0: "2"}[float(value)]


def _analysis_instance(
    measurement: Mapping[str, Any], family_id: str
) -> str:
    """Return the exact Rivet analysis identifier written into YODA paths."""

    name = str(measurement["analysis"]["name"])
    options = measurement["families"][family_id].get(
        "analysis_options", {}
    )
    if not options:
        return name
    return name + ":" + ":".join(
        f"{key}={value}" for key, value in sorted(options.items())
    )


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
    if process_kind == "polarized_sidis":
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
                "STAR jet generator cuts must be selected from 3, 4, and 5 GeV"
            )
    else:
        if requested_jet_cut is not None:
            raise CampaignError(
                "--jet-kt-min-gev is only valid for polarized pp jet measurements"
            )
        jet_kt_min_gev = None

    options = {
        "profile": args.profile,
        "families": _family_selector(getattr(args, "families", None), measurement),
        "jobs": int(args.jobs or defaults["default_jobs"]),
        "shards": int(args.shards or defaults["default_shards"]),
        "seed_base": int(args.seed_base or defaults["default_seed_base"]),
        "lo_events": lo_events,
        "events_by_contribution": events_by_contribution,
        "smoke": smoke,
        # Keep the manifest representation JSON-native so a resumed campaign
        # compares equal after the on-disk JSON has been reloaded.
        "variation_points": [list(point) for point in _variation_points(args)],
    }
    if jet_kt_min_gev is not None:
        options["jet_kt_min_gev"] = jet_kt_min_gev
    if (options["jobs"] <= 0 or options["shards"] <= 0 or
            any(int(value) <= 0 for value in events_by_contribution.values())):
        raise CampaignError("Jobs, shards, and event counts must be positive")
    return options


def build_job_matrix(measurement: Mapping[str, Any], options: Mapping[str, Any]) -> list[dict[str, Any]]:
    families = measurement["families"]
    selected_families = list(options.get("families", ["nominal"]))
    jobs: list[dict[str, Any]] = []
    seed_slot = 0
    for family_id in selected_families:
        family = families[family_id]
        variation_points = options["variation_points"] if family_id == "nominal" else [(0, 0, 1.0)]
        for channel in measurement["channels"]:
            for polarized_member, unpolarized_member, scale in variation_points:
                targets = (measurement["cards"].get("target_components", {"none": "none"})
                           if measurement["process_kind"] == "polarized_sidis"
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
                                jobs.append({
                                    "id": job_id, "measurement": measurement["id"],
                                    "channel": channel, "target_component": target_component,
                                    "observable_level": level, "family": family_id,
                                    "family_label": family["label"], "helicity": helicity,
                                    "analysis_instance": analysis_instance,
                                    "contribution": contribution, "order": contribution,
                                    "polarized_pdf_member": polarized_member,
                                    "unpolarized_pdf_member": unpolarized_member,
                                    "scale": scale, "mpi": family["mpi"],
                                    "jet_kt_min_gev": options.get("jet_kt_min_gev"),
                                    "shard": shard, "shards": len(event_splits),
                                    "events": events,
                                    "seed": int(options["seed_base"])+seed_slot,
                                    "initial_seed": int(options["seed_base"])+seed_slot,
                                    "attempt": 0, "status": "planned", "stem": run_stem,
                                    "card_input": f"{run_stem}.in",
                                    "run_file": f"runs/{run_stem}.run",
                                    "output_yoda": f"yoda/{job_id}.yoda",
                                })
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
    if measurement["process_kind"] == "polarized_sidis":
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


def _signature(measurement: Mapping[str, Any]) -> str:
    digest = hashlib.sha256(experimental.canonical_json_bytes(
        {key: value for key, value in measurement.items() if not key.startswith("_")}))
    files = [DISPOL_ROOT/measurement["analysis"][key] for key in ("source", "info", "plot")]
    files.extend(
        DISPOL_ROOT/str(path)
        for path in measurement["analysis"].get("support_files", [])
    )
    files.append(DISPOL_ROOT/measurement["reference"]["snapshot"])
    for key in ("source_manifest", "raw_snapshot"):
        if measurement["reference"].get(key):
            files.append(DISPOL_ROOT/str(measurement["reference"][key]))
    card_dir = DISPOL_ROOT/measurement["cards"]["directory"]
    files.extend(sorted(card_dir.glob("*.in")))
    for path in files:
        digest.update(str(path.relative_to(DISPOL_ROOT)).encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


def _assert_manifest_signature_current(
    manifest: Mapping[str, Any], measurement: Mapping[str, Any]
) -> None:
    """Refuse to reinterpret shards produced by another measurement definition."""

    recorded = manifest.get("configuration", {}).get("measurement_signature")
    current = _signature(measurement)
    if recorded != current:
        raise CampaignError(
            f"{measurement['id']} campaign products were generated with a "
            "different analysis, reference, card, or metadata signature. "
            "Postprocessing and plotting are intentionally refused; use a "
            "new immutable tag and rerun the event campaign."
        )


def _card_text(measurement: Mapping[str, Any], job: Mapping[str, Any]) -> str:
    cards = measurement["cards"]
    source_helicity = job["helicity"] if job["helicity"] != "00" else "PP"
    source_stem = cards["stem_pattern"].format(
        channel=job["channel"], helicity=source_helicity,
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
    if not math.isclose(float(job["scale"]), 1.0):
        if measurement["process_kind"] == "polarized_sidis":
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
    family = measurement["families"][job["family"]]
    analysis_options = family.get("analysis_options", {})
    if analysis_options:
        analysis_name = str(measurement["analysis"]["name"])
        analysis_instance = _analysis_instance(
            measurement, str(job["family"])
        )
        text = re.sub(
            rf"(insert\s+/Herwig/Analysis/Rivet:Analyses\s+\d+\s+)"
            rf"{re.escape(analysis_name)}(?::\S+)?",
            rf"\g<1>{analysis_instance}",
            text,
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
    return {"measurement": measurement["id"], "tag": args.tag,
            "measurement_signature": _signature(measurement), **plan["options"]}


def prepare_pp(args: argparse.Namespace, measurement: Mapping[str, Any]) -> Path:
    plan = _plan(measurement, args)
    campaign_dir = _campaign_dir(measurement["id"], args.tag)
    if args.dry_run:
        print(json.dumps(plan, indent=2, sort_keys=True))
        return campaign_dir
    runtime = _runtime(measurement)
    configuration = _manifest_configuration(measurement, args, plan)
    manifest_path = campaign_dir/experimental.MANIFEST_NAME
    existing = _load_json(manifest_path) if manifest_path.exists() else None
    if existing is not None and existing.get("configuration") != configuration:
        raise CampaignError("An incompatible manifest exists; use a new immutable tag")
    for directory in ("build", "cards", "runs", "yoda", "logs", "work", "postprocess", "plots"):
        (campaign_dir/directory).mkdir(parents=True, exist_ok=True)
    snapshot = validate_vendored(str(measurement["id"]))
    reference_yoda = DISPOL_ROOT/measurement["analysis"]["reference_yoda"]
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
    analysis_instances = {
        str(job.get("analysis_instance", analysis)) for job in jobs
    }
    if len(analysis_instances) != 1:
        raise CampaignError(
            f"Shards disagree on their Rivet analysis instance: "
            f"{sorted(analysis_instances)}"
        )
    analysis_instance = next(iter(analysis_instances))
    shards = [experimental.read_histogram_series(
        campaign_dir/str(job["output_yoda"]),
        f"/{analysis_instance}/{object_name}",
    )
        for job in jobs]
    return experimental.combine_shard_series(shards,
        [int(job["events"]) for job in jobs])


def _logical_sidis_groups(
    manifest: Mapping[str, Any], campaign_dir: Path
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
    for helicities in groups.values():
        if set(helicities) != set(sidis.HELICITIES):
            raise CampaignError("SIDIS postprocessing requires PP, PM, MP, and MM")
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
) -> dict[str, experimental.BinSeries]:
    """Load shards and add normalized NLO components for every helicity."""

    output: dict[str, experimental.BinSeries] = {}
    for helicity in sidis.HELICITIES:
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
    if measurement["postprocessor"] == "hermes_sidis":
        dataset = snapshot.get("datasets", {}).get(observable)
        if dataset and dataset.get("observable") == "A_parallel":
            return str(dataset["rivet_path"])
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
    snapshot = validate_vendored(str(measurement["id"]))
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
    if measurement is not None:
        expected_points = measurement.get(
            "comparison_policy", {}
        ).get("covariance_points")
        if expected_points is not None and len(retained_indices) != int(expected_points):
            raise CampaignError(
                f"STAR primary covariance retained {len(retained_indices)} points, "
                f"expected {int(expected_points)}"
            )
    if not retained_indices:
        return {"points": 0, "status": "no finite theory bins"}

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
    return True


def postprocess_sidis(
    args: argparse.Namespace, measurement: Mapping[str, Any]
) -> Path:
    campaign_dir = _campaign_dir(measurement["id"], args.tag)
    manifest_path = campaign_dir / experimental.MANIFEST_NAME
    if not manifest_path.exists():
        raise CampaignError(f"No prepared campaign at {campaign_dir}")
    manifest = _load_json(manifest_path)
    _assert_manifest_signature_current(manifest, measurement)
    groups = _logical_sidis_groups(manifest, campaign_dir)
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

    snapshot = validate_vendored(str(measurement["id"]))
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
    snapshot = validate_vendored(str(measurement["id"]))
    analysis = str(measurement["analysis"]["name"])
    predictions: dict[tuple[Any, ...], dict[str, Any]] = {}
    raw_samples: dict[tuple[Any, ...], dict[str, Mapping[str, experimental.BinSeries]]] = {}
    for key, helicity_jobs in groups.items():
        family, channel, polarized, unpolarized, scale, mpi = key
        required = set(measurement["families"][family]["helicities"])
        if set(helicity_jobs) != required:
            raise CampaignError(f"Incomplete helicity matrix for {_variation_id(key)}")
        channel_spec = measurement["channels"][channel]
        object_names = (channel_spec.get("raw_objects") or
                        {"yield": channel_spec["raw_object"]})
        objects: dict[str, Mapping[str, experimental.BinSeries]] = {}
        for observable, object_name in object_names.items():
            objects[observable] = {
                helicity: _load_series(jobs, campaign_dir, analysis, str(object_name))
                for helicity, jobs in helicity_jobs.items()
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
        if required != set(DENOMINATOR):
            continue
        if measurement["postprocessor"] == "star_weak_bosons":
            predictions[key] = _star_prediction(channel, objects["yield"])
        elif measurement["postprocessor"] == "phenix_prompt_photon":
            predictions[key] = _phenix_prediction(objects)
        elif measurement["postprocessor"] == "star_jet_all":
            predictions[key] = _star_jet_prediction(objects)
            _apply_star_display_binning(predictions[key], snapshot)
        else:
            raise CampaignError(f"Unknown postprocessor {measurement['postprocessor']}")

    nominal_mpi = measurement["families"]["nominal"]["mpi"]
    central_keys = [("nominal", channel, 0, 0, 1.0, nominal_mpi)
                    for channel in measurement["channels"]]
    missing_central = [key for key in central_keys if key not in predictions]
    if missing_central:
        raise CampaignError(f"Missing central predictions: {missing_central}")
    bands = aggregate_uncertainties(predictions, measurement)
    correlated_goodness_of_fit = None
    if measurement["postprocessor"] == "star_jet_all":
        correlated_goodness_of_fit = _star_correlated_goodness_of_fit(
            predictions[central_keys[0]], snapshot, measurement
        )
    yoda = experimental._import_yoda()
    include_diagnostics = bool(
        getattr(args, "include_diagnostics", False)
    )
    objects_by_family: dict[str, list[Any]] = {
        family: [] for family in manifest["configuration"]["families"]
    }
    objects_by_family.setdefault("nominal", [])
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
                "JetKtMinGeV": manifest["configuration"].get(
                    "jet_kt_min_gev"
                ),
            }
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
    prediction_path = output_dir/"prediction.yoda"
    prediction_entries: list[dict[str, str]] = []
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
    if not experimental._nonempty(prediction_path):
        raise CampaignError("Nominal postprocessing produced no prediction objects")
    summary = {
        "measurement": measurement["id"], "tag": args.tag,
        "hard_process_accuracy": "LO",
        "central_sample_label": measurement["families"]["nominal"]["label"],
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
    }
    experimental.atomic_write_json(output_dir/"summary.json", summary)
    if central_rows:
        with (output_dir/"central.csv").open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(central_rows[0]))
            writer.writeheader(); writer.writerows(central_rows)
    manifest["postprocess"] = {
        "created_at": experimental.utc_now(),
        "prediction": str(prediction_path.relative_to(campaign_dir)),
        "predictions": prediction_entries,
        "summary": "postprocess/summary.json",
        "include_diagnostics": include_diagnostics,
    }
    manifest["updated_at"] = experimental.utc_now()
    manifest["history"].append({"at": experimental.utc_now(), "action": "postprocess"})
    experimental.atomic_write_json(manifest_path, manifest)
    print(f"Wrote normalized-helicity predictions to {prediction_path}")
    return prediction_path


def plot_pp(args: argparse.Namespace, measurement: Mapping[str, Any]) -> Path:
    campaign_dir = _campaign_dir(measurement["id"], args.tag)
    manifest_path = campaign_dir/experimental.MANIFEST_NAME
    if not manifest_path.exists():
        raise CampaignError(f"No prepared campaign at {campaign_dir}")
    manifest = _load_json(manifest_path)
    _assert_manifest_signature_current(manifest, measurement)
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
    runtime = manifest.get("runtime") or _runtime(measurement)
    output = campaign_dir/"plots"/"html"
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
            snapshot = validate_vendored(str(measurement["id"]))
            for dataset in snapshot["datasets"].values():
                if dataset.get("alternate_projection"):
                    command.extend(
                        ["-M", re.escape(str(dataset["rivet_path"]))]
                    )
    for prediction, label, _family_id in predictions:
        command.append(f"{prediction}:Title={label}")
    if args.dry_run:
        print(" ".join(command))
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
    experimental._run_logged(command, campaign_dir, environment,
                             campaign_dir/"logs"/"rivet-mkhtml-generate.log")
    plot_scripts = sorted(
        path for path in output.rglob("*.py")
        if not path.name.endswith("__data.py")
    )
    if not plot_scripts:
        raise CampaignError(f"rivet-mkhtml generated no plot scripts below {output}")
    script_log = campaign_dir/"logs"/"rivet-plot-scripts.log"
    rendered_scripts: list[Path] = []
    snapshot = validate_vendored(str(measurement["id"]))
    summary_path = campaign_dir/"postprocess"/"summary.json"
    summary = _load_json(summary_path) if summary_path.is_file() else {}
    with script_log.open("w", encoding="utf-8") as log:
        for script in plot_scripts:
            log.write(f"script: {script}\n")
            log.flush()
            experimental.add_theory_uncertainty_overlay(
                script,
                _pp_theory_uncertainty_bands(
                    measurement, snapshot, summary, script.stem
                ),
                nominal_prediction.name,
            )
            experimental.add_experimental_error_overlay(
                script,
                _pp_reference_overlay_points(measurement, snapshot, script.stem),
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
    manifest["plots"] = {"created_at": experimental.utc_now(),
                          "index": str(index.relative_to(campaign_dir))}
    manifest["updated_at"] = experimental.utc_now()
    manifest["history"].append({"at": experimental.utc_now(), "action": "plot"})
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
            "STAR jet generator cut in GeV; only the pinned 3, 4, and 5 GeV "
            "scan points are accepted"
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
            return experimental.main(_legacy_arguments(args))
        if getattr(args, "comparisons", False):
            args.profile = "paper"
            if getattr(args, "families", None) is None:
                args.families = "all"
        if args.command == "fetch-data":
            outputs = fetch_and_validate(
                args.measurement,
                source_file=getattr(args, "source_file", None),
            )
            print(f"Checksum- and schema-validated {len(outputs)} source(s); refreshed reference YODA")
        elif args.command == "prepare":
            prepare_pp(args, measurement)
        elif args.command == "campaign":
            run_pp(args, measurement)
        elif args.command == "postprocess":
            if measurement["process_kind"] == "polarized_sidis":
                postprocess_sidis(args, measurement)
            else:
                postprocess_pp(args, measurement)
        elif args.command == "plot":
            plot_pp(args, measurement)
        elif args.command == "full":
            prepare_pp(args, measurement)
            if not args.dry_run:
                run_pp(args, measurement)
                if measurement["process_kind"] == "polarized_sidis":
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
