#!/usr/bin/env python3
"""Immutable pilot/central/paper controller for the 2026 SIDIS tranches."""

from __future__ import annotations

import argparse
import csv
import fcntl
import hashlib
import json
import math
import os
import shlex
import subprocess
import sys
import tarfile
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


CONTROL_DIR = Path(__file__).resolve().parent
REPOSITORY = CONTROL_DIR.parents[2]
RUNTIME_DIR = CONTROL_DIR / "runtime"
CONFIG_PATH = CONTROL_DIR / "campaigns.json"
LOCK_PATH = CONTROL_DIR / "runtime-lock.json"
RUNNER = REPOSITORY / "scripts" / "run_phenomenology_campaign.py"
REFERENCE = REPOSITORY / "scripts" / "phenomenology_reference_data.py"
PROFILES = ("pilot", "pilot2", "central", "paper")


class ControllerError(RuntimeError):
    """A user-facing campaign-control failure."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ControllerError(f"Could not read {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise ControllerError(f"Expected a JSON object in {path}")
    return value


def canonical_json(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def immutable_json(path: Path, payload: Mapping[str, Any]) -> str:
    encoded = canonical_json(payload)
    digest = hashlib.sha256(encoded).hexdigest()
    if path.exists():
        if path.read_bytes() != encoded:
            raise ControllerError(f"Immutable plan already exists with different content: {path}")
    else:
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_suffix(path.suffix + ".tmp")
        temporary.write_bytes(encoded)
        os.replace(temporary, path)
    checksum = path.with_suffix(path.suffix + ".sha256")
    checksum_payload = f"{digest}  {path.name}\n"
    if checksum.exists() and checksum.read_text() != checksum_payload:
        raise ControllerError(f"Immutable checksum changed: {checksum}")
    if not checksum.exists():
        checksum.write_text(checksum_payload)
    return digest


def verify_immutable_json(path: Path) -> dict[str, Any]:
    checksum = path.with_suffix(path.suffix + ".sha256")
    if not path.is_file() or not checksum.is_file():
        raise ControllerError(f"Missing checksum-pinned event plan {path}")
    expected = checksum.read_text().split()[0]
    actual = sha256_file(path)
    if actual != expected:
        raise ControllerError(f"Event-plan checksum mismatch: {actual} != {expected}")
    return load_json(path)


def command_output(command: Sequence[str], *, cwd: Path = REPOSITORY) -> str:
    try:
        completed = subprocess.run(
            list(command), cwd=cwd, check=True, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        detail = exc.stderr.strip() if isinstance(exc, subprocess.CalledProcessError) else str(exc)
        raise ControllerError(f"Command failed: {shlex.join(command)}: {detail}") from exc
    return completed.stdout.strip()


def run(command: Sequence[str]) -> None:
    print(f"+ {shlex.join(command)}", flush=True)
    try:
        subprocess.run(list(command), cwd=REPOSITORY, check=True)
    except (OSError, subprocess.CalledProcessError) as exc:
        raise ControllerError(f"Command failed: {shlex.join(command)}") from exc


def git(*arguments: str) -> str:
    return command_output(["git", "-C", str(REPOSITORY), *arguments])


def configuration() -> dict[str, Any]:
    config = load_json(CONFIG_PATH)
    analyses = config.get("analyses")
    if not isinstance(analyses, list) or len(analyses) != 5:
        raise ControllerError("Campaign configuration must contain five analyses")
    pilot_total = 0
    central_total = 0
    central_shards = 0
    for spec in analyses:
        logical = int(spec["logical_jobs"])
        if logical % 2:
            raise ControllerError("Each SIDIS matrix must contain POSNLO/NEGNLO pairs")
        target_helicity_jobs = logical // 2
        pilot_total += target_helicity_jobs * (
            int(spec["pilot"]["posnlo_events"]) + int(spec["pilot"]["negnlo_events"])
        )
        central_total += target_helicity_jobs * (
            int(spec["central_floor"]["posnlo_events"])
            + int(spec["central_floor"]["negnlo_events"])
        )
        central_shards += logical * int(spec["central_floor"]["shards"])
    expected = config["expected"]
    if (pilot_total, central_total, central_shards) != (
        int(expected["pilot_events"]), int(expected["central_floor_events"]),
        int(expected["central_floor_shards"]),
    ):
        raise ControllerError("Configured campaign totals do not match the frozen expectations")
    return config


def select_analyses(selector: str) -> list[dict[str, Any]]:
    analyses = [dict(item) for item in configuration()["analyses"]]
    if selector == "all":
        return analyses
    selected = [item for item in analyses if item["measurement"] == selector]
    if not selected:
        raise ControllerError(f"Unknown SIDIS analysis selector {selector!r}")
    return selected


def plan_path(profile: str) -> Path:
    return RUNTIME_DIR / "plans" / f"event-plan-{profile}.json"


def _static_entries(profile: str) -> list[dict[str, Any]]:
    entries: list[dict[str, Any]] = []
    for spec in configuration()["analyses"]:
        source = spec["pilot"] if profile == "pilot" else spec["central_floor"]
        runner_profile = "paper" if profile == "paper" else "central"
        polarized = "all" if profile == "paper" and spec["measurement"] == "COMPASS_2010_I862410" else "central"
        unpolarized = "all" if profile == "paper" else "central"
        scales = "all" if profile == "paper" else "central"
        offset = {"pilot": 0, "central": 10_000_000, "paper": 20_000_000}[profile]
        entries.append({
            "measurement": spec["measurement"],
            "tag": f"sidis-{profile}-20260828-v1",
            "seed_base": int(spec["seed_base"]) + offset,
            "posnlo_events": int(source["posnlo_events"]),
            "negnlo_events": int(source["negnlo_events"]),
            "shards": int(source["shards"]),
            "jobs": int(configuration()["jobs"]),
            "profile": runner_profile,
            "polarized_pdf_members": polarized,
            "unpolarized_pdf_members": unpolarized,
            "scales": scales,
            "logical_jobs": int(spec["logical_jobs"]),
        })
    return entries


def event_plan(profile: str, *, require_gate: bool = False) -> dict[str, Any]:
    if profile == "pilot2":
        return verify_immutable_json(plan_path(profile))
    if profile == "central" and plan_path(profile).is_file():
        plan = verify_immutable_json(plan_path(profile))
        if require_gate:
            assessment = load_json(RUNTIME_DIR / "assessments" / "latest-pilot.json")
            if not assessment.get("extrapolation_accepted"):
                raise ControllerError("The latest pilot assessment did not permit extrapolation")
            if assessment.get("event_plan_sha256") != sha256_file(plan_path(profile)):
                raise ControllerError("Pilot assessment does not pin the current central event plan")
        return plan
    if profile == "central" and require_gate:
        raise ControllerError("Run 'assess pilot all' (or pilot2) before central preparation")
    entries = _static_entries(profile)
    return {
        "schema_version": 1,
        "profile": profile,
        "configuration_sha256": sha256_file(CONFIG_PATH),
        "entries": entries,
    }


def selected_plan(profile: str, selector: str, *, require_gate: bool = False) -> list[dict[str, Any]]:
    allowed = {item["measurement"] for item in select_analyses(selector)}
    entries = event_plan(profile, require_gate=require_gate).get("entries", [])
    selected = [dict(item) for item in entries if item.get("measurement") in allowed]
    if len(selected) != len(allowed):
        raise ControllerError(f"Event plan {profile} does not cover selector {selector}")
    return selected


def campaign_directory(spec: Mapping[str, Any]) -> Path:
    return REPOSITORY / "campaigns" / "phenomenology" / str(spec["measurement"]) / str(spec["tag"])


def generation_arguments(spec: Mapping[str, Any]) -> list[str]:
    return [
        "--measurement", str(spec["measurement"]), "--tag", str(spec["tag"]),
        "--posnlo-events", str(spec["posnlo_events"]),
        "--negnlo-events", str(spec["negnlo_events"]),
        "--shards", str(spec["shards"]), "--jobs", str(spec["jobs"]),
        "--seed-base", str(spec["seed_base"]), "--profile", str(spec["profile"]),
        "--polarized-pdf-members", str(spec["polarized_pdf_members"]),
        "--unpolarized-pdf-members", str(spec["unpolarized_pdf_members"]),
        "--scales", str(spec["scales"]),
    ]


def stage_command(stage: str, spec: Mapping[str, Any], *, dry_run: bool = False,
                  recover: bool = False) -> list[str]:
    command = [sys.executable, str(RUNNER), stage]
    if stage in {"prepare", "campaign", "full"}:
        command.extend(generation_arguments(spec))
        if dry_run:
            command.append("--dry-run")
        if recover:
            command.append("--recover-failed")
    else:
        command.extend(["--measurement", str(spec["measurement"]), "--tag", str(spec["tag"])])
    return command


def _pdf_inventory(data_directory: Path, set_name: str) -> dict[str, Any]:
    root = data_directory / set_name
    if not root.is_dir():
        raise ControllerError(f"Missing PDF set {root}")
    inventory = hashlib.sha256()
    files = members = 0
    for path in sorted(candidate for candidate in root.rglob("*") if candidate.is_file()):
        relative = path.relative_to(root).as_posix()
        inventory.update(f"{relative}\0{path.stat().st_size}\0{sha256_file(path)}\n".encode())
        files += 1
        members += relative.endswith(".dat")
    return {"file_count": files, "member_file_count": members,
            "inventory_sha256": inventory.hexdigest()}


def verify() -> dict[str, Any]:
    lock = load_json(LOCK_PATH)
    if git("status", "--porcelain"):
        raise ControllerError("The checkout is dirty; runtime-locked work is refused")
    origin = git("remote", "get-url", "origin")
    if origin != lock["repository"]["origin_url"]:
        raise ControllerError(f"Unexpected origin remote {origin!r}")
    commit = git("rev-parse", "HEAD")
    remote_line = command_output(["git", "ls-remote", "origin", lock["repository"]["required_remote_ref"]])
    remote_commit = remote_line.split()[0] if remote_line else ""
    if commit != remote_commit:
        raise ControllerError(f"HEAD {commit} differs from GitHub main {remote_commit or 'missing'}")

    herwig = lock["herwig"]
    executable = Path(command_output(["bash", "-c", "command -v Herwig"])).resolve()
    if executable.parent.parent.resolve() != Path(herwig["prefix"]).resolve():
        raise ControllerError(f"Active Herwig is outside locked prefix: {executable}")
    artifacts: dict[str, Any] = {}
    for label, record in herwig["artifacts"].items():
        path = Path(record["path"])
        actual = sha256_file(path) if path.is_file() else "missing"
        if actual != record["sha256"]:
            raise ControllerError(f"Runtime hash mismatch for {label}: {actual}")
        artifacts[label] = {"path": str(path), "sha256": actual, "size": path.stat().st_size}
    rivet_version = command_output(["rivet", "--version"])
    if rivet_version != herwig["rivet_version"]:
        raise ControllerError(f"Rivet version mismatch: {rivet_version!r}")
    pdf_directory = Path(command_output(["lhapdf-config", "--datadir"]))
    for set_name, expected in lock["pdf_sets"].items():
        actual = _pdf_inventory(pdf_directory, set_name)
        if any(actual[key] != expected[key] for key in actual):
            raise ControllerError(f"PDF inventory mismatch for {set_name}")

    measurements = [item["measurement"] for item in configuration()["analyses"]]
    for measurement in measurements:
        command = [sys.executable, str(REFERENCE), measurement]
        if measurement == "HERMES_2013_I1208547":
            command.append("--full-covariance")
        command_output(command)
    listing = command_output([sys.executable, str(RUNNER), "list"])
    missing = [measurement for measurement in measurements if measurement not in listing]
    if missing:
        raise ControllerError(f"Registry listing is missing {missing}")
    report = {
        "passed": True, "verified_at": utc_now(), "repository": str(REPOSITORY),
        "source_commit": commit, "origin_commit": remote_commit,
        "loaded_modules": os.environ.get("LOADEDMODULES", ""),
        "herwig_executable": str(executable), "rivet_version": rivet_version,
        "artifacts": artifacts, "runtime_lock_sha256": sha256_file(LOCK_PATH),
        "source_measurements": measurements,
    }
    atomic_json(RUNTIME_DIR / "verify.json", report)
    print(f"Verified synchronized source, data, registry, Herwig, Rivet and HwMEDIS at {commit}")
    return report


def write_profile_plan(profile: str) -> dict[str, Any]:
    if plan_path(profile).is_file():
        existing = verify_immutable_json(plan_path(profile))
        if profile == "central":
            event_plan(profile, require_gate=True)
        if existing.get("source_commit") not in (None, git("rev-parse", "HEAD")):
            raise ControllerError("Existing event plan was created for another source commit")
        return existing
    plan = event_plan(profile, require_gate=profile == "central")
    payload = dict(plan)
    payload["source_commit"] = git("rev-parse", "HEAD")
    payload["runtime_lock_sha256"] = sha256_file(LOCK_PATH)
    immutable_json(plan_path(profile), payload)
    return payload


def action_dry_run(profile: str, selector: str) -> None:
    expected_paper = configuration()["expected"]["paper_jobs"]
    output = []
    for spec in selected_plan(profile, selector):
        logical_jobs = (
            int(expected_paper[spec["measurement"]])
            if profile == "paper" else int(spec["logical_jobs"])
        )
        variation_factor = logical_jobs // int(spec["logical_jobs"])
        target_helicity_jobs = int(spec["logical_jobs"]) // 2
        total_events = variation_factor * target_helicity_jobs * (
            int(spec["posnlo_events"]) + int(spec["negnlo_events"])
        )
        output.append({
            "measurement": spec["measurement"], "profile": profile,
            "logical_jobs": logical_jobs,
            "shard_jobs": logical_jobs * int(spec["shards"]),
            "total_events": total_events,
            "prepare_command": shlex.join(stage_command("prepare", spec)),
            "launch_command": shlex.join(stage_command("campaign", spec)),
        })
    print(json.dumps({"dry_run": True, "campaigns": output}, indent=2))


def action_prepare(profile: str, selector: str) -> None:
    if profile == "paper" and selector == "all":
        raise ControllerError("Paper profiles must be prepared one analysis at a time")
    verification = verify()
    plan = write_profile_plan(profile)
    if plan.get("source_commit") != verification["source_commit"]:
        raise ControllerError("Event plan and verified source commit differ")
    for spec in selected_plan(profile, selector, require_gate=profile == "central"):
        run(stage_command("prepare", spec))


def _manifest(spec: Mapping[str, Any]) -> dict[str, Any]:
    path = campaign_directory(spec) / "manifest.json"
    if not path.is_file():
        raise ControllerError(f"Missing campaign manifest {path}")
    return load_json(path)


def action_launch(profile: str, selector: str, *, recover: bool = False) -> None:
    if profile == "paper" and selector == "all":
        raise ControllerError("Paper profiles must be launched one analysis at a time")
    verify()
    for spec in selected_plan(profile, selector, require_gate=profile == "central"):
        manifest = _manifest(spec)
        if not recover and manifest.get("status") != "prepared":
            raise ControllerError(f"Campaign is not prepared: {spec['measurement']}")
        run(stage_command("campaign", spec, recover=recover))


def action_status(profile: str, selector: str) -> None:
    for spec in selected_plan(profile, selector):
        path = campaign_directory(spec) / "manifest.json"
        if not path.is_file():
            print(f"{spec['measurement']:28s} missing")
            continue
        manifest = load_json(path)
        counts: dict[str, int] = {}
        for job in manifest.get("jobs", []):
            status = str(job.get("status", "planned"))
            counts[status] = counts.get(status, 0) + 1
        print(f"{spec['measurement']:28s} {manifest.get('status')} {json.dumps(counts, sort_keys=True)}")


def action_stage(stage: str, profile: str, selector: str) -> None:
    verify()
    for spec in selected_plan(profile, selector, require_gate=profile == "central"):
        manifest = _manifest(spec)
        if manifest.get("status") != "complete":
            raise ControllerError(f"Cannot {stage} incomplete campaign {spec['measurement']}")
        run(stage_command(stage, spec))


def _finite_float(value: Any) -> float | None:
    if value in (None, "", "None", "nan", "NaN"):
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def assess_one(spec: Mapping[str, Any]) -> dict[str, Any]:
    directory = campaign_directory(spec)
    summary_path = directory / "postprocess" / "summary.json"
    csv_path = directory / "postprocess" / "central.csv"
    if not summary_path.is_file() or not csv_path.is_file():
        raise ControllerError(f"Postprocess {spec['measurement']} before assessment")
    summary = load_json(summary_path)
    declared_masked = sum(
        len(value) for value in summary.get("masked_bins", {}).values()
    )
    csv_masked = 0
    zero_information_bins = 0
    ratios: list[float] = []
    rows = 0
    with csv_path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            rows += 1
            theory = _finite_float(row.get("theory"))
            mc = _finite_float(row.get("mc_stat"))
            data = _finite_float(row.get("data"))
            stat = _finite_float(row.get("data_stat"))
            syst = _finite_float(row.get("data_systematic"))
            if None in (theory, mc, stat, syst):
                csv_masked += 1
                continue
            # A populated experimental multiplicity with a zero prediction and
            # zero variance means the pilot observed no contributing hadron,
            # not that its Monte Carlo precision is infinite. Such a bin has
            # no basis for a 1/sqrt(N) extrapolation and must trigger pilot2.
            if data not in (None, 0.0) and theory == 0.0 and mc == 0.0:
                csv_masked += 1
                zero_information_bins += 1
                continue
            experimental = math.hypot(float(stat), float(syst))
            if experimental <= 0.0:
                csv_masked += 1
                continue
            ratios.append(float(mc) / experimental)
    if rows == 0:
        raise ControllerError(f"No primary rows in {csv_path}")
    ratios.sort()
    if ratios:
        q90 = ratios[max(0, math.ceil(0.9 * len(ratios)) - 1)]
        maximum = ratios[-1]
        fraction_half = sum(value <= 0.5 for value in ratios) / len(ratios)
        all_unit = maximum <= 1.0
        scale = max(1.0, (q90 / 0.5) ** 2, maximum**2)
    else:
        q90 = maximum = fraction_half = None
        all_unit = False
        scale = math.inf
    # The summary and CSV describe the same primary bins.  Taking their maximum
    # catches either representation without counting one masked bin twice.
    masked = max(declared_masked, csv_masked)
    return {
        "measurement": spec["measurement"], "rows": rows,
        "finite_primary_bins": len(ratios), "masked_primary_bins": masked,
        "summary_masked_primary_bins": declared_masked,
        "csv_masked_primary_bins": csv_masked,
        "zero_information_primary_bins": zero_information_bins,
        "fraction_mc_le_half_exp": fraction_half, "maximum_mc_over_exp": maximum,
        "q90_mc_over_exp": q90, "current_sample_passes": (
            masked == 0 and fraction_half is not None and fraction_half >= 0.9 and all_unit
        ),
        "recommended_scale_from_this_sample": scale,
        "manifest_sha256": sha256_file(directory / "manifest.json"),
        "summary_sha256": sha256_file(summary_path), "central_csv_sha256": sha256_file(csv_path),
    }


def _recommended_central_entries(profile_specs: Sequence[Mapping[str, Any]],
                                 results: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    floors = {item["measurement"]: item for item in _static_entries("central")}
    recommendations: list[dict[str, Any]] = []
    for pilot, result in zip(profile_specs, results):
        floor = dict(floors[str(pilot["measurement"])])
        scale = float(result["recommended_scale_from_this_sample"])
        shards = int(floor["shards"])
        floor_neg = int(floor["negnlo_events"])
        proposed_neg = int(math.ceil(int(pilot["negnlo_events"]) * scale / shards) * shards)
        neg = max(floor_neg, proposed_neg)
        floor["negnlo_events"] = neg
        floor["posnlo_events"] = 10 * neg
        recommendations.append(floor)
    return recommendations


def _doubled_pilot_entries(profile_specs: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    entries = []
    for source in profile_specs:
        item = dict(source)
        item["tag"] = "sidis-pilot2-20260828-v1"
        item["seed_base"] = int(item["seed_base"]) + 30_000_000
        item["posnlo_events"] = 2 * int(item["posnlo_events"])
        item["negnlo_events"] = 2 * int(item["negnlo_events"])
        item["shards"] = 2 * int(item["shards"])
        entries.append(item)
    return entries


def action_assess(profile: str, selector: str) -> None:
    if profile not in {"pilot", "pilot2", "central"}:
        raise ControllerError("Only pilot, pilot2, or central campaigns can be assessed")
    if profile in {"pilot", "pilot2"} and selector != "all":
        raise ControllerError("Pilot extrapolation requires the complete five-analysis selector 'all'")
    specs = selected_plan(profile, selector, require_gate=profile == "central")
    results = [assess_one(spec) for spec in specs]
    masked = sum(int(item["masked_primary_bins"]) for item in results)
    sample_passed = all(bool(item["current_sample_passes"]) for item in results)
    report: dict[str, Any] = {
        "schema_version": 1, "created_at": utc_now(), "profile": profile,
        "selector": selector, "source_commit": git("rev-parse", "HEAD"),
        "criteria": {
            "masked_primary_bins": 0,
            "fraction_mc_le_half_exp": 0.9,
            "maximum_mc_over_exp": 1.0,
        },
        "results": results, "sample_passed": sample_passed,
        "extrapolation_accepted": False,
    }
    if profile in {"pilot", "pilot2"}:
        if masked:
            if profile == "pilot2":
                raise ControllerError("Pilot2 still has masked bins; no further automatic extrapolation is allowed")
            payload = {
                "schema_version": 1, "profile": "pilot2",
                "configuration_sha256": sha256_file(CONFIG_PATH),
                "source_commit": git("rev-parse", "HEAD"),
                "runtime_lock_sha256": sha256_file(LOCK_PATH),
                "reason": "masked pilot bins; statistics doubled without extrapolation",
                "inputs": [{"measurement": item["measurement"], "summary_sha256": item["summary_sha256"]} for item in results],
                "entries": _doubled_pilot_entries(specs),
            }
            digest = immutable_json(plan_path("pilot2"), payload)
            report["next_action"] = "prepare pilot2 all"
            report["pilot2_plan_sha256"] = digest
        else:
            payload = {
                "schema_version": 1, "profile": "central",
                "configuration_sha256": sha256_file(CONFIG_PATH),
                "source_commit": git("rev-parse", "HEAD"),
                "runtime_lock_sha256": sha256_file(LOCK_PATH),
                "source_profile": profile,
                "inputs": [{"measurement": item["measurement"], "summary_sha256": item["summary_sha256"], "central_csv_sha256": item["central_csv_sha256"]} for item in results],
                "entries": _recommended_central_entries(specs, results),
            }
            digest = immutable_json(plan_path("central"), payload)
            report["extrapolation_accepted"] = True
            report["event_plan_sha256"] = digest
            report["next_action"] = "prepare central all"
    report_path = RUNTIME_DIR / "assessments" / f"{profile}-{selector}.json"
    atomic_json(report_path, report)
    if profile in {"pilot", "pilot2"}:
        atomic_json(RUNTIME_DIR / "assessments" / "latest-pilot.json", report)
    print(json.dumps(report, indent=2, sort_keys=True))


def action_package(profile: str, selector: str) -> None:
    verify()
    specs = selected_plan(profile, selector, require_gate=profile == "central")
    output_dir = RUNTIME_DIR / "packages"
    output_dir.mkdir(parents=True, exist_ok=True)
    archive = output_dir / f"sidis-{profile}-{selector}-{git('rev-parse', '--short=12', 'HEAD')}.tar.gz"
    if archive.exists():
        raise ControllerError(f"Refusing to overwrite package {archive}")
    included: list[dict[str, Any]] = []
    with tarfile.open(archive, "w:gz") as stream:
        for spec in specs:
            root = campaign_directory(spec)
            for relative in (Path("manifest.json"), Path("postprocess"), Path("plots")):
                path = root / relative
                if not path.exists():
                    continue
                stream.add(path, arcname=Path(spec["measurement"]) / relative, recursive=True)
            for path in sorted(root.glob("postprocess/*")):
                if path.is_file():
                    included.append({"path": str(path), "sha256": sha256_file(path), "size": path.stat().st_size})
    inventory = {"created_at": utc_now(), "archive": str(archive),
                 "archive_sha256": sha256_file(archive), "files": included}
    atomic_json(archive.with_suffix(archive.suffix + ".json"), inventory)
    print(f"Wrote {archive} ({inventory['archive_sha256']})")


@contextmanager
def action_lock() -> Iterable[None]:
    RUNTIME_DIR.mkdir(parents=True, exist_ok=True)
    with (RUNTIME_DIR / "controller.lock").open("w") as stream:
        try:
            fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise ControllerError("Another SIDIS tranche controller is active") from exc
        yield


def action_commands() -> None:
    print("Controller commands are documented in README.md. Equivalent raw runner commands:\n")
    for profile in ("pilot", "central", "paper"):
        print(f"# {profile}")
        for spec in selected_plan(profile, "all"):
            print(shlex.join(stage_command("prepare", spec)))
            print(shlex.join(stage_command("campaign", spec)))


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("action", choices=(
        "verify", "commands", "dry-run", "prepare", "launch", "status",
        "recover", "postprocess", "assess", "plot", "package",
    ))
    result.add_argument("profile", nargs="?", choices=PROFILES)
    result.add_argument("selector", nargs="?", default="all")
    return result


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        if args.action == "verify":
            verify()
            return 0
        if args.action == "commands":
            action_commands()
            return 0
        if args.profile is None:
            raise ControllerError(f"{args.action} requires a profile and selector")
        if args.action == "dry-run":
            action_dry_run(args.profile, args.selector)
        elif args.action == "status":
            action_status(args.profile, args.selector)
        else:
            with action_lock():
                if args.action == "prepare":
                    action_prepare(args.profile, args.selector)
                elif args.action == "launch":
                    action_launch(args.profile, args.selector)
                elif args.action == "recover":
                    action_launch(args.profile, args.selector, recover=True)
                elif args.action in {"postprocess", "plot"}:
                    action_stage(args.action, args.profile, args.selector)
                elif args.action == "assess":
                    action_assess(args.profile, args.selector)
                elif args.action == "package":
                    action_package(args.profile, args.selector)
                else:
                    raise ControllerError(f"Unsupported action {args.action}")
    except ControllerError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
