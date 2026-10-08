#!/usr/bin/env python3
"""Discover current campaign results and build a searchable offline plot browser.

Uses the standard library and the existing build_sidis_results_browser helper;
rendering registered reference-only measurements additionally uses matplotlib.
Reads campaign manifests and finished plots; never runs Herwig or changes input
campaigns. By default scans both registries, excludes smoke runs, and selects
the newest compatible, fully postprocessed result by campaign creation time.
"""
from __future__ import annotations

import argparse
from collections import Counter
import copy
from datetime import datetime, timezone
import fnmatch
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

import build_sidis_results_browser as common
import compass_open_charm_reference as open_charm

ROOT = Path(__file__).resolve().parents[1]
ASSETS = Path(__file__).with_name("results_browser_assets")
ACTIVE_STATES = {"running", "preparing", "prepared", "postprocessing", "plotting"}
SIGNATURE_POLICY = (
    "Only current-definition plots are included. A recorded presentation-only "
    "refresh is accepted only after its historical plot bytes reproduce the "
    "original generation signature."
)


def canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode()


def source_file(repo: Path, relative: str) -> Path:
    path = repo / relative
    if not common.regular_file(path, repo):
        raise common.ExportError(f"Missing source file: {path}")
    return path


def registry(repo: Path) -> dict:
    result = {}
    for kind in ("experimental", "phenomenology", "reference"):
        for path in sorted((repo / "config" / kind).glob("*.json")):
            source_file(repo, path.relative_to(repo).as_posix())
            descriptor = common.read_json(path)
            if kind == "reference":
                descriptor = open_charm.discover(repo)[path.stem]
            identifier = descriptor.get("id")
            if identifier != path.stem or not common.TOKEN.fullmatch(identifier):
                raise common.ExportError(f"Invalid registry identity: {path}")
            if identifier in result:
                raise common.ExportError(f"Duplicate registry identity: {identifier}")
            result[identifier] = (kind, descriptor)
    if not result:
        raise common.ExportError(f"No measurement registries found below {repo}")
    return dict(sorted(result.items()))


def pinned_reference_paths(repo: Path, entries: object) -> list[str]:
    """Validate the runner's additional-reference/source-file checksum records."""
    if not isinstance(entries, list):
        raise common.ExportError("Reference checksum records must be a list")
    paths = []
    for entry in entries:
        if (not isinstance(entry, dict) or set(entry) != {"path", "sha256"}
                or not isinstance(entry["path"], str)
                or not re.fullmatch(r"[0-9a-f]{64}", str(entry["sha256"]))):
            raise common.ExportError("Reference files require a path and SHA-256 checksum")
        relative = entry["path"]
        if relative in paths:
            raise common.ExportError(f"Duplicate reference source path: {relative}")
        if common.digest(source_file(repo, relative)) != entry["sha256"]:
            raise common.ExportError(f"Reference source checksum mismatch: {relative}")
        paths.append(relative)
    return paths


def measurement_signature(repo: Path, descriptor: dict, kind: str,
                          plot_bytes: bytes | None = None) -> str:
    """Keep both runner signature contracts, without importing a physics runtime."""
    if kind == "reference":
        return open_charm.signature(repo, descriptor)
    payload = copy.deepcopy({k: v for k, v in descriptor.items() if not k.startswith("_")})
    if kind == "experimental":
        payload.get("campaign", {}).pop("comparison_profile", None)
    digest = hashlib.sha256(canonical(payload))
    analysis = descriptor["analysis"]
    paths = []
    for spec in [analysis, *analysis.get("companions", [])]:
        paths.extend(spec[k] for k in ("source", "info", "plot"))
        paths.extend(spec.get("support_files", []))
    reference = descriptor["reference"]
    paths.append(reference["snapshot"])
    if kind == "experimental":
        additional = pinned_reference_paths(repo, reference.get("additional_snapshots", []))
        paths.extend(additional)
        paths.extend(dict.fromkeys(
            relative for snapshot in additional
            for relative in pinned_reference_paths(
                repo, common.read_json(source_file(repo, snapshot)).get("source_files", []))
        ))
    for key in (("source_manifest", "raw_snapshot") if kind == "phenomenology" else ("raw_snapshot",)):
        if reference.get(key):
            paths.append(reference[key])
    directory = repo / descriptor["cards"]["directory"]
    if not common.contained(directory, repo) or directory.is_symlink():
        raise common.ExportError(f"Invalid card directory: {directory}")
    paths.extend(p.relative_to(repo).as_posix() for p in sorted(directory.glob("*.in")))
    for relative in paths:
        path = source_file(repo, relative)
        digest.update(relative.encode())
        digest.update(plot_bytes if plot_bytes is not None and relative == analysis["plot"] else path.read_bytes())
    return digest.hexdigest()


def comparison_signature(repo: Path, descriptor: dict) -> str:
    profile = descriptor["campaign"]["comparison_profile"]
    digest = hashlib.sha256(canonical(profile))
    paths = set()
    for family in profile["families"].values():
        for component in family.get("components", descriptor["cards"].get("components", {"": {}})):
            for helicity in family["helicities"]:
                for order in family["orders"]:
                    stem = family["stem_pattern"].format(component=component, helicity=helicity, order=order)
                    paths.add(Path(family["card_directory"]) / (stem + ".in"))
    for relative in sorted(paths):
        path = source_file(repo, relative.as_posix())
        digest.update(relative.as_posix().encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


def compatible(repo: Path, descriptor: dict, kind: str, manifest: dict, current: str) -> tuple[bool, str]:
    config = manifest.get("configuration", {})
    if kind == "experimental" and config.get("comparisons"):
        if config.get("comparison_signature") != comparison_signature(repo, descriptor):
            return False, "Comparison definitions or cards have changed."
    recorded = config.get("measurement_signature")
    if recorded == current:
        return True, "Current analysis, reference and card definitions match."
    refresh = (manifest.get("plots") or {}).get("presentation_only_refresh") or {}
    plot = descriptor["analysis"]["plot"]
    commit = str(refresh.get("generation_source_commit", ""))
    if (refresh.get("mode") == "presentation_only_plot_metadata_refresh"
            and refresh.get("current_measurement_signature") == current
            and refresh.get("generation_measurement_signature") == recorded
            and refresh.get("plot_path") == plot
            and refresh.get("current_plot_sha256") == common.digest(source_file(repo, plot))
            and re.fullmatch(r"[0-9a-f]{40,64}", commit)):
        try:
            historical = subprocess.check_output(
                ["git", "-C", str(repo), "show", f"{commit}:{plot}"], stderr=subprocess.DEVNULL,
                timeout=15,
            )
            if (hashlib.sha256(historical).hexdigest() == refresh.get("generation_plot_sha256")
                    and measurement_signature(repo, descriptor, kind, historical) == recorded):
                return True, "Verified presentation-only plot refresh; generation definitions match."
        except (OSError, subprocess.SubprocessError):
            pass
    return False, "Analysis/reference/card signature differs from the current source; plots are withheld."


def info_field(text: str, key: str) -> str:
    match = re.search(r"^" + re.escape(key) + r":\s*([^\n]*)(?:\n((?:[ \t]+[^\n]*\n?)*))?", text, re.M)
    if not match:
        return ""
    inline, block = match.group(1), match.group(2) or ""
    return " ".join((inline if inline not in {">", "|"} else "", block.strip())).strip()


def point_groups(value: object) -> list:
    """Find numerical groups in flat, dataset and nested channel snapshots."""
    groups = []
    if isinstance(value, dict):
        if isinstance(value.get("points"), list):
            groups.append(value["points"])
        for key, child in value.items():
            if key != "points":
                groups.extend(point_groups(child))
    elif isinstance(value, list):
        for child in value:
            groups.extend(point_groups(child))
    return groups


def metadata(repo: Path, identifier: str, kind: str, descriptor: dict) -> dict:
    reference = descriptor["reference"]
    snapshot = common.read_json(source_file(repo, reference["snapshot"]))
    if snapshot.get("measurement") != identifier:
        raise common.ExportError(f"Reference identity mismatch for {identifier}")
    analysis = descriptor["analysis"]
    info = source_file(repo, analysis["info"]).read_text() if analysis.get("info") else ""
    if identifier in common.CATALOG:
        _, item = common.metadata(repo, identifier)
    else:
        records = point_groups(snapshot)
        provenance = snapshot.get("provenance", {})
        sources = common.read_json(source_file(repo, reference["source_manifest"])) if reference.get("source_manifest") else {"reference": reference, "provenance": provenance}
        diagnostic = reference.get("kind") == "internal_observable_definition"
        hepdata = reference.get("format") == "hepdata-table-json" or provenance.get("hepdata_record_doi")
        authority = ("Generator diagnostic" if diagnostic else snapshot.get("provenance", {}).get("source")
                     or ("Official HERMES archive" if reference.get("format") == "tar-five-column"
                         else "HEPData + publication" if hepdata and provenance.get("publication_pdf_url")
                         else "HEPData" if hepdata else "Registered experimental reference"))
        item = {
            "id": identifier, "experiment": identifier.split("_")[0], "title": descriptor["title"],
            "data_type": info_field(info, "Summary") or analysis.get("summary") or descriptor["title"], "authority": authority,
            "authority_note": snapshot.get("provenance", {}).get("notes", ""),
            "data_available": not diagnostic,
            "reference_entries": sum(len(points) for points in records),
            "dataset_count": len(records), "links": [], "source_manifest": sources,
            "physics": descriptor.get("physics", {}), "pdf_ensembles": descriptor.get("pdf_ensembles", {}),
            "reference_metadata": {k: v for k, v in snapshot.items() if k not in {"datasets", "points", "primary_covariance", "observables", "bin_edges"}},
            "rivet_status": info_field(info, "Status") or analysis.get("status", "unspecified"),
            "metadata_sha256": {"descriptor": common.digest(repo / "config" / kind / f"{identifier}.json"), "reference_snapshot": common.digest(repo / reference["snapshot"])},
        }
    links = {url: label for label, url in item["links"]}
    def add(label: str, url: object) -> None:
        if isinstance(url, str) and re.match(r"https?://", url):
            links.setdefault(url, label)
    for key, label in (("source_url", "Numerical source"), ("hepdata_table", "HEPData table")):
        add(label, reference.get(key))
    for key, label in (("publication", "Publication"), ("official_archive", "Official numerical archive"), ("hepdata_record", "HEPData record")):
        add(label, snapshot.get("provenance", {}).get(key))
    for url in reference.get("sources", []):
        add("Reference source", url)
    provenance = snapshot.get("provenance", {})
    add("Publication", provenance.get("publication_pdf_url"))
    if provenance.get("hepdata_record_doi"):
        add("HEPData record", "https://doi.org/" + provenance["hepdata_record_doi"])
    for doi in re.findall(r"(?:DOI:|doi\s*=\s*[\"{])([^\s\"}]+)", info):
        add("Publication DOI", "https://doi.org/" + doi)
    if reference.get("record_doi"):
        add("HEPData record", "https://doi.org/" + reference["record_doi"])
    item.update({
        "links": [[label, url] for url, label in links.items()], "registry": kind,
        "process_kind": descriptor.get("process_kind", "fixed_target_dis"),
        "description": info_field(info, "Description") or analysis.get("description", ""), "run_info": info_field(info, "RunInfo"),
        "simulation": descriptor.get("simulation", {}),
        "families": descriptor.get("families", descriptor.get("campaign", {}).get("comparison_profile", {}).get("families", {})),
        "generator_cuts": descriptor.get("generator_cuts", {}),
        "global_uncertainties": descriptor.get("global_uncertainties", {}),
        "scales": descriptor.get("scales", {}),
    })
    return item


def timestamp(value: object) -> float | None:
    if not isinstance(value, str):
        return None
    try:
        dt = datetime.fromisoformat(value.replace("Z", "+00:00"))
        return (dt if dt.tzinfo else dt.replace(tzinfo=timezone.utc)).timestamp()
    except (ValueError, OverflowError):
        return None


def processes(repo: Path, proc: Path = Path("/proc")) -> dict:
    """A local observation, never inferred from a manifest's recorded state."""
    result = {"available": False, "controllers": [], "worker_directories": []}
    if not proc.is_dir():
        return result
    readable = 0
    try:
        for path in proc.iterdir():
            if not path.name.isdigit():
                continue
            try:
                if path.stat().st_uid != os.getuid():
                    continue
                args = (path / "cmdline").read_bytes().decode(errors="replace").split("\0")
                readable += 1
                names = {Path(arg).name for arg in args}
                if "Herwig" in names:
                    result["worker_directories"].append(os.readlink(path / "cwd"))
                if not names.intersection({"run_phenomenology_campaign.py", "run_validation_campaign.py", "run_experimental_campaign.py"}):
                    continue
                if not any(str(repo) in arg for arg in args):
                    continue
                if "--measurement" in args and "--tag" in args:
                    result["controllers"].append([args[args.index("--measurement") + 1], args[args.index("--tag") + 1]])
            except (OSError, IndexError):
                continue
    except OSError:
        return result
    result["available"] = bool(readable)
    return result


def activity(candidate: dict, observation: dict) -> str:
    if candidate["generation_complete"]:
        return "Generation complete"
    if not observation["available"]:
        return "Activity not checked on this host"
    if [candidate["measurement"], candidate["tag"]] in observation["controllers"]:
        return "Controller active"
    prefix = candidate["path"] + "/"
    if any(p.startswith(prefix) for p in observation["worker_directories"]):
        return "Workers active; no controller observed"
    return "Stopped — no active controller or workers observed"


def scan_campaigns(repo: Path, roots: list[Path], identifier: str, kind: str,
                   descriptor: dict, current: str, include_smoke: bool) -> tuple[list, list]:
    candidates, rejected = [], []
    seen = set()
    for root in roots:
        analysis_dir = root / identifier
        if not analysis_dir.is_dir():
            continue
        for path in sorted(analysis_dir.glob("*/manifest.json")):
            try:
                if not common.regular_file(path, root):
                    continue
                if path.resolve() in seen:
                    continue
                seen.add(path.resolve())
                d = common.read_json(path)
                tag = path.parent.name
                if d.get("measurement") != identifier or d.get("tag") != tag or not common.TOKEN.fullmatch(tag):
                    raise common.ExportError("Manifest identity does not match its directory")
                config = d.get("configuration") or {}
                if config.get("smoke") and not include_smoke:
                    rejected.append({"tag": tag, "reason": "Smoke test excluded"})
                    continue
                jobs = d.get("jobs") or []
                if not isinstance(jobs, list) or any(not isinstance(j, dict) for j in jobs):
                    raise common.ExportError("Invalid shard inventory")
                ok, reason = compatible(repo, descriptor, kind, d, current)
                created = timestamp(d.get("created_at"))
                fallback = created is None
                created = path.stat().st_mtime if fallback else created
                counts = Counter(j.get("status", "unknown") for j in jobs)
                provenance = (d.get("runtime") or {}).get("provenance") or {}
                candidates.append({
                    "measurement": identifier, "tag": tag, "path": str(path.parent),
                    "created_at": d.get("created_at"), "updated_at": d.get("updated_at"),
                    "ordering_time": created, "ordering_basis": "manifest mtime (creation time unavailable)" if fallback else "campaign created_at",
                    "recorded_status": d.get("status", "unknown"), "configuration": config,
                    "compatible": ok, "compatibility_note": reason,
                    "total_shards": len(jobs), "successful_shards": counts.get("success", 0), "shards": dict(counts),
                    "generation_complete": bool(jobs) and counts.get("success", 0) == len(jobs),
                    "events_in_successful_shards": sum(j.get("events", 0) for j in jobs if j.get("status") == "success"),
                    "source_control": provenance.get("source_control", {}),
                    "versions": {k: (d.get("runtime") or {}).get(k) for k in ("herwig_version", "rivet_version")},
                    "manifest_sha256": common.digest(path), "_manifest": d,
                })
            except (common.ExportError, OSError, ValueError, TypeError, KeyError) as exc:
                rejected.append({"path": str(path), "reason": str(exc)})
    return sorted(candidates, key=lambda c: (c["ordering_time"], c["tag"], c["path"]), reverse=True), rejected


def plot_rules(repo: Path, descriptor: dict) -> list:
    rules = []
    for spec in [descriptor["analysis"], *descriptor["analysis"].get("companions", [])]:
        current = None
        for line in source_file(repo, spec["plot"]).read_text().splitlines():
            line = line.strip().removeprefix("# ").strip()
            if line.startswith("BEGIN PLOT "):
                current = [line[len("BEGIN PLOT "):], {}]
            elif line.startswith("END PLOT"):
                if current:
                    rules.append(current)
                current = None
            elif current and "=" in line:
                key, value = line.split("=", 1)
                if key in {"Title", "XLabel", "YLabel"}:
                    current[1][key] = value.strip()
    return rules


def plot_labels(name: str, rules: list) -> dict:
    labels = {}
    for pattern, values in rules:
        match = fnmatch.fnmatchcase("/" + name, pattern)
        if not match:
            try:
                match = bool(re.fullmatch(pattern, "/" + name))
            except re.error:
                pass
        if match:
            labels.update(values)
    return labels


def postprocess_paths(candidate: dict) -> dict:
    """Honor the recorded layouts of both the legacy and native runners."""
    record = candidate["_manifest"].get("postprocess") or {}
    return {
        "prediction": record.get("prediction", record.get("yoda", "postprocess/prediction.yoda")),
        "summary": record.get("summary", record.get("summary_json", "postprocess/summary.json")),
        "csv": record.get("central_csv", record.get("summary_csv", "postprocess/central.csv")),
    }


def rendered(candidate: dict, formats: str, rules: list) -> tuple[list, list, list]:
    campaign = Path(candidate["path"])
    index = (candidate["_manifest"].get("plots") or {}).get("index")
    issues, plots, files = [], [], []
    if not index or not common.regular_file(campaign / index, campaign / "plots"):
        return plots, files, ["No finished plot index is recorded."]
    plot_root = (campaign / index).parent
    products = postprocess_paths(candidate)
    required = {products["prediction"], products["summary"]}
    for prediction in (candidate["_manifest"].get("postprocess") or {}).get("predictions", []):
        required.add(prediction.get("path", prediction.get("yoda", products["prediction"])))
    for name in sorted(required):
        if not common.regular_file(campaign / name, campaign / "postprocess"):
            issues.append(f"Postprocessing output missing: {name}.")
    pairs = {}
    for path in sorted(plot_root.rglob("*")):
        fmt = path.suffix.lower().removeprefix(".")
        if fmt not in {"png", "pdf"} or not common.regular_file(path, plot_root):
            continue
        name = path.relative_to(plot_root).with_suffix("").as_posix()
        pairs.setdefault(name, {})[fmt] = path
    for name, pair in pairs.items():
        if set(pair) != {"png", "pdf"}:
            issues.append(f"Missing PNG/PDF sibling: {name}")
        plot = {"name": name, "labels": plot_labels(name, rules), "files": {}}
        for fmt, path in pair.items():
            if formats != "both" and fmt != formats:
                continue
            relative = Path("analyses") / candidate["measurement"] / "plots" / path.relative_to(plot_root)
            plot["files"][fmt] = relative.as_posix()
            files.append((path, relative))
        if plot["files"]:
            plots.append(plot)
    if not plots:
        issues.append("No rendered plots in the requested format.")
    return plots, files, issues


def compact(value: object, limit: int = 24, depth: int = 0) -> object:
    """Keep diagnostic scalars, explicitly summarizing large arrays/mappings."""
    if depth > 7:
        return {"omitted": "nested diagnostic payload"}
    if isinstance(value, list):
        if len(value) > limit:
            return {"entries": len(value), "preview": [compact(v, limit, depth + 1) for v in value[:3]], "truncated": True}
        return [compact(v, limit, depth + 1) for v in value]
    if isinstance(value, dict):
        result = {k: compact(v, limit, depth + 1) for k, v in list(value.items())[:limit]}
        if len(value) > limit:
            result["_omitted_entries"] = len(value) - limit
        return result
    if isinstance(value, float) and not math.isfinite(value):
        return "non-finite (recorded in source summary)"
    return value


def diagnostics(path: Path) -> dict:
    data = common.read_json(path)
    masks = data.get("masked_bins")
    mask_counts = {k: len(v) for k, v in masks.items() if isinstance(v, list)} if isinstance(masks, dict) else None
    keys = ("measurement", "tag", "hard_process_accuracy", "central_sample_label", "diagnostic_only",
            "prediction_object_count", "nlo_combination", "fit_policy", "systematic_model", "correction_policy",
            "goodness_of_fit", "correlated_goodness_of_fit", "global_uncertainties", "prediction_assumption",
            "assumptions", "sign_convention", "target_combination", "uncertainties", "prediction_families")
    return {"summary_sha256": common.digest(path), "masked_bins_by_dataset": mask_counts,
            "total_masked_bins": sum(mask_counts.values()) if mask_counts is not None else None,
            "details": {k: compact(data[k]) for k in keys if k in data}}


def public_candidate(candidate: dict | None) -> dict | None:
    return {k: v for k, v in candidate.items() if not k.startswith("_")} if candidate else None


def inspect_analysis(repo: Path, roots: list[Path], identifier: str, kind: str, descriptor: dict,
                     observation: dict, formats: str, include_smoke: bool, tag: str | None) -> tuple[dict, list]:
    item = metadata(repo, identifier, kind, descriptor)
    current = measurement_signature(repo, descriptor, kind)
    if kind == "reference":
        if tag is not None:
            raise common.ExportError(f"{identifier} is reference-only; campaign selection is unavailable")
        snapshot = open_charm.validate(repo)
        item.update({"selected": None, "latest": None, "history": [], "rejected": [],
                     "plots": open_charm.plot_records(snapshot, formats),
                     "issues": [descriptor["simulation"]["reason"]], "availability": "reference",
                     "current_signature": current, "diagnostics": {}, "downloads": {}})
        files = []
        for key, label, filename in [("raw_snapshot", "Published numerical data (CSV)", "reference.csv"),
                                      ("snapshot", "Complete reference and bin definitions (JSON)", "reference.json")]:
            relative = Path("analyses") / identifier / filename
            files.append((source_file(repo, descriptor["reference"][key]), relative))
            item["downloads"][label] = relative.as_posix()
        return item, files
    candidates, rejected = scan_campaigns(repo, roots, identifier, kind, descriptor, current, include_smoke)
    for c in candidates:
        c["activity"] = activity(c, observation)
    latest = candidates[0] if candidates else None
    rules = plot_rules(repo, descriptor)
    selected, plots, files, issues = None, [], [], []
    eligible = [c for c in candidates if c["compatible"] and (tag is None or c["tag"] == tag)]
    if tag is not None and not eligible:
        raise common.ExportError(f"No compatible campaign found for explicit selection {identifier}={tag}")
    # Prefer complete result products. A newer unfinished attempt is still shown.
    partial = None
    for candidate in eligible:
        try:
            found, payload, missing = rendered(candidate, formats, rules)
        except (common.ExportError, OSError) as exc:
            rejected.append({"tag": candidate["tag"], "reason": str(exc)})
            continue
        if candidate["generation_complete"] and found and not missing:
            selected, plots, files, issues = candidate, found, payload, missing
            break
        if found and partial is None:
            partial = (candidate, found, payload, missing)
    if selected is None and partial:
        selected, plots, files, issues = partial
        issues.insert(0, "Partial result products; generation or postprocessing is incomplete.")
    if selected and latest and selected["path"] != latest["path"]:
        if tag is not None:
            issues.append(f"Explicitly selected {selected['tag']}; the latest campaign is {latest['tag']}.")
        else:
            issues.append(f"Showing {selected['tag']}; the newer campaign {latest['tag']} has no complete compatible result.")
    if not selected:
        issues.append("No current-definition plot results are available.")
    availability = "ready" if selected and selected["generation_complete"] and not issues else "partial" if selected else "unavailable"
    item.update({"selected": public_candidate(selected), "latest": public_candidate(latest),
                 "history": [public_candidate(c) for c in candidates], "rejected": rejected,
                 "plots": plots, "issues": issues, "availability": availability,
                 "current_signature": current, "diagnostics": {}, "downloads": {}})
    if selected:
        campaign = Path(selected["path"])
        products = postprocess_paths(selected)
        path = campaign / products["csv"]
        if common.regular_file(path, campaign / "postprocess"):
            relative = Path("analyses") / identifier / path.name
            files.append((path, relative))
            item["downloads"]["Numerical results (CSV)"] = relative.as_posix()
        summary = campaign / products["summary"]
        if common.regular_file(summary, campaign / "postprocess"):
            try:
                item["diagnostics"] = diagnostics(summary)
            except common.ExportError as exc:
                item["issues"].append(str(exc))
                item["availability"] = "partial"
    return item, files


def index_html(data: dict) -> str:
    payload = json.dumps(data, ensure_ascii=False, separators=(",", ":"), allow_nan=False)
    payload = payload.replace("&", "\\u0026").replace("<", "\\u003c").replace(">", "\\u003e").replace("\u2028", "\\u2028").replace("\u2029", "\\u2029")
    return '''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Herwig Polarized · Results</title><link rel="stylesheet" href="assets/browser.css"></head>
<body><a class="skip" href="#results">Skip to results</a>
<header><div class="brand"><span class="brand-mark">H</span><span>HERWIG POLARIZED <small>RESULTS LIBRARY</small></span></div>
<a href="bundle.json" class="quiet-link">Snapshot inventory ↗</a></header>
<main><section class="intro"><div><p class="eyebrow">EXPERIMENTAL COMPARISONS & GENERATOR DIAGNOSTICS</p>
<h1>Your latest results,<br><em>in one place.</em></h1><p class="lead">Browse plots, inspect analysis definitions, and see which campaigns still need attention.</p></div>
<div class="snapshot"><span class="live-dot"></span><strong>Offline snapshot</strong><p id="created"></p><small>Campaigns are checked when this browser is built.</small></div></section>
<section class="metrics" id="metrics" aria-label="Snapshot totals"></section>
<section class="toolbar" aria-label="Browse and filter results"><div class="view-switch"><button id="analyses-tab" aria-pressed="true">Analyses</button><button id="plots-tab" aria-pressed="false">All plots</button></div>
<label class="search"><span>Search</span><input id="search" type="search" placeholder="Observable, analysis, plot or campaign…"></label>
<label><span>Experiment</span><select id="experiment"><option value="">All experiments</option></select></label>
<label><span>Availability</span><select id="availability"><option value="">All results</option><option value="ready">Ready</option><option value="reference">Reference data only</option><option value="partial">Partial / newer run pending</option><option value="unavailable">No current plots</option></select></label>
<button id="reset" class="text-button">Reset</button></section>
<div class="result-heading"><h2 id="result-title">Analysis overview</h2><span id="result-count" aria-live="polite"></span></div>
<div id="active-filter"></div><section id="results" class="analysis-grid"></section>
<nav id="pagination" aria-label="Plot pages"></nav>
<p class="method-note">Ready means the recorded generation and plot products are present. It does not establish agreement with experimental data. Reference data only means published points are available but no generator comparison is enabled. Smoke tests are excluded by default; incompatible historical definitions remain visible in campaign history.</p>
<details class="scan-notes"><summary>How this snapshot was selected</summary><div id="scan-policy"></div></details>
</main><dialog id="info-dialog" aria-labelledby="info-title"><div class="dialog-header"><p class="eyebrow">ANALYSIS INFORMATION</p><button id="close-info" class="close" aria-label="Close analysis information">×</button></div><div id="info-content"></div></dialog>
<footer>Herwig Polarized Phenomenology · Local, portable, and built from recorded results.<span>PNG / PDF plots · Numerical tables · Analysis provenance</span></footer>
<noscript>This browser needs JavaScript for search and galleries. Compact per-analysis metadata are also available in the bundle inventory.</noscript>
<script id="results-data" type="application/json">''' + payload + '''</script><script src="assets/browser.js"></script></body></html>'''


def build(repo: Path, roots: list[Path], output: Path, *, identifiers: list[str] | None = None,
          selections: dict | None = None, formats: str = "both", include_smoke: bool = False,
          archive: bool = False, dry_run: bool = False, observation: dict | None = None) -> dict:
    repo, output = repo.resolve(), output.absolute()
    roots = list(dict.fromkeys(p.resolve() for p in roots))
    if output.exists() or output.is_symlink():
        raise common.ExportError(f"Refusing to overwrite existing output: {output}")
    for root in roots:
        if common.contained(output, root) or common.contained(root, output):
            raise common.ExportError("Output must be separate from every campaign input tree.")
    if formats not in {"both", "png", "pdf"}:
        raise common.ExportError("Plot format must be both, png or pdf")
    tar_path = output.with_name(output.name + ".tar.gz")
    if archive and (tar_path.exists() or tar_path.is_symlink()):
        raise common.ExportError(f"Refusing to overwrite archive: {tar_path}")
    registered = registry(repo)
    selections = selections or {}
    chosen = sorted(set(identifiers or registered) | set(selections))
    unknown = set(chosen) - set(registered)
    if unknown:
        raise common.ExportError(f"Unknown measurements: {', '.join(sorted(unknown))}")
    if any(not common.TOKEN.fullmatch(str(tag)) for tag in selections.values()):
        raise common.ExportError("Invalid explicit campaign tag")
    observation = processes(repo) if observation is None else observation
    items, payload = [], []
    for identifier in chosen:
        kind, descriptor = registered[identifier]
        item, files = inspect_analysis(repo, roots, identifier, kind, descriptor, observation, formats, include_smoke, selections.get(identifier))
        items.append(item)
        payload.extend(files)
    report = {"output": str(output), "analyses": len(items), "plots": sum(len(i["plots"]) for i in items),
              "availability": dict(Counter(i["availability"] for i in items)),
              "payload_bytes": sum(p.stat().st_size for p, _ in payload),
              "reference_plots": sum(len(i["plots"]) for i in items if i["availability"] == "reference"),
              "selected_campaigns": {i["id"]: i["selected"]["tag"] if i["selected"] else None for i in items}}
    if dry_run:
        return report
    output.parent.mkdir(parents=True, exist_ok=True)
    created = datetime.now(timezone.utc).isoformat(timespec="seconds")
    data = {"created_at": created, "repository": str(repo), "campaign_roots": [str(p) for p in roots],
            "selection_policy": "Newest compatible completed plots by campaign created_at; newer unfinished attempts are shown separately. File mtime is a labelled fallback only when created_at is absent. Partial plots are used only when no complete result is available.",
            "signature_policy": SIGNATURE_POLICY, "include_smoke": include_smoke,
            "process_observation_available": observation["available"], "items": items}
    with tempfile.TemporaryDirectory(prefix=".results-browser-", dir=output.parent) as temporary:
        staging = Path(temporary) / output.name
        staging.mkdir()
        for source, relative in payload:
            destination = staging / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, destination)
        for item in items:
            if item["availability"] == "reference":
                open_charm.render(open_charm.validate(repo), staging, formats)
                report["payload_bytes"] += sum((staging / relative).stat().st_size
                                               for plot in item["plots"] for relative in plot["files"].values())
            folder = staging / "analyses" / item["id"]
            folder.mkdir(parents=True, exist_ok=True)
            sources = item.pop("source_manifest")
            common.write_json(folder / "source-manifest.json", sources)
            item["downloads"]["Analysis metadata"] = f"analyses/{item['id']}/metadata.json"
            item["downloads"]["Reference provenance"] = f"analyses/{item['id']}/source-manifest.json"
            common.write_json(folder / "metadata.json", {k: v for k, v in item.items() if k != "plots"})
        (staging / "assets").mkdir()
        for name in ("browser.css", "browser.js"):
            shutil.copyfile(ASSETS / name, staging / "assets" / name)
        (staging / "index.html").write_text(index_html(data), encoding="utf-8")
        (staging / "README.txt").write_text(
            "Open index.html in a browser; no web server or Herwig environment is required.\n"
            "Re-run build_results_browser.py on the campaign host to discover newer outputs.\n"
            "Contains PNG/PDF figures, CSV tables, compact diagnostics and source metadata.\n"
            "Campaign logs, events, raw YODA, run files, libraries and dense covariance arrays are not copied.\n"
            "This is a viewing snapshot, not a campaign backup.\n", encoding="utf-8")
        inventory = [{"path": p.relative_to(staging).as_posix(), "bytes": p.stat().st_size, "sha256": common.digest(p)} for p in sorted(staging.rglob("*")) if p.is_file()]
        common.write_json(staging / "bundle.json", {"created_at": created, **report,
                          "selection_policy": data["selection_policy"], "signature_policy": SIGNATURE_POLICY,
                          "files": inventory})
        staging.rename(output)
    report["index"] = str(output / "index.html")
    if archive:
        result = common.write_archive(output)
        report.update(archive=str(result), archive_bytes=result.stat().st_size, archive_sha256=common.digest(result))
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=ROOT)
    parser.add_argument("--campaign-root", type=Path, action="append", help="Repeat for additional campaign roots; defaults to campaigns/{experimental,phenomenology}")
    parser.add_argument("--output", type=Path, help="New output directory; default: campaigns/exports/results-TIMESTAMP")
    parser.add_argument("--measurement", action="append", help="Repeat to limit the analysis inventory")
    parser.add_argument("--campaign", action="append", default=[], metavar="ANALYSIS=TAG", help="Pin an exact compatible campaign instead of automatic selection")
    parser.add_argument("--include-smoke", action="store_true", help="Also consider explicitly recorded smoke tests")
    parser.add_argument("--plot-format", choices=("both", "png", "pdf"), default="both")
    parser.add_argument("--tar", action="store_true", help="Create a portable .tar.gz beside the result directory")
    parser.add_argument("--dry-run", action="store_true", help="Report selections and sizes without writing")
    args = parser.parse_args(argv)
    try:
        selections = {}
        for value in args.campaign:
            identifier, separator, tag = value.partition("=")
            if not separator or not common.TOKEN.fullmatch(identifier) or not common.TOKEN.fullmatch(tag):
                raise common.ExportError(f"Expected ANALYSIS=TAG, got {value!r}")
            if identifier in selections:
                raise common.ExportError(f"Duplicate campaign override: {identifier}")
            selections[identifier] = tag
        roots = args.campaign_root or [args.repo / "campaigns" / kind for kind in ("phenomenology", "experimental")]
        output = args.output or args.repo / "campaigns/exports" / datetime.now(timezone.utc).strftime("results-%Y%m%dT%H%M%S%fZ")
        report = build(args.repo, roots, output, identifiers=args.measurement, selections=selections,
                       formats=args.plot_format, include_smoke=args.include_smoke, archive=args.tar, dry_run=args.dry_run)
        print(json.dumps(report, indent=2))
        return 0
    except (common.ExportError, OSError, ValueError, KeyError, TypeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
