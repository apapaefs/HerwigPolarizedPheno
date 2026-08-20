#!/usr/bin/env python3
"""Locked validation and production controller for the compatibility rerun."""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import math
import os
import shlex
import subprocess
import sys
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
STABILITY = REPOSITORY / "scripts" / "check_star_generator_cut_stability.py"
PACKAGER = REPOSITORY / "scripts" / "package_corrected_comparison_plots.py"


class ControllerError(RuntimeError):
    """A user-facing validation or production-control failure."""


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


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def inventory_digest(records: Mapping[str, Mapping[str, Any]]) -> str:
    payload = json.dumps(records, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def atomic_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(temporary, path)


def command_output(command: Sequence[str], *, cwd: Path = REPOSITORY) -> str:
    try:
        completed = subprocess.run(
            list(command), cwd=cwd, check=True, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        detail = exc.stderr.strip() if isinstance(exc, subprocess.CalledProcessError) else str(exc)
        raise ControllerError(
            f"Command failed: {shlex.join(command)}: {detail}"
        ) from exc
    return completed.stdout.strip()


def run(command: Sequence[str], *, cwd: Path = REPOSITORY) -> None:
    print(f"+ {shlex.join(command)}", flush=True)
    try:
        subprocess.run(list(command), cwd=cwd, check=True)
    except (OSError, subprocess.CalledProcessError) as exc:
        raise ControllerError(f"Command failed: {shlex.join(command)}") from exc


def git(*arguments: str) -> str:
    return command_output(["git", "-C", str(REPOSITORY), *arguments])


def _pdf_inventory(data_directory: Path, set_name: str) -> dict[str, Any]:
    root = data_directory / set_name
    if not root.is_dir():
        raise ControllerError(f"Missing PDF set {root}")
    inventory = hashlib.sha256()
    files = 0
    members = 0
    for path in sorted(candidate for candidate in root.rglob("*") if candidate.is_file()):
        relative = path.relative_to(root).as_posix()
        digest = sha256_file(path)
        size = path.stat().st_size
        inventory.update(f"{relative}\0{size}\0{digest}\n".encode("utf-8"))
        files += 1
        members += relative.endswith(".dat")
    return {
        "directory": str(root.resolve()),
        "file_count": files,
        "member_file_count": members,
        "inventory_sha256": inventory.hexdigest(),
    }


def verify(*, production: bool) -> dict[str, Any]:
    lock = load_json(LOCK_PATH)
    repository_lock = lock["repository"]
    branch = git("symbolic-ref", "--quiet", "--short", "HEAD")
    expected_branch = (
        repository_lock["production_branch"]
        if production else repository_lock["validation_branch"]
    )
    if branch != expected_branch:
        raise ControllerError(
            f"Expected branch {expected_branch!r}, found {branch!r}"
        )
    if production and REPOSITORY.resolve() != Path(
        repository_lock["canonical_path"]
    ).resolve():
        raise ControllerError(
            f"Production requires the canonical checkout {repository_lock['canonical_path']}"
        )
    if git("status", "--porcelain"):
        raise ControllerError("The checkout is dirty; production and validation are refused")
    if git("remote", "get-url", "origin") != repository_lock["origin_url"]:
        raise ControllerError("The origin remote is not the canonical GitHub repository")
    commit = git("rev-parse", "HEAD")
    remote_line = command_output(
        ["git", "ls-remote", "origin", f"refs/heads/{expected_branch}"]
    )
    remote_commit = remote_line.split()[0] if remote_line else ""
    if remote_commit != commit:
        raise ControllerError(
            f"HEAD {commit} differs from origin/{expected_branch} {remote_commit or 'missing'}"
        )

    herwig_lock = lock["herwig"]
    executable = Path(command_output(["bash", "-c", "command -v Herwig"])).resolve()
    prefix = executable.parent.parent.resolve()
    if str(prefix) != herwig_lock["prefix"]:
        raise ControllerError(
            f"Active Herwig prefix {prefix} differs from {herwig_lock['prefix']}"
        )
    artifacts: dict[str, Any] = {}
    for label, expected in herwig_lock["artifacts"].items():
        path = Path(expected["path"]).resolve()
        if not path.is_file():
            raise ControllerError(f"Missing locked runtime artifact {path}")
        digest = sha256_file(path)
        if digest != expected["sha256"]:
            raise ControllerError(
                f"Runtime hash mismatch for {label}: {digest} != {expected['sha256']}"
            )
        artifacts[label] = {
            "path": str(path), "size": path.stat().st_size, "sha256": digest
        }
    rivet_version = command_output(["rivet", "--version"])
    if rivet_version != herwig_lock["rivet_version"]:
        raise ControllerError(
            f"Rivet version {rivet_version!r} differs from the lock"
        )
    pdf_directory = Path(command_output(["lhapdf-config", "--datadir"]))
    pdfs: dict[str, Any] = {}
    for set_name, expected in lock["pdf_sets"].items():
        actual = _pdf_inventory(pdf_directory, set_name)
        for key in ("file_count", "member_file_count", "inventory_sha256"):
            if actual[key] != expected[key]:
                raise ControllerError(
                    f"PDF inventory mismatch for {set_name}/{key}: "
                    f"{actual[key]} != {expected[key]}"
                )
        pdfs[set_name] = actual
    report = {
        "verified_at": utc_now(),
        "production_mode": production,
        "repository": {
            "path": str(REPOSITORY.resolve()), "branch": branch,
            "commit": commit, "origin_commit": remote_commit,
            "origin_url": repository_lock["origin_url"], "clean": True,
        },
        "herwig_prefix": str(prefix),
        "rivet_version": rivet_version,
        "artifacts": artifacts,
        "pdf_sets": pdfs,
        "runtime_lock_sha256": sha256_file(LOCK_PATH),
    }
    atomic_json(RUNTIME_DIR / "verify.json", report)
    print(f"Verified Git/runtime lock at {commit}")
    return report


def configuration() -> dict[str, Any]:
    return load_json(CONFIG_PATH)


def specs(group: str) -> list[dict[str, Any]]:
    value = configuration().get(group, [])
    if not isinstance(value, list):
        raise ControllerError(f"Campaign group {group} is not a list")
    return [dict(item) for item in value]


def campaign_directory(spec: Mapping[str, Any]) -> Path:
    return (
        REPOSITORY / "campaigns" / str(spec["storage"])
        / str(spec["measurement"]) / str(spec["tag"])
    )


def generation_arguments(spec: Mapping[str, Any]) -> list[str]:
    arguments = [
        "--measurement", str(spec["measurement"]),
        "--tag", str(spec["tag"]),
        "--profile", "central",
        "--jobs", str(spec["jobs"]),
        "--shards", str(spec["shards"]),
        "--seed-base", str(spec["seed_base"]),
    ]
    for key, option in (
        ("posnlo_events", "--posnlo-events"),
        ("negnlo_events", "--negnlo-events"),
        ("lo_events", "--lo-events"),
        ("jet_kt_min_gev", "--jet-kt-min-gev"),
    ):
        if key in spec:
            arguments.extend([option, str(spec[key])])
    return arguments


def stage_command(
    stage: str, spec: Mapping[str, Any], *, recover: bool = False
) -> list[str]:
    command = [sys.executable, str(RUNNER), stage]
    if stage in {"prepare", "campaign", "full"}:
        command.extend(generation_arguments(spec))
        if recover:
            command.append("--recover-failed")
    else:
        command.extend(
            ["--measurement", str(spec["measurement"]), "--tag", str(spec["tag"])]
        )
    return command


def _manifest(spec: Mapping[str, Any]) -> dict[str, Any]:
    path = campaign_directory(spec) / "manifest.json"
    if not path.is_file():
        raise ControllerError(f"Missing prepared manifest {path}")
    return load_json(path)


def _verify_recorded_file(record: Mapping[str, Any]) -> None:
    path = Path(str(record["path"]))
    if not path.is_file() or path.stat().st_size != int(record["size"]):
        raise ControllerError(f"Prepared artifact is missing or changed: {path}")
    if sha256_file(path) != record["sha256"]:
        raise ControllerError(f"Prepared artifact hash changed: {path}")


def verify_manifest_provenance(manifest: Mapping[str, Any], commit: str) -> None:
    recorded_runtime = manifest.get("runtime", {}).get("provenance", {})
    source = recorded_runtime.get("source_control", {})
    if source.get("commit") != commit or not source.get("tracked_clean"):
        raise ControllerError("Manifest source commit is not the clean current checkout")
    runtime_lock = load_json(LOCK_PATH)
    recorded_artifacts = recorded_runtime.get("artifacts", {})
    for label, expected in runtime_lock["herwig"]["artifacts"].items():
        actual = recorded_artifacts.get(label, {})
        if (
            actual.get("sha256") != expected["sha256"]
            or Path(str(actual.get("path", ""))).resolve()
            != Path(expected["path"]).resolve()
        ):
            raise ControllerError(
                f"Manifest runtime artifact does not match the lock: {label}"
            )
    recorded_pdfs = recorded_runtime.get("pdf_inventories", {})
    for set_name, expected in runtime_lock["pdf_sets"].items():
        actual = recorded_pdfs.get(set_name, {})
        for key in ("file_count", "member_file_count", "inventory_sha256"):
            if actual.get(key) != expected[key]:
                raise ControllerError(
                    f"Manifest PDF inventory does not match the lock: "
                    f"{set_name}/{key}"
                )
    plugin = manifest.get("plugin_provenance")
    if not isinstance(plugin, Mapping):
        raise ControllerError("Manifest does not pin the Rivet plugin")
    _verify_recorded_file(plugin)
    prepared = manifest.get("prepared_artifacts", {}).get("files")
    if not isinstance(prepared, Mapping) or not prepared:
        raise ControllerError("Manifest does not contain a prepared-artifact inventory")
    for record in prepared.values():
        _verify_recorded_file(record)
    recorded_digest = manifest.get("prepared_artifacts", {}).get(
        "inventory_sha256"
    )
    if recorded_digest != inventory_digest(prepared):
        raise ControllerError("Prepared-artifact inventory digest is inconsistent")


def require_validation(commit: str) -> dict[str, Any]:
    report_path = RUNTIME_DIR / "validation-report.json"
    report = load_json(report_path)
    if not report.get("passed") or report.get("source_commit") != commit:
        raise ControllerError(
            f"Validation report {report_path} does not certify current commit {commit}"
        )
    if report.get("runtime_lock_sha256") != sha256_file(LOCK_PATH):
        raise ControllerError("The runtime lock changed after validation")
    return report


def assert_prepared(spec: Mapping[str, Any], commit: str) -> dict[str, Any]:
    manifest = _manifest(spec)
    if manifest.get("status") != "prepared":
        raise ControllerError(
            f"{spec['measurement']}/{spec['tag']} is {manifest.get('status')}, not prepared"
        )
    jobs = manifest.get("jobs", [])
    if len(jobs) != int(spec["expected_shards"]):
        raise ControllerError(
            f"{spec['measurement']}/{spec['tag']} has {len(jobs)} shards, "
            f"expected {spec['expected_shards']}"
        )
    if any(job.get("status") != "planned" for job in jobs):
        raise ControllerError("A prepared production manifest contains a non-pending job")
    for job in jobs:
        output = campaign_directory(spec) / str(job["output_yoda"])
        if output.is_file() and output.stat().st_size:
            raise ControllerError(f"Production YODA exists before launch: {output}")
    verify_manifest_provenance(manifest, commit)
    return manifest


def _production_processes() -> list[str]:
    completed = subprocess.run(
        ["pgrep", "-af", "Herwig run"], text=True,
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
    )
    tags = {
        str(spec["tag"])
        for group in ("fixed", "star_cut_scan", "star510")
        for spec in specs(group)
    }
    return [
        line for line in completed.stdout.splitlines()
        if any(tag in line for tag in tags)
    ]


def assert_handoff(commit: str) -> dict[str, Any]:
    summaries = []
    production_shards = 0
    scan_shards = 0
    for group in ("fixed", "star_cut_scan", "star510"):
        for spec in specs(group):
            manifest = assert_prepared(spec, commit)
            count = len(manifest["jobs"])
            if group == "star_cut_scan":
                scan_shards += count
            else:
                production_shards += count
            summaries.append(
                {
                    "group": group, "measurement": spec["measurement"],
                    "tag": spec["tag"], "status": manifest["status"],
                    "shards": count,
                    "prepared_artifacts_sha256": manifest["prepared_artifacts"]["inventory_sha256"],
                    "plugin_sha256": manifest["plugin_provenance"]["sha256"],
                }
            )
    processes = _production_processes()
    if processes:
        raise ControllerError("A production Herwig process is already running")
    if production_shards != 7600 or scan_shards != 1200:
        raise ControllerError(
            f"Unexpected handoff counts: {production_shards} production, {scan_shards} scan"
        )
    report = {
        "passed": True, "created_at": utc_now(), "source_commit": commit,
        "production_shards": production_shards, "cut_scan_shards": scan_shards,
        "manifests": summaries, "production_herwig_processes": processes,
        "hard_stop": "No production campaign command was executed during preparation",
    }
    atomic_json(RUNTIME_DIR / "handoff-report.json", report)
    return report


def _finite_count(value: Any) -> int:
    if isinstance(value, Mapping):
        return sum(_finite_count(item) for item in value.values())
    if isinstance(value, list):
        return sum(_finite_count(item) for item in value)
    if isinstance(value, (int, float)) and math.isfinite(float(value)):
        return 1
    return 0


def smoke_spec(measurement: str, index: int, commit: str) -> dict[str, Any]:
    storage = "experimental" if (REPOSITORY / "config" / "experimental" / f"{measurement}.json").is_file() else "phenomenology"
    result: dict[str, Any] = {
        "measurement": measurement,
        "storage": storage,
        "tag": f"compat-smoke-20260819-{commit[:12]}",
        "seed_base": 280_000_000 + index * 10_000,
        "shards": 1,
        "jobs": 16,
    }
    if measurement == "STAR_2022_I1949588":
        result["jet_kt_min_gev"] = 4.0
    return result


def validate_smokes(verification: Mapping[str, Any]) -> dict[str, Any]:
    commit = str(verification["repository"]["commit"])
    measurements = configuration()["smoke_measurements"]
    reports = []
    for index, measurement in enumerate(measurements):
        spec = smoke_spec(str(measurement), index, commit)
        command = [sys.executable, str(RUNNER), "full"] + generation_arguments(spec)
        command.append("--smoke")
        run(command)
        manifest = _manifest(spec)
        if manifest.get("status") != "complete":
            raise ControllerError(f"Smoke campaign did not complete: {measurement}")
        jobs = manifest.get("jobs", [])
        if not jobs or any(job.get("status") != "success" for job in jobs):
            raise ControllerError(f"Smoke campaign has an incomplete matrix: {measurement}")
        for job in jobs:
            output = campaign_directory(spec) / str(job["output_yoda"])
            if not output.is_file() or output.stat().st_size == 0:
                raise ControllerError(f"Smoke campaign has an empty YODA: {output}")
        verify_manifest_provenance(manifest, commit)
        summary_path = campaign_directory(spec) / "postprocess" / "summary.json"
        summary = load_json(summary_path)
        finite = _finite_count(summary)
        if finite == 0:
            raise ControllerError(f"Smoke campaign has no finite summary values: {measurement}")
        reports.append(
            {
                "measurement": measurement, "tag": spec["tag"],
                "manifest_sha256": sha256_file(campaign_directory(spec) / "manifest.json"),
                "summary_sha256": sha256_file(summary_path), "finite_summary_values": finite,
                "shards": len(jobs), "diagnostic_only": measurement == "HERMES_2007_I726689_LEGACY",
            }
        )
    return {"smokes": reports}


@contextmanager
def action_lock() -> Iterable[None]:
    RUNTIME_DIR.mkdir(parents=True, exist_ok=True)
    with (RUNTIME_DIR / "controller.lock").open("w", encoding="utf-8") as stream:
        try:
            fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise ControllerError("Another compatibility controller is active") from exc
        yield


def action_validate() -> None:
    verification = verify(production=False)
    run(["make", "test"])
    run(["make", "check-rivet"])
    smoke_report = validate_smokes(verification)
    report = {
        "passed": True, "created_at": utc_now(),
        "source_commit": verification["repository"]["commit"],
        "runtime_lock_sha256": verification["runtime_lock_sha256"],
        "python_suite": "passed", "rivet_compile_and_registration": "passed",
        **smoke_report,
        "contracts": [
            "exact COMPASS 2010 x bins", "prompt fixed-target leptons",
            "HERMES deuteron and sparse-3D mapping", "STAR covariance and nuisance profiling",
            "complete helicity/NLO matrices", "nonempty YODA and finite summaries",
        ],
    }
    atomic_json(RUNTIME_DIR / "validation-report.json", report)
    print(f"Validated commit {report['source_commit']}")


def action_prepare() -> None:
    verification = verify(production=True)
    commit = str(verification["repository"]["commit"])
    require_validation(commit)
    for group in ("fixed", "star_cut_scan", "star510"):
        for spec in specs(group):
            run(stage_command("prepare", spec))
    report = assert_handoff(commit)
    print(
        f"Prepared {report['production_shards']} production and "
        f"{report['cut_scan_shards']} cut-scan shards; none were launched"
    )


def _require_prepared_group(group: str, commit: str) -> None:
    for spec in specs(group):
        assert_prepared(spec, commit)


def action_campaign(group: str, *, postprocess: bool = False, recover: bool = False) -> None:
    verification = verify(production=True)
    commit = str(verification["repository"]["commit"])
    require_validation(commit)
    if not recover:
        _require_prepared_group(group, commit)
    for spec in specs(group):
        run(stage_command("campaign", spec, recover=recover))
        if postprocess:
            run(stage_command("postprocess", spec))


def action_check_star_cut() -> None:
    verification = verify(production=True)
    commit = str(verification["repository"]["commit"])
    require_validation(commit)
    scan = sorted(specs("star_cut_scan"), key=lambda item: item["jet_kt_min_gev"])
    command = [sys.executable, str(STABILITY)]
    for spec in scan:
        command.extend(
            [f"--campaign-{int(spec['jet_kt_min_gev'])}", str(campaign_directory(spec))]
        )
    command.extend(
        [
            "--output-json", str(RUNTIME_DIR / "star-cut-stability.json"),
            "--output-markdown", str(RUNTIME_DIR / "star-cut-stability.md"),
        ]
    )
    run(command)


def _require_star_gate(commit: str) -> None:
    report_path = RUNTIME_DIR / "star-cut-stability.json"
    report = load_json(report_path)
    if not report.get("gate_passed") or report.get("source_commit") != commit:
        raise ControllerError("The STAR 3-vs-4 GeV generator-cut gate has not passed")
    for spec in specs("star_cut_scan"):
        label = f"{float(spec['jet_kt_min_gev']):g}gev"
        recorded = report["inputs"][label]
        directory = campaign_directory(spec)
        if sha256_file(directory / "manifest.json") != recorded["manifest_sha256"]:
            raise ControllerError(f"STAR cut manifest changed after the gate: {label}")
        if sha256_file(directory / "postprocess" / "summary.json") != recorded["summary_sha256"]:
            raise ControllerError(f"STAR cut summary changed after the gate: {label}")


def action_star510() -> None:
    verification = verify(production=True)
    commit = str(verification["repository"]["commit"])
    require_validation(commit)
    _require_star_gate(commit)
    _require_prepared_group("star510", commit)
    for spec in specs("star510"):
        run(stage_command("campaign", spec))


def action_postprocess_or_plot(stage: str) -> None:
    verification = verify(production=True)
    commit = str(verification["repository"]["commit"])
    require_validation(commit)
    for group in ("fixed", "star510"):
        for spec in specs(group):
            manifest = _manifest(spec)
            if manifest.get("status") != "complete":
                raise ControllerError(f"Cannot {stage} incomplete campaign {spec['measurement']}")
            run(stage_command(stage, spec))


def action_package() -> None:
    verification = verify(production=True)
    commit = str(verification["repository"]["commit"])
    require_validation(commit)
    run(
        [
            sys.executable, str(PACKAGER), "--campaign-config", str(CONFIG_PATH),
            "--output-root", str(REPOSITORY / "figures"),
            "--source-commit", commit,
        ]
    )


def action_status() -> None:
    try:
        verification = verify(production=True)
        sync = f"ok ({verification['repository']['commit']})"
    except ControllerError as exc:
        sync = f"FAILED: {exc}"
    print(f"Repository/runtime synchronization: {sync}")
    for group in ("fixed", "star_cut_scan", "star510"):
        for spec in specs(group):
            path = campaign_directory(spec) / "manifest.json"
            if not path.is_file():
                print(f"{group:14s} {spec['measurement']:28s} missing")
                continue
            manifest = load_json(path)
            counts: dict[str, int] = {}
            for job in manifest.get("jobs", []):
                status = str(job.get("status", "planned"))
                counts[status] = counts.get(status, 0) + 1
            print(
                f"{group:14s} {spec['measurement']:28s} "
                f"{manifest.get('status')} {json.dumps(counts, sort_keys=True)}"
            )
    processes = _production_processes()
    print(f"Production Herwig processes: {len(processes)}")
    for line in processes:
        print(f"  {line}")


def action_commands() -> None:
    print("User-facing controller commands:")
    for command in (
        '"$CONTROL/run.sh" verify', '"$CONTROL/run.sh" commands',
        '"$CONTROL/run.sh" status', '"$CONTROL/run.sh" fixed',
        '"$CONTROL/run.sh" star-cut-scan', '"$CONTROL/run.sh" check-star-cut',
        '"$CONTROL/run.sh" star510', '"$CONTROL/run.sh" recover fixed',
        '"$CONTROL/run.sh" recover star-cut-scan',
        '"$CONTROL/run.sh" recover star510', '"$CONTROL/run.sh" postprocess',
        '"$CONTROL/run.sh" plot', '"$CONTROL/run.sh" package',
    ):
        print(f"  {command}")
    print("\nUnderlying immutable campaign commands:")
    for group in ("fixed", "star_cut_scan", "star510"):
        for spec in specs(group):
            print(f"  # {group}: {spec['measurement']}")
            print(f"  {shlex.join(stage_command('prepare', spec))}")
            print(f"  {shlex.join(stage_command('campaign', spec))}")


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "action",
        choices=(
            "verify", "commands", "status", "validate", "prepare", "fixed",
            "star-cut-scan", "check-star-cut", "star510", "recover",
            "postprocess", "plot", "package",
        ),
    )
    parser.add_argument(
        "recover_group", nargs="?", choices=("fixed", "star-cut-scan", "star510")
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = make_parser().parse_args(argv)
    try:
        if args.action == "commands":
            action_commands()
        elif args.action == "status":
            action_status()
        elif args.action == "verify":
            verify(production=True)
        else:
            with action_lock():
                if args.action == "validate":
                    action_validate()
                elif args.action == "prepare":
                    action_prepare()
                elif args.action == "fixed":
                    action_campaign("fixed")
                elif args.action == "star-cut-scan":
                    action_campaign("star_cut_scan", postprocess=True)
                elif args.action == "check-star-cut":
                    action_check_star_cut()
                elif args.action == "star510":
                    action_star510()
                elif args.action == "recover":
                    if args.recover_group is None:
                        raise ControllerError("recover requires fixed, star-cut-scan, or star510")
                    group = {
                        "fixed": "fixed", "star-cut-scan": "star_cut_scan",
                        "star510": "star510",
                    }[args.recover_group]
                    action_campaign(
                        group, postprocess=group == "star_cut_scan", recover=True
                    )
                elif args.action in {"postprocess", "plot"}:
                    action_postprocess_or_plot(args.action)
                elif args.action == "package":
                    action_package()
                else:
                    raise ControllerError(f"Unsupported action {args.action}")
    except ControllerError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
