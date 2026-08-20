#!/usr/bin/env python3
"""File-level provenance for Herwig Polarized production campaigns.

The campaign manifests need more than executable paths: a rebuilt source tree
does not establish which shared objects an already-installed ``Herwig`` loads.
This module records cryptographic identities for the executable, repository,
runtime libraries, Rivet plugin, and complete PDF-set inventories.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any, Mapping, Sequence


class ProvenanceError(RuntimeError):
    """Raised when a runtime or source-control identity cannot be established."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_record(path: Path) -> dict[str, Any]:
    resolved = path.resolve()
    if not resolved.is_file():
        raise ProvenanceError(f"Required runtime artifact is missing: {resolved}")
    return {
        "path": str(resolved),
        "size": resolved.stat().st_size,
        "sha256": sha256_file(resolved),
    }


def _git(repository: Path, *arguments: str) -> str:
    try:
        completed = subprocess.run(
            ["git", "-C", str(repository), *arguments],
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        raise ProvenanceError(
            f"Could not establish Git provenance for {repository}: {exc}"
        ) from exc
    return completed.stdout.strip()


def source_control_record(repository: Path) -> dict[str, Any]:
    root = Path(_git(repository, "rev-parse", "--show-toplevel")).resolve()
    branch = _git(root, "symbolic-ref", "--quiet", "--short", "HEAD")
    tracked_status = _git(root, "status", "--porcelain", "--untracked-files=no")
    origin_url = _git(root, "remote", "get-url", "origin")
    return {
        "repository": str(root),
        "commit": _git(root, "rev-parse", "HEAD"),
        "branch": branch,
        "tracked_clean": not bool(tracked_status),
        "origin_url": origin_url,
    }


def pdf_set_inventory(data_directory: Path, set_name: str) -> dict[str, Any]:
    root = (data_directory / set_name).resolve()
    if not root.is_dir():
        raise ProvenanceError(f"Required PDF set is missing: {root}")
    files: list[dict[str, Any]] = []
    inventory = hashlib.sha256()
    for path in sorted(candidate for candidate in root.rglob("*") if candidate.is_file()):
        relative = path.relative_to(root).as_posix()
        size = path.stat().st_size
        digest = sha256_file(path)
        files.append({"path": relative, "size": size, "sha256": digest})
        inventory.update(f"{relative}\0{size}\0{digest}\n".encode("utf-8"))
    if not files:
        raise ProvenanceError(f"PDF set contains no files: {root}")
    member_files = [item for item in files if item["path"].endswith(".dat")]
    return {
        "set": set_name,
        "directory": str(root),
        "file_count": len(files),
        "member_file_count": len(member_files),
        "inventory_sha256": inventory.hexdigest(),
        "files": files,
    }


def runtime_record(
    *,
    repository: Path,
    tools: Mapping[str, str],
    herwig_prefix: Path,
    artifact_paths: Mapping[str, Path],
    pdf_sets: Sequence[str],
    lhapdf_data_directory: Path,
) -> dict[str, Any]:
    artifacts = {
        label: file_record(path) for label, path in sorted(artifact_paths.items())
    }
    inventories = {
        name: pdf_set_inventory(lhapdf_data_directory, name)
        for name in sorted(set(pdf_sets))
    }
    return {
        "source_control": source_control_record(repository),
        "herwig_prefix": str(herwig_prefix.resolve()),
        "tools": dict(sorted(tools.items())),
        "artifacts": artifacts,
        "pdf_inventories": inventories,
    }


def inventory_digest(records: Mapping[str, Mapping[str, Any]]) -> str:
    """Return a stable digest for a set of already-recorded artifacts."""

    payload = json.dumps(records, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()
