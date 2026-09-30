"""Automatic, provenance-checked HERMES projection gallery for campaign plots."""
from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path
from typing import Any, Mapping

from hermes_apar_integrated import MEASUREMENT, REFERENCE, ROOT, reconstruct


CACHE_MANIFEST = "cache-manifest.json"
SETTINGS = {"model": "GD11", "quadrature_order": 32, "plots": True, "controls": True}
SOURCE_INPUTS = (
    REFERENCE,
    ROOT / "data/experimental/HERMES_2007_I726689/unpolarized-fit.json",
    *(ROOT / "data/experimental/HERMES_2007_I726689" / name for name in (
        "raw-born-apar-table-7.json", "raw-born-apar-table-8.json",
        "raw-born-apar-table-18.json", "raw-born-apar-table-19.json",
        "raw-born-apar-fig4.eps.gz")),
    *(ROOT / "scripts" / name for name in (
        "hermes_apar_integrated.py", "hermes_apar_integrated_plots.py",
        "hermes_unpolarized_fit.py", "hermes_born_reference_data.py",
        "hermes_apar_integrated_campaign.py")),
    ROOT / "config/experimental/HERMES_2007_I726689.json",
)
# Changing the analysis/card selection requires reviewing the fitted-UU cuts,
# not merely giving an old acceptance calculation a new cache key.
FIDUCIAL_SOURCES = {
    ROOT / "analyses/rivet/dis/HERMES_2007_I726689.cc":
        "29ea7fa2db2276f56a35bacb32c9b7a570c0c25876af97582251fad04beea4be",
    ROOT / "cards/experimental/HERMES_2007_I726689/HERMES_2007_I726689-Common.in":
        "9657649425143b93e3fafc8368b5c9c92cd1f47aa4453114f7d2b63416a49202",
}
REQUIRED_OUTPUTS = {"integrated.json", "integrated.csv", "index.html"} | {
    f"{stem}-{selection}.{suffix}"
    for stem in ("integrated-overview", "fit-domain-sensitivity")
    for selection in ("Q2GT1", "Q2GT4") for suffix in ("png", "pdf")
} | {
    f"AParallel_{target}_vs_{axis}_{selection}_reconstructed.{suffix}"
    for target in ("P", "D") for axis in ("x", "q2")
    for selection in ("Q2GT1", "Q2GT4") for suffix in ("png", "pdf")
}


class CacheError(RuntimeError):
    """An automatic reconstruction cannot be published as complete."""


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _source_inputs(summary_path: Path) -> dict[str, str]:
    return {str(path.resolve()): _digest(path)
            for path in (summary_path, *SOURCE_INPUTS, *FIDUCIAL_SOURCES)}


def _validate_fiducial_sources() -> None:
    for path, expected in FIDUCIAL_SOURCES.items():
        if _digest(path) != expected:
            raise CacheError(
                f"HERMES fiducial source changed: {path}. Review the independent "
                "fit's DEFAULT_CUTS and update its validated acceptance pins before reconstructing."
            )


def _valid_cache(directory: Path, key: str) -> bool:
    try:
        record = json.loads((directory / CACHE_MANIFEST).read_text(encoding="utf-8"))
        outputs = record["outputs"]
        if record["cache_key"] != key or not REQUIRED_OUTPUTS.issubset(outputs):
            return False
        for name, digest in outputs.items():
            relative = Path(name)
            if relative.is_absolute() or ".." in relative.parts:
                return False
            path = directory / relative
            if path.is_symlink() or not path.is_file() or path.stat().st_size == 0:
                return False
            if _digest(path) != digest:
                return False
        return True
    except (OSError, ValueError, KeyError, TypeError):
        return False


def ensure_campaign_plots(campaign_dir: Path) -> dict[str, Any]:
    """Reuse a complete cache or write a new revision, preserving older outputs."""
    summary_path = campaign_dir / "postprocess/summary.json"
    if not summary_path.is_file():
        raise CacheError(
            "HERMES reconstructed plots need postprocess/summary.json with the "
            "Born cells; run postprocess on the completed campaign first."
        )
    if json.loads(summary_path.read_text(encoding="utf-8")).get("measurement") != MEASUREMENT:
        raise CacheError("The reconstructed gallery requires a HERMES_2007_I726689 summary")
    _validate_fiducial_sources()
    from hermes_unpolarized_fit import MODEL_METADATA, SNAPSHOT_PATH
    if MODEL_METADATA["snapshot_sha256"] != _digest(SNAPSHOT_PATH):
        raise CacheError("The unpolarized fit changed in this Python process; restart the plot command")
    inputs = _source_inputs(summary_path)
    provenance = {"inputs": inputs, "settings": SETTINGS}
    key = hashlib.sha256(json.dumps(provenance, sort_keys=True).encode("utf-8")).hexdigest()
    cache_root = campaign_dir / "derived-integrated-data/automatic"
    prefix = key[:16]
    for directory in sorted(cache_root.glob(f"{prefix}*")):
        if directory.is_dir() and not directory.is_symlink() and _valid_cache(directory, key):
            return _gallery_info(campaign_dir, directory, key, created=False)

    cache_root.mkdir(parents=True, exist_ok=True)
    attempt = 0
    while True:
        directory = cache_root / (prefix if attempt == 0 else f"{prefix}-{attempt:03d}")
        try:
            directory.mkdir()
            break
        except FileExistsError:
            attempt += 1
    # Failure leaves this incomplete revision available for inspection; retry
    # creates a sibling. No prior reconstruction is deleted or overwritten.
    reconstruct(summary_path, directory, **SETTINGS)
    outputs = {path.relative_to(directory).as_posix(): _digest(path)
               for path in directory.rglob("*") if path.is_file()}
    if not REQUIRED_OUTPUTS.issubset(outputs) or any(
        (directory / name).stat().st_size == 0 for name in outputs
    ):
        raise CacheError(f"Incomplete reconstructed HERMES gallery at {directory}")
    if inputs != _source_inputs(summary_path):
        raise CacheError("Reconstruction inputs changed during rendering; retry the plot stage")
    record = {"schema_version": 1, "cache_key": key, **provenance, "outputs": outputs}
    (directory / CACHE_MANIFEST).write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    return _gallery_info(campaign_dir, directory, key, created=True)


def _gallery_info(campaign_dir: Path, directory: Path, key: str, *, created: bool) -> dict[str, Any]:
    return {"directory": directory.relative_to(campaign_dir).as_posix(),
            "index": (directory / "index.html").relative_to(campaign_dir).as_posix(),
            "cache_key": key, "created": created,
            "estimator": "GD11-weighted published-cell asymmetries; model-assisted",
            "theory_family": "nominal"}


def publish_gallery(campaign_dir: Path, output_dir: Path, info: Mapping[str, Any]) -> Path:
    """Copy the validated gallery inside the analysis folder for portable viewing."""
    source = campaign_dir / str(info["directory"])
    key = str(info["cache_key"])
    if source.is_symlink() or not _valid_cache(source, key):
        raise CacheError(f"Cannot publish an incomplete reconstructed gallery: {source}")
    record_bytes = (source / CACHE_MANIFEST).read_bytes()
    record = json.loads(record_bytes)
    root = output_dir / MEASUREMENT / "reconstructed-born"
    for candidate in sorted(root.glob(f"{source.name}*")):
        if (candidate.is_dir() and not candidate.is_symlink()
                and _valid_cache(candidate, key)
                and (candidate / CACHE_MANIFEST).read_bytes() == record_bytes):
            return candidate
    root.mkdir(parents=True, exist_ok=True)
    attempt = 0
    while True:
        destination = root / (source.name if attempt == 0 else f"{source.name}-view-{attempt:03d}")
        try:
            destination.mkdir()
            break
        except FileExistsError:
            attempt += 1
    # Publish only the recorded files. Keep failed or modified copies as history.
    for name in record["outputs"]:
        path = destination / name
        path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source / name, path)
    if ((source / CACHE_MANIFEST).read_bytes() != record_bytes
            or not _valid_cache(source, key)
            or any(_digest(destination / name) != expected
                   for name, expected in record["outputs"].items())):
        raise CacheError("Reconstructed gallery changed during publication; retry the plot stage")
    (destination / CACHE_MANIFEST).write_bytes(record_bytes)
    return destination
