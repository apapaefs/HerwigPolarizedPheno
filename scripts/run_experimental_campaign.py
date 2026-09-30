#!/usr/bin/env python3
"""Registry-driven fixed-target experimental campaigns for polarized DIS.

The runner deliberately keeps generator execution separate from observable
construction. Every physical helicity and signed NLO contribution is generated
independently; postprocessing combines already-normalized bins and propagates
the ordinary/1-D within-sample covariance recorded by the Rivet analysis.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import copy
import csv
import gzip
import hashlib
import html
import io
import json
import math
import os
import re
import runpy
import shutil
import subprocess
import sys
import tarfile
import tempfile
import threading
import time
import urllib.request
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Sequence

import runtime_provenance as provenance


SCRIPT_PATH = Path(__file__).resolve()
DISPOL_ROOT = SCRIPT_PATH.parents[1]
REGISTRY_DIR = DISPOL_ROOT / "config" / "experimental"
CAMPAIGN_ROOT = DISPOL_ROOT / "campaigns" / "experimental"
MANIFEST_NAME = "manifest.json"
TAG_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
PROGRESS_MARKER_RE = re.compile(
    r"event>\s+(?P<current>init|\d+)(?:\s+(?P<total>\d+)|/(?P<total_alt>\d+))"
)

try:
    from prettytable import PrettyTable
except ImportError:
    PrettyTable = None


class CampaignError(RuntimeError):
    """A user-facing campaign configuration or execution error."""


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_json_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")


def load_json(path: Path) -> dict[str, Any]:
    try:
        with path.open("r", encoding="utf-8") as stream:
            value = json.load(stream)
    except (OSError, json.JSONDecodeError) as exc:
        raise CampaignError(f"Could not read JSON file {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise CampaignError(f"Expected a JSON object in {path}")
    return value


def atomic_write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as stream:
            json.dump(value, stream, indent=2, sort_keys=True)
            stream.write("\n")
        os.replace(temporary, path)
    except Exception:
        try:
            os.unlink(temporary)
        except OSError:
            pass
        raise


def atomic_write_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as stream:
            stream.write(value)
        os.replace(temporary, path)
    except Exception:
        try:
            os.unlink(temporary)
        except OSError:
            pass
        raise


def _configured_components(measurement: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
    configured = measurement.get("cards", {}).get("components")
    if configured is None:
        return {"": {"label": "target"}}
    if not isinstance(configured, Mapping) or not configured:
        raise CampaignError(f"{measurement.get('id', 'measurement')} cards.components must be non-empty")
    components: dict[str, dict[str, Any]] = {}
    for component, metadata in configured.items():
        token = str(component)
        if not token or not TAG_PATTERN.fullmatch(token):
            raise CampaignError(f"Invalid target-component token {token!r}")
        if not isinstance(metadata, Mapping):
            raise CampaignError(f"Target component {token!r} metadata must be an object")
        components[token] = dict(metadata)
    return components


def _validate_target_component_maps(
    configured: Any, components: Iterable[str], context: str
) -> None:
    """Require explicit, finite coefficients for each independent target sample."""
    if not isinstance(configured, Mapping) or set(configured) != {"unpolarized", "longitudinal"}:
        raise CampaignError(f"{context} must define unpolarized and longitudinal target coefficients")
    labels = set(components)
    for channel, coefficients in configured.items():
        if not isinstance(coefficients, Mapping) or set(coefficients) != labels:
            raise CampaignError(f"{context} {channel} coefficients must match cards.components")
        for component, value in coefficients.items():
            try:
                coefficient = float(value)
            except (TypeError, ValueError, OverflowError) as exc:
                raise CampaignError(f"{context} {channel}/{component} coefficient must be finite") from exc
            if isinstance(value, bool) or not math.isfinite(coefficient):
                raise CampaignError(f"{context} {channel}/{component} coefficient must be finite")


def _order_combination_coefficients(
    measurement: Mapping[str, Any], orders: Iterable[str]
) -> dict[str, float]:
    """Resolve signs of stored order bins; absent coefficients preserve legacy sums."""
    configured = measurement["combination"].get("order_coefficients")
    if "order_coefficients" in measurement["combination"]:
        nominal_orders = set(measurement["cards"]["orders"])
        if not isinstance(configured, Mapping) or set(configured) != nominal_orders:
            raise CampaignError("combination.order_coefficients must match the nominal cards.orders")
        for order, value in configured.items():
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise CampaignError(f"Order coefficient {order!r} must be a finite number")
            try:
                coefficient = float(value)
            except OverflowError as exc:
                raise CampaignError(f"Order coefficient {order!r} must be a finite number") from exc
            if not math.isfinite(coefficient):
                raise CampaignError(f"Order coefficient {order!r} must be a finite number")
    return {
        str(order): float(configured.get(str(order), 1.0)) if configured is not None else 1.0
        for order in orders
    }


def _order_combination_formula(coefficients: Mapping[str, float]) -> str:
    terms: list[str] = []
    for order, coefficient in coefficients.items():
        magnitude = abs(float(coefficient))
        term = str(order) if magnitude == 1.0 else f"{magnitude:g}*{order}"
        sign = "-" if coefficient < 0.0 else ("+" if terms else "")
        terms.append(sign + term)
    return "".join(terms)


def _validate_measurement(measurement: dict[str, Any], source: Path) -> None:
    required = {"id", "analysis", "reference", "cards", "campaign", "raw_observables", "outputs", "combination"}
    missing = sorted(required - set(measurement))
    if missing:
        raise CampaignError(f"Registry entry {source} is missing: {', '.join(missing)}")
    if measurement.get("schema_version") not in {1, 2}:
        raise CampaignError(f"Unsupported schema_version in {source}")
    if source.stem != measurement["id"]:
        raise CampaignError(f"Registry filename {source.name} must match id {measurement['id']}")
    helicities = measurement["cards"].get("helicities", {})
    if set(helicities) != {"PP", "PM", "MP", "MM"}:
        raise CampaignError(f"{measurement['id']} must define PP, PM, MP, and MM")
    orders = measurement["cards"].get("orders", {})
    if set(orders) != {"POSNLO", "NEGNLO"}:
        raise CampaignError(f"{measurement['id']} must define POSNLO and NEGNLO")
    _order_combination_coefficients(measurement, orders)
    if set(measurement["raw_observables"]) != set(measurement["outputs"]):
        raise CampaignError(f"{measurement['id']} output selections must match raw-observable selections")
    components = _configured_components(measurement)
    target_combination = measurement["combination"].get("target_components")
    if target_combination is not None or len(components) > 1:
        _validate_target_component_maps(
            target_combination, components, f"{measurement['id']} combination.target_components"
        )
    for selection, output in measurement["outputs"].items():
        if not isinstance(output, Mapping):
            raise CampaignError(f"{measurement['id']} output {selection!r} must be an object")
        estimator = str(output.get("estimator", "a1"))
        if estimator not in {"a1", "a_parallel"}:
            raise CampaignError(f"Unknown estimator {estimator!r} for output {selection!r}")
        if str(output.get("axis", "x")) not in {"x", "q2"}:
            raise CampaignError(f"Unknown axis for output {selection!r}")
        names = measurement["raw_observables"][selection]
        if estimator == "a_parallel":
            if set(names) != {"ordinary"} or "a1" in output:
                raise CampaignError(
                    f"Direct A_parallel output {selection!r} must use only an ordinary histogram and no A1 path"
                )
            if "apar" not in output:
                raise CampaignError(f"Direct A_parallel output {selection!r} requires an apar path")
        if "supported_bins" in output:
            bins = output["supported_bins"]
            if (not isinstance(bins, list) or any(type(value) is not int or value < 1 for value in bins)
                    or len(set(bins)) != len(bins)):
                raise CampaignError(f"Output {selection!r} supported_bins must be unique positive bin numbers")
        if "target_components" in output:
            _validate_target_component_maps(
                output["target_components"], components,
                f"{measurement['id']} outputs.{selection}.target_components",
            )

    comparison_profile = measurement.get("campaign", {}).get("comparison_profile")
    if comparison_profile is None:
        return
    if not isinstance(comparison_profile, dict) or not isinstance(
        comparison_profile.get("families"), dict
    ):
        raise CampaignError(f"{measurement['id']} comparison_profile must define a families object")
    allowed_event_options = {"posnlo_events", "negnlo_events", "lo_events"}
    allowed_postprocessing = {"helicity_asymmetry", "direct_unpolarized"}
    for family_id, family in comparison_profile["families"].items():
        if str(family_id) == "nominal":
            raise CampaignError(
                f"{measurement['id']} comparison family id 'nominal' is reserved"
            )
        required_family = {
            "label",
            "id_prefix",
            "card_directory",
            "stem_pattern",
            "helicities",
            "orders",
            "event_options",
            "postprocess",
        }
        if not isinstance(family, dict):
            raise CampaignError(
                f"{measurement['id']} comparison family {family_id!r} must be an object"
            )
        if not required_family.issubset(family):
            missing_family = sorted(required_family - set(family))
            raise CampaignError(
                f"{measurement['id']} comparison family {family_id!r} is missing: "
                + ", ".join(missing_family)
            )
        helicity_map = family["helicities"]
        order_map = family["orders"]
        event_options = family["event_options"]
        if not isinstance(helicity_map, dict) or not helicity_map:
            raise CampaignError(f"Comparison family {family_id!r} has no helicities")
        if not isinstance(order_map, dict) or not order_map:
            raise CampaignError(f"Comparison family {family_id!r} has no orders")
        if not isinstance(event_options, dict):
            raise CampaignError(
                f"Comparison family {family_id!r} event_options must be an object"
            )
        if set(event_options) != set(order_map):
            raise CampaignError(
                f"Comparison family {family_id!r} event_options must match its orders"
            )
        if not set(event_options.values()).issubset(allowed_event_options):
            raise CampaignError(f"Comparison family {family_id!r} has an unknown event option")
        postprocess = str(family["postprocess"])
        if postprocess not in allowed_postprocessing:
            raise CampaignError(
                f"Comparison family {family_id!r} has unsupported postprocess mode {postprocess!r}"
            )
        if postprocess == "helicity_asymmetry" and set(helicity_map) != {"PP", "PM", "MP", "MM"}:
            raise CampaignError(
                f"Comparison family {family_id!r} must define PP, PM, MP, and MM"
            )
        if postprocess == "direct_unpolarized" and set(helicity_map) != {"00"}:
            raise CampaignError(
                f"Comparison family {family_id!r} must contain only the 00 helicity"
            )


def discover_registry(registry_dir: Path = REGISTRY_DIR) -> dict[str, dict[str, Any]]:
    registry: dict[str, dict[str, Any]] = {}
    for path in sorted(registry_dir.glob("*.json")):
        measurement = load_json(path)
        _validate_measurement(measurement, path)
        identifier = str(measurement["id"])
        if identifier in registry:
            raise CampaignError(f"Duplicate measurement id {identifier}")
        measurement["_registry_path"] = str(path.resolve())
        registry[identifier] = measurement
    return registry


def get_measurement(identifier: str, registry_dir: Path = REGISTRY_DIR) -> dict[str, Any]:
    registry = discover_registry(registry_dir)
    if identifier not in registry:
        available = ", ".join(sorted(registry)) or "none"
        raise CampaignError(f"Unknown measurement {identifier!r}; available: {available}")
    return registry[identifier]


def resolve_dispol_path(relative: str) -> Path:
    path = (DISPOL_ROOT / relative).resolve()
    try:
        path.relative_to(DISPOL_ROOT.resolve())
    except ValueError as exc:
        raise CampaignError(f"Registry path escapes DISPOL: {relative}") from exc
    return path


def campaign_directory(measurement_id: str, tag: str) -> Path:
    if not TAG_PATTERN.fullmatch(tag):
        raise CampaignError("Campaign tags may contain only letters, digits, '.', '_', and '-'")
    return CAMPAIGN_ROOT / measurement_id / tag


def render_table(
    headers: Sequence[str],
    rows: Sequence[Sequence[str]],
    aligns: Sequence[str] | None = None,
) -> list[str]:
    """Render the same progress tables used by the validation campaign."""

    if not rows:
        return []
    string_rows = [[str(cell) for cell in row] for row in rows]
    alignments = list(aligns) if aligns is not None else ["l"] * len(headers)
    if PrettyTable is not None:
        table = PrettyTable()
        table.field_names = list(headers)
        for header, align in zip(headers, alignments):
            table.align[header] = "r" if align == "r" else "l"
        for row in string_rows:
            table.add_row(row)
        return table.get_string().splitlines()

    widths: list[int] = []
    for index, header in enumerate(headers):
        cell_width = max((len(row[index]) for row in string_rows), default=0)
        widths.append(max(len(header), cell_width))

    def format_row(row: Sequence[str], is_header: bool = False) -> str:
        cells: list[str] = []
        for index, cell in enumerate(row):
            align = alignments[index] if index < len(alignments) else "l"
            if align == "r" and not is_header:
                cells.append(cell.rjust(widths[index]))
            else:
                cells.append(cell.ljust(widths[index]))
        return "| " + " | ".join(cells) + " |"

    separator = "+-" + "-+-".join("-" * width for width in widths) + "-+"
    output = [separator, format_row(list(headers), is_header=True), separator]
    output.extend(format_row(row) for row in string_rows)
    output.append(separator)
    return output


def campaign_monitor_dir(campaign_dir: Path) -> Path:
    return campaign_dir / "monitor"


def campaign_status_json_path(campaign_dir: Path) -> Path:
    return campaign_monitor_dir(campaign_dir) / "status.json"


def campaign_status_txt_path(campaign_dir: Path) -> Path:
    return campaign_monitor_dir(campaign_dir) / "status.txt"


def fmt_seconds_compact(value: float | int | None) -> str:
    if value is None or not isinstance(value, (int, float)) or not math.isfinite(value):
        return "n/a"
    total = int(round(float(value)))
    hours, remainder = divmod(total, 3600)
    minutes, seconds = divmod(remainder, 60)
    return f"{hours:d}:{minutes:02d}:{seconds:02d}"


def build_campaign_monitor_payload(
    *,
    measurement: Mapping[str, Any],
    tag: str,
    phase: str,
    started_at: float,
    message: str | None = None,
    herwig: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    physics = measurement.get("physics", {})
    payload: dict[str, Any] = {
        "measurement": measurement["id"],
        "tag": tag,
        "phase": phase,
        "elapsed_s": max(0.0, time.time() - started_at),
        "pdfs": {
            "unpolarized": physics.get("unpolarized_pdf", "unspecified"),
            "polarized_diff": physics.get("polarized_pdf", "unspecified"),
        },
        "prediction_assumption": physics.get("a2_g2_assumption", "unspecified"),
    }
    if message:
        payload["message"] = message
    if herwig is not None:
        payload["herwig"] = dict(herwig)
    return payload


def render_campaign_monitor_text(payload: Mapping[str, Any]) -> str:
    pdfs = payload.get("pdfs", {})
    lines = [
        f"Measurement: {payload['measurement']}",
        f"Tag: {payload['tag']}",
        f"Phase: {payload['phase']}",
        f"Elapsed: {fmt_seconds_compact(payload.get('elapsed_s'))}",
        f"Unpolarized PDF: {pdfs.get('unpolarized', 'unspecified')}",
        f"Polarized diff PDF: {pdfs.get('polarized_diff', 'unspecified')}",
        f"Prediction assumption: {payload.get('prediction_assumption', 'unspecified')}",
    ]
    message = payload.get("message")
    if isinstance(message, str) and message:
        lines.append(f"Message: {message}")

    herwig = payload.get("herwig")
    if isinstance(herwig, Mapping):
        lines.extend(
            [
                "",
                "Herwig",
                "------",
                f"Shards: completed {herwig.get('completed', 0)}/{herwig.get('total', 0)} | "
                f"running {herwig.get('running', 0)} | pending {herwig.get('pending', 0)} | "
                f"failed {herwig.get('failed', 0)}",
            ]
        )
        logical_rows = herwig.get("logical_rows")
        if isinstance(logical_rows, list) and logical_rows:
            lines.append("Logical runs with work remaining:")
            lines.extend(
                render_table(
                    ["Run", "Running", "Pending", "Total"],
                    logical_rows,
                    aligns=("l", "r", "r", "r"),
                )
            )
        active_rows = herwig.get("active_rows")
        if isinstance(active_rows, list) and active_rows:
            lines.append("Active shards:")
            lines.extend(
                render_table(
                    ["Run", "Tag", "Progress", "Events", "Seed", "Runtime"],
                    active_rows,
                    aligns=("l", "l", "r", "r", "r", "r"),
                )
            )
    return "\n".join(lines).rstrip() + "\n"


def write_campaign_monitor_files(campaign_dir: Path, payload: Mapping[str, Any]) -> None:
    atomic_write_json(campaign_status_json_path(campaign_dir), payload)
    atomic_write_text(campaign_status_txt_path(campaign_dir), render_campaign_monitor_text(payload))


def existing_herwig_monitor(campaign_dir: Path) -> Mapping[str, Any] | None:
    path = campaign_status_json_path(campaign_dir)
    if not path.is_file():
        return None
    try:
        payload = load_json(path)
    except CampaignError:
        return None
    herwig = payload.get("herwig")
    return herwig if isinstance(herwig, Mapping) else None


def emit_progress(lines: Sequence[str], interactive: bool) -> None:
    text = "\n".join(lines)
    if interactive:
        sys.stdout.write("\x1b[2J\x1b[H")
        sys.stdout.write(text + "\n")
    else:
        sys.stdout.write(text + "\n\n")
    sys.stdout.flush()


# ---------------------------------------------------------------------------
# Reference extraction and depolarization model
# ---------------------------------------------------------------------------


def r1990(x: float, q2: float) -> float:
    """Whitlow R1990 fit B, with Q2 expressed in GeV^2."""
    if x <= 0.0 or q2 <= 0.04:
        raise ValueError("R1990 requires x > 0 and Q2 > 0.04 GeV2")
    scale = 0.125**2
    theta = 1.0 + 12.0 * q2 / (q2 + 1.0) * scale / (scale + x * x)
    return 0.0635 / math.log(q2 / 0.04) * theta + 0.5747 / q2 - 0.3534 / (q2 * q2 + 0.09)


def depolarization(x: float, y: float, q2: float, proton_mass: float = 0.9382720813) -> float:
    """Longitudinal virtual-photon depolarization factor D.

    This is the conventional factor in A_parallel = D (A1 + eta A2). The
    inverse-D estimator neglects eta*A2 when forming A1. The ordinary
    A_parallel estimator never invokes this conversion.
    """
    if not (0.0 < y < 1.0):
        raise ValueError("Depolarization requires 0 < y < 1")
    gamma2 = 4.0 * proton_mass * proton_mass * x * x / q2
    numerator = y * (2.0 - y) * (1.0 + 0.5 * gamma2 * y)
    denominator = y * y * (1.0 + gamma2) + 2.0 * (1.0 + r1990(x, q2)) * (
        1.0 - y - 0.25 * gamma2 * y * y
    )
    if denominator <= 0.0:
        raise ValueError("Non-positive depolarization denominator")
    return numerator / denominator


def r1998(x: float, q2: float) -> float:
    """E143 R1998 average of the published Ra, Rb, and Rc fits."""
    if x <= 0.0 or q2 <= 0.04:
        raise ValueError("R1998 requires x > 0 and Q2 > 0.04 GeV2")
    a = (0.0485, 0.5470, 2.0621, -0.3804, 0.5090, -0.0285)
    b = (0.0481, 0.6114, -0.3509, -0.4611, 0.7172, -0.0317)
    c = (0.0577, 0.4644, 1.8288, 12.3708, -43.1043, 41.7415)
    theta = 1.0 + 12.0 * q2 / (q2 + 1.0) * 0.125**2 / (0.125**2 + x * x)
    logarithm = math.log(q2 / 0.04)
    ra = a[0] / logarithm * theta
    ra += a[1] / (q2**4 + a[2] ** 4) ** 0.25 * (1.0 + a[3] * x + a[4] * x * x) * x ** a[5]
    rb = b[0] / logarithm * theta
    rb += (b[1] / q2 + b[2] / (q2 * q2 + 0.3**2)) * (
        1.0 + b[3] * x + b[4] * x * x
    ) * x ** b[5]
    q2_threshold = c[3] * x + c[4] * x * x + c[5] * x**3
    rc = c[0] / logarithm * theta
    rc += c[1] / math.sqrt((q2 - q2_threshold) ** 2 + c[2] ** 2)
    return (ra + rb + rc) / 3.0


def compass_eta(
    x: float,
    y: float,
    q2: float,
    nucleon_mass: float = 0.938918754,
    muon_mass: float = 0.1056583755,
) -> float:
    """COMPASS eta factor including the finite incoming-muon mass."""
    if x <= 0.0 or q2 <= 0.0 or not (0.0 < y < 1.0):
        raise ValueError("COMPASS eta requires x,Q2 > 0 and 0 < y < 1")
    gamma = 2.0 * nucleon_mass * x / math.sqrt(q2)
    mass_term = y * y * muon_mass * muon_mass / q2
    numerator = gamma * (1.0 - y - 0.25 * gamma * gamma * y * y - mass_term)
    denominator = (1.0 + 0.5 * gamma * gamma * y) * (1.0 - 0.5 * y) - mass_term
    if denominator == 0.0:
        raise ValueError("Zero COMPASS eta denominator")
    return numerator / denominator


def compass_depolarization(
    x: float,
    y: float,
    q2: float,
    nucleon_mass: float = 0.938918754,
    muon_mass: float = 0.1056583755,
) -> float:
    """COMPASS virtual-photon depolarization factor with R1998 and muon mass."""
    if x <= 0.0 or q2 <= 0.0 or not (0.0 < y < 1.0):
        raise ValueError("COMPASS depolarization requires x,Q2 > 0 and 0 < y < 1")
    gamma2 = 4.0 * nucleon_mass * nucleon_mass * x * x / q2
    muon2_over_q2 = muon_mass * muon_mass / q2
    numerator = y * (
        (1.0 + 0.5 * gamma2 * y) * (2.0 - y) - 2.0 * y * y * muon2_over_q2
    )
    denominator = y * y * (1.0 - 2.0 * muon2_over_q2) * (1.0 + gamma2)
    denominator += 2.0 * (1.0 + r1998(x, q2)) * (
        1.0 - y - 0.25 * gamma2 * y * y
    )
    if denominator <= 0.0:
        raise ValueError("Non-positive COMPASS depolarization denominator")
    return numerator / denominator


def parse_five_column_member(payload: bytes, expected_rows: int = 15) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    for raw_line in payload.decode("ascii", errors="strict").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("*"):
            continue
        fields = line.split()
        if len(fields) != 5:
            raise CampaignError(f"Unexpected five-column reference-data row: {raw_line!r}")
        x_mean, q2_mean, value, stat, systematic = map(float, fields)
        rows.append(
            {
                "x_mean": x_mean,
                "q2_mean": q2_mean,
                "value": value,
                "stat": stat,
                "systematic_combined": systematic,
            }
        )
    if len(rows) != expected_rows:
        raise CampaignError(f"Expected {expected_rows} reference rows, found {len(rows)}")
    return rows


def extract_tar_reference(
    payload: bytes, member: str, member_sha256: str, expected_rows: int = 15
) -> list[dict[str, float]]:
    try:
        with tarfile.open(fileobj=io.BytesIO(payload), mode="r:*") as archive:
            extracted = archive.extractfile(member)
            if extracted is None:
                raise CampaignError(f"Archive member {member!r} is not a regular file")
            member_payload = extracted.read()
    except (tarfile.TarError, KeyError, OSError) as exc:
        raise CampaignError(f"Could not extract {member!r} from reference archive: {exc}") from exc
    actual = sha256_bytes(member_payload)
    if actual != member_sha256:
        raise CampaignError(f"Checksum mismatch for {member}: expected {member_sha256}, got {actual}")
    return parse_five_column_member(member_payload, expected_rows=expected_rows)


def _hepdata_symmetric_error(cell: Mapping[str, Any], label: str) -> float:
    for error in cell.get("errors", []):
        if str(error.get("label")) == label and "symerror" in error:
            return float(error["symerror"])
    raise CampaignError(f"HEPData cell does not contain a symmetric {label!r} uncertainty")


def parse_hepdata_table_json(
    payload: bytes, reference: Mapping[str, Any]
) -> list[dict[str, Any]]:
    """Validate and normalize the official HEPData table-display JSON."""
    try:
        table = json.loads(payload.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise CampaignError(f"Could not decode HEPData table JSON: {exc}") from exc
    if not isinstance(table, Mapping):
        raise CampaignError("HEPData table payload is not a JSON object")
    expected_name = str(reference["table_name"])
    expected_doi = str(reference["table_doi"])
    if str(table.get("name")) != expected_name or str(table.get("doi")) != expected_doi:
        raise CampaignError(
            f"HEPData table identity mismatch: expected {expected_name} ({expected_doi})"
        )
    headers = [str(header.get("name")) for header in table.get("headers", [])]
    if headers != [str(value) for value in reference["expected_headers"]]:
        raise CampaignError(f"Unexpected HEPData headers: {headers}")
    qualifiers = table.get("qualifiers", {})
    for label, expected in reference.get("expected_qualifiers", {}).items():
        entries = qualifiers.get(label, []) if isinstance(qualifiers, Mapping) else []
        actual = [str(entry.get("value")) for entry in entries]
        if str(expected) not in actual:
            raise CampaignError(
                f"HEPData qualifier {label!r} mismatch: expected {expected!r}, found {actual}"
            )
    values = table.get("values", [])
    if not isinstance(values, list) or len(values) != int(reference["expected_rows"]):
        raise CampaignError(
            f"Expected {reference['expected_rows']} HEPData rows, found {len(values) if isinstance(values, list) else 0}"
        )
    columns = reference["columns"]
    overrides = {
        (int(item["row"]), str(item["field"])): item
        for item in reference.get("normalization_overrides", [])
    }
    rows: list[dict[str, Any]] = []
    for row_index, entry in enumerate(values, start=1):
        independent = entry.get("x", [])
        dependent = entry.get("y", [])
        try:
            x_cell = independent[int(columns["x"])]
            q2_cell = independent[int(columns["q2"])]
            a1_cell = next(
                cell for cell in dependent if int(cell.get("group", -1)) == int(columns["a1_group"])
            )
        except (IndexError, KeyError, StopIteration, TypeError, ValueError) as exc:
            raise CampaignError(f"Malformed HEPData row {row_index}") from exc
        x_low = float(x_cell["low"])
        x_high = float(x_cell["high"])
        for field in ("low", "high"):
            override = overrides.get((row_index, field))
            if override is None:
                continue
            raw_value = x_low if field == "low" else x_high
            if not math.isclose(raw_value, float(override["from"]), rel_tol=0.0, abs_tol=1.0e-15):
                raise CampaignError(
                    f"HEPData normalization override no longer matches row {row_index} {field}"
                )
            if field == "low":
                x_low = float(override["to"])
            else:
                x_high = float(override["to"])
        x_mean_source = "published" if "value" in x_cell else "bin midpoint"
        x_mean = float(x_cell["value"]) if "value" in x_cell else 0.5 * (x_low + x_high)
        row: dict[str, Any] = {
            "x_low": x_low,
            "x_high": x_high,
            "x_mean": x_mean,
            "x_mean_source": x_mean_source,
            "q2_mean": float(q2_cell["value"]),
            "value": float(a1_cell["value"]),
            "stat": _hepdata_symmetric_error(a1_cell, "stat"),
            "systematic_combined": _hepdata_symmetric_error(a1_cell, "sys"),
        }
        if "g1_group" in columns:
            try:
                g1_cell = next(
                    cell
                    for cell in dependent
                    if int(cell.get("group", -1)) == int(columns["g1_group"])
                )
            except (StopIteration, TypeError, ValueError) as exc:
                raise CampaignError(f"HEPData row {row_index} has no configured g1 column") from exc
            row.update(
                {
                    "g1_value": float(g1_cell["value"]),
                    "g1_stat": _hepdata_symmetric_error(g1_cell, "stat"),
                    "g1_systematic_combined": _hepdata_symmetric_error(g1_cell, "sys"),
                }
            )
        rows.append(row)
    return rows


def validate_reference_snapshot(snapshot: Mapping[str, Any], rows: Sequence[Mapping[str, float]]) -> None:
    edges = snapshot.get("bin_edges", [])
    points = snapshot.get("points", [])
    if not points or len(edges) != len(points) + 1 or len(rows) != len(points):
        raise CampaignError("A reference snapshot must contain one more edge than points")
    if any(float(edges[index]) >= float(edges[index + 1]) for index in range(len(points))):
        raise CampaignError("Reference x-bin edges are not strictly increasing")
    keys = ("x_mean", "q2_mean", "value", "stat", "systematic_combined")
    for index, (point, row) in enumerate(zip(points, rows), start=1):
        if int(point.get("bin", -1)) != index:
            raise CampaignError(f"Unexpected bin number at reference point {index}")
        for key in keys:
            if not math.isclose(float(point[key]), float(row[key]), rel_tol=0.0, abs_tol=5.0e-7):
                raise CampaignError(f"Snapshot/source mismatch in bin {index} for {key}")
        if "x_low" in row and not math.isclose(
            float(edges[index - 1]), float(row["x_low"]), rel_tol=0.0, abs_tol=5.0e-12
        ):
            raise CampaignError(f"Snapshot/source mismatch in bin {index} for x_low")
        if "x_high" in row and not math.isclose(
            float(edges[index]), float(row["x_high"]), rel_tol=0.0, abs_tol=5.0e-12
        ):
            raise CampaignError(f"Snapshot/source mismatch in bin {index} for x_high")
        for key in ("g1_value", "g1_stat", "g1_systematic_combined"):
            if key in row and key in point and not math.isclose(
                float(point[key]), float(row[key]), rel_tol=0.0, abs_tol=5.0e-7
            ):
                raise CampaignError(f"Snapshot/source mismatch in bin {index} for {key}")
        systematics = point.get("systematics", {})
        if not isinstance(systematics, Mapping) or not systematics:
            raise CampaignError(f"Missing systematic components in reference bin {index}")
        if any(not math.isfinite(float(value)) or float(value) < 0.0 for value in systematics.values()):
            raise CampaignError(f"Invalid systematic component in reference bin {index}")
        quadrature = math.sqrt(sum(float(value) ** 2 for value in systematics.values()))
        if not math.isclose(quadrature, float(point["systematic_combined"]), rel_tol=0.0, abs_tol=1.0e-4):
            raise CampaignError(f"Systematic quadrature mismatch in reference bin {index}")


def _import_yoda() -> Any:
    try:
        import yoda  # type: ignore
    except (ImportError, OSError) as exc:
        raise CampaignError(
            "YODA Python bindings are unavailable. Load the same Herwig/Rivet module used for the campaign."
        ) from exc
    return yoda


def _write_yoda_objects(yoda: Any, objects: Sequence[Any], destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.suffix == ".gz":
        with tempfile.TemporaryDirectory(prefix="hermes-ref-", dir=destination.parent) as temporary:
            plain = Path(temporary) / destination.with_suffix("").name
            yoda.write(list(objects), str(plain))
            compressed = Path(temporary) / destination.name
            with plain.open("rb") as source, compressed.open("wb") as target:
                with gzip.GzipFile(filename="", mode="wb", fileobj=target, mtime=0) as zipped:
                    shutil.copyfileobj(source, zipped)
            os.replace(compressed, destination)
        return
    yoda.write(list(objects), str(destination))


def _additional_reference_paths(measurement: Mapping[str, Any]) -> list[Path]:
    configured = measurement["reference"].get("additional_snapshots", [])
    if not isinstance(configured, list):
        raise CampaignError("reference.additional_snapshots must be a list")
    paths: list[Path] = []
    for item in configured:
        if (not isinstance(item, Mapping) or set(item) != {"path", "sha256"}
                or not isinstance(item["path"], str)
                or not re.fullmatch(r"[0-9a-f]{64}", str(item["sha256"]))):
            raise CampaignError("Each additional reference snapshot requires path and a SHA-256 checksum")
        path = resolve_dispol_path(item["path"])
        if path in paths:
            raise CampaignError("Additional reference snapshot paths must be unique")
        if not path.is_file() or sha256_file(path) != item["sha256"]:
            raise CampaignError(f"Additional reference snapshot checksum mismatch: {item['path']}")
        paths.append(path)
    return paths


def _reference_source_paths(snapshot: Mapping[str, Any]) -> list[Path]:
    configured = snapshot.get("source_files", [])
    if not isinstance(configured, list):
        raise CampaignError("Reference source_files must be a list")
    paths: list[Path] = []
    for item in configured:
        if (not isinstance(item, Mapping) or set(item) != {"path", "sha256"}
                or not isinstance(item["path"], str)
                or not re.fullmatch(r"[0-9a-f]{64}", str(item["sha256"]))):
            raise CampaignError("Each reference source file requires path and a SHA-256 checksum")
        path = resolve_dispol_path(item["path"])
        if path in paths:
            raise CampaignError("Reference source file paths must be unique")
        if not path.is_file() or sha256_file(path) != item["sha256"]:
            raise CampaignError(f"Reference source file checksum mismatch: {item['path']}")
        paths.append(path)
    return paths


def load_reference_snapshot(measurement: Mapping[str, Any]) -> dict[str, Any]:
    """Merge checksum-pinned reference panels without guessing their target or axis."""
    snapshot = load_json(resolve_dispol_path(measurement["reference"]["snapshot"]))
    additional = _additional_reference_paths(measurement)
    if not additional:
        return snapshot
    merged = copy.deepcopy(snapshot)
    if not isinstance(merged.get("datasets"), list):
        raise CampaignError("Supplementary reference snapshots require a primary datasets list")
    provenance: list[Any] = []
    for path in additional:
        supplement = load_json(path)
        _reference_source_paths(supplement)
        if (supplement.get("measurement") != measurement["id"]
                or not isinstance(supplement.get("datasets"), list)):
            raise CampaignError(f"Supplementary reference snapshot has a different measurement or no datasets: {path}")
        merged["datasets"].extend(copy.deepcopy(supplement["datasets"]))
        provenance.append(supplement.get("provenance", {}))
    paths = [str(dataset.get("rivet_path", "")) for dataset in merged["datasets"]]
    ids = [str(dataset.get("id", "")) for dataset in merged["datasets"]]
    if (not all(paths) or len(set(paths)) != len(paths)
            or not all(ids) or len(set(ids)) != len(ids)):
        raise CampaignError("Merged reference datasets require unique explicit IDs and Rivet paths")
    merged.setdefault("provenance", {})["additional_references"] = provenance
    return merged


def write_reference_yoda(snapshot: Mapping[str, Any], destination: Path) -> None:
    if "datasets" in snapshot:
        datasets = snapshot["datasets"]
        if isinstance(datasets, list) and datasets and all("bin_edges" in item for item in datasets):
            yoda = _import_yoda()
            objects = [
                _reference_estimate(yoda, {**snapshot, **dataset})
                for dataset in datasets
            ]
            _write_yoda_objects(yoda, objects, destination)
            return
        try:
            from phenomenology_reference_data import write_reference_yoda as write_multi_reference
            generated = write_multi_reference(str(snapshot["measurement"]), snapshot)
        except Exception as exc:
            raise CampaignError(f"Could not write multi-dataset reference YODA: {exc}") from exc
        if generated.resolve() != destination.resolve():
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(generated, destination)
        return
    yoda = _import_yoda()
    _write_yoda_objects(yoda, [_reference_estimate(yoda, snapshot)], destination)


def ensure_reference_yoda(
    measurement: Mapping[str, Any], campaign_dir: Path,
    snapshot: Mapping[str, Any] | None = None,
) -> Path:
    """Materialize supplementary references in the ignored campaign cache."""
    supplementary = bool(measurement["reference"].get("additional_snapshots"))
    destination = (
        campaign_dir / "reference" / Path(str(measurement["analysis"]["reference_yoda"])).name
        if supplementary else resolve_dispol_path(measurement["analysis"]["reference_yoda"])
    )
    if supplementary or not destination.is_file() or destination.stat().st_size == 0:
        write_reference_yoda(
            load_reference_snapshot(measurement) if snapshot is None else snapshot,
            destination,
        )
    return destination


def _reference_estimate(yoda: Any, snapshot: Mapping[str, Any]) -> Any:
    edges = [float(value) for value in snapshot["bin_edges"]]
    path = str(snapshot.get("rivet_path", f"/REF/{snapshot['measurement']}/d01-x01-y01"))
    estimate = yoda.BinnedEstimate1D(edges, path)
    title = str(snapshot.get("provenance", {}).get("hepdata_table", snapshot.get("id", snapshot["measurement"])))
    try:
        estimate.setTitle(title)
    except Exception:
        estimate.setAnnotation("Title", title)
    estimate.setAnnotation("IsRef", 1)
    estimate.setAnnotation("Observable", str(snapshot.get("observable", "")))
    estimate.setAnnotation("Selection", str(snapshot.get("selection", "")))
    if "target" in snapshot:
        estimate.setAnnotation("Target", str(snapshot["target"]))
    estimate.setAnnotation("PlotAxis", str(snapshot.get("plot_axis", "x_mean")))
    estimate.setAnnotation("PublishedXMeans", json.dumps([point["x_mean"] for point in snapshot["points"]]))
    estimate.setAnnotation("PublishedQ2MeansGeV2", json.dumps([point["q2_mean"] for point in snapshot["points"]]))
    for index, point in enumerate(snapshot["points"], start=1):
        bin_object = estimate.bin(index)
        bin_object.setVal(float(point["value"]))
        errors = {"stat": float(point["stat"])}
        errors.update({str(label): float(error) for label, error in point["systematics"].items()})
        for label, error in errors.items():
            bin_object.setErr(-error, error, label)
    return estimate


def _multi_tar_reference_rows(
    measurement: Mapping[str, Any], snapshot: Mapping[str, Any], payload: bytes
) -> dict[str, list[dict[str, float]]]:
    """Verify each pinned archive member against its corresponding snapshot."""
    if snapshot.get("schema_version") != 2 or snapshot.get("measurement") != measurement["id"]:
        raise CampaignError("Multi-dataset reference snapshot schema or measurement does not match")
    sources = measurement["reference"].get("datasets")
    datasets = snapshot.get("datasets")
    if not isinstance(sources, list) or not sources or not isinstance(datasets, list) or not datasets:
        raise CampaignError("Multi-dataset references require non-empty source and snapshot datasets")
    source_ids = [str(item.get("id", "")) for item in sources]
    dataset_ids = [str(item.get("id", "")) for item in datasets]
    if (not all(source_ids) or len(set(source_ids)) != len(source_ids)
            or len(set(dataset_ids)) != len(dataset_ids) or set(source_ids) != set(dataset_ids)):
        raise CampaignError("Multi-dataset reference source IDs must match unique snapshot dataset IDs")
    by_id = {str(item["id"]): item for item in datasets}
    result: dict[str, list[dict[str, float]]] = {}
    for source in sources:
        dataset_id = str(source["id"])
        rows = extract_tar_reference(
            payload, str(source["archive_member"]), str(source["member_sha256"]),
            expected_rows=int(source["expected_rows"]),
        )
        validate_reference_snapshot(by_id[dataset_id], rows)
        result[dataset_id] = rows
    return result


def fetch_reference_data(measurement: Mapping[str, Any]) -> Path:
    reference = measurement["reference"]
    if str(reference.get("format")) == "hepdata-table-json-multidataset":
        try:
            from phenomenology_reference_data import fetch_and_validate
            cache = CAMPAIGN_ROOT / "_data_cache" / str(measurement["id"])
            outputs = fetch_and_validate(str(measurement["id"]), cache)
        except Exception as exc:
            raise CampaignError(f"Could not validate multi-dataset reference: {exc}") from exc
        if not outputs:
            raise CampaignError("Multi-dataset reference fetch produced no cached source")
        return outputs[0]
    snapshot_path = resolve_dispol_path(reference["snapshot"])
    snapshot = load_json(snapshot_path)
    request = urllib.request.Request(
        str(reference["source_url"]),
        headers={"User-Agent": "HerwigPol-experimental-reference/1.0"},
    )
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            payload = response.read()
    except OSError as exc:
        raise CampaignError(f"Could not download reference archive: {exc}") from exc
    actual = sha256_bytes(payload)
    if actual != reference["source_sha256"]:
        raise CampaignError(f"Source checksum mismatch: expected {reference['source_sha256']}, got {actual}")
    reference_format = str(reference.get("format"))
    cache_name = "official-source.dat"
    if reference_format == "tar-five-column-multidataset":
        rows = _multi_tar_reference_rows(measurement, snapshot, payload)
    elif reference_format == "tar-five-column":
        rows = extract_tar_reference(
            payload,
            reference["archive_member"],
            reference["member_sha256"],
            expected_rows=int(reference["expected_rows"]),
        )
    elif reference_format == "hepdata-table-json":
        rows = parse_hepdata_table_json(payload, reference)
        cache_name = "official-source.json"
        raw_snapshot = resolve_dispol_path(str(reference["raw_snapshot"]))
        if sha256_file(raw_snapshot) != actual or raw_snapshot.read_bytes() != payload:
            raise CampaignError(
                f"Vendored raw HEPData snapshot {raw_snapshot} does not match the checksum-validated source"
            )
    else:
        raise CampaignError(f"Unsupported reference format {reference_format!r}")
    if reference_format != "tar-five-column-multidataset":
        validate_reference_snapshot(snapshot, rows)

    cache = CAMPAIGN_ROOT / "_data_cache" / str(measurement["id"])
    cache.mkdir(parents=True, exist_ok=True)
    archive_path = cache / cache_name
    archive_path.write_bytes(payload)
    atomic_write_json(
        cache / "validation.json",
        {
            "measurement": measurement["id"],
            "verified_at": utc_now(),
            "source_url": reference["source_url"],
            "source_sha256": actual,
            "snapshot_sha256": sha256_file(snapshot_path),
            "rows": rows,
        },
    )
    if reference.get("additional_snapshots"):
        ensure_reference_yoda(measurement, cache)
    else:
        write_reference_yoda(load_reference_snapshot(measurement), resolve_dispol_path(measurement["analysis"]["reference_yoda"]))
    return archive_path


# ---------------------------------------------------------------------------
# Job matrix, runtime preflight, and preparation
# ---------------------------------------------------------------------------


def split_events(total: int, shards: int) -> list[int]:
    if total <= 0 or shards <= 0:
        raise CampaignError("Event counts and shard counts must be positive")
    if shards > total:
        raise CampaignError("The number of shards cannot exceed the events in a logical job")
    quotient, remainder = divmod(total, shards)
    return [quotient + (1 if index < remainder else 0) for index in range(shards)]


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
    if not scales or any(scale not in {0.5, 1.0, 2.0} for scale in scales):
        raise CampaignError("Hard scales must be selected from 0.5, 1, and 2")
    return scales


def _variation_points(
    profile: str,
    polarized_members: str | None = None,
    unpolarized_members: str | None = None,
    scales: str | None = None,
) -> list[tuple[int, int, float]]:
    polarized = _member_selector(polarized_members, profile)
    unpolarized = _member_selector(unpolarized_members, profile)
    scale_values = _scale_selector(scales, profile)
    # The two PDF ensembles are independent uncertainty sources.  Varying
    # their Cartesian product would not implement the prescribed replica
    # combination and would multiply the campaign size needlessly.
    points = {(0, 0, 1.0)}
    points.update((member, 0, 1.0) for member in polarized)
    points.update((0, member, 1.0) for member in unpolarized)
    points.update((0, 0, scale) for scale in scale_values)
    return sorted(points, key=lambda item: (item[2] != 1.0, item[0] != 0,
                                             item[1] != 0, item))


def _scale_token(value: float) -> str:
    return {0.5: "0p5", 1.0: "1", 2.0: "2"}[float(value)]


def _comparison_profile(measurement: Mapping[str, Any]) -> Mapping[str, Any]:
    profile = measurement.get("campaign", {}).get("comparison_profile", {})
    return profile if isinstance(profile, Mapping) else {}


def campaign_family_specs(
    measurement: Mapping[str, Any], include_comparisons: bool
) -> dict[str, dict[str, Any]]:
    nominal = {
        "label": "NLO+PS (polarized; full spin)",
        "id_prefix": "",
        "card_directory": str(measurement["cards"]["directory"]),
        "materialized_directory": "",
        "stem_pattern": str(measurement["cards"]["stem_pattern"]),
        "components": _configured_components(measurement),
        "helicities": dict(measurement["cards"]["helicities"]),
        "orders": dict(measurement["cards"]["orders"]),
        "event_options": {"POSNLO": "posnlo_events", "NEGNLO": "negnlo_events"},
        "postprocess": "helicity_asymmetry",
        "perturbative_order": "NLO+PS",
        "real_emission_spin_density": "enabled",
        "shower_spin_correlations": "enabled",
        "beam_polarization": "physical PP, PM, MP, MM",
    }
    families = {"nominal": nominal}
    if include_comparisons:
        configured = _comparison_profile(measurement).get("families", {})
        if not isinstance(configured, Mapping) or not configured:
            raise CampaignError(
                f"{measurement['id']} does not define registry comparison families"
            )
        for family_id, family in configured.items():
            copied = copy.deepcopy(dict(family))
            copied.setdefault("components", _configured_components(measurement))
            families[str(family_id)] = copied
    return families


def build_job_matrix(
    measurement: Mapping[str, Any],
    posnlo_events: int,
    negnlo_events: int,
    shards: int,
    seed_base: int,
    include_comparisons: bool = False,
    lo_events: int | None = None,
    variation_points: Sequence[Sequence[float | int]] | None = None,
) -> list[dict[str, Any]]:
    event_options = {
        "posnlo_events": int(posnlo_events),
        "negnlo_events": int(negnlo_events),
        "lo_events": int(posnlo_events if lo_events is None else lo_events),
    }
    jobs: list[dict[str, Any]] = []
    slot = 0
    nominal_variations = [
        (int(point[0]), int(point[1]), float(point[2]))
        for point in (variation_points or [(0, 0, 1.0)])
    ]
    if not nominal_variations or (0, 0, 1.0) not in nominal_variations:
        raise CampaignError("The DIS variation matrix must contain the central PDF/scale point")
    for family_id, family in campaign_family_specs(measurement, include_comparisons).items():
        stem_pattern = str(family["stem_pattern"])
        card_directory = str(family["card_directory"])
        materialized_directory = str(family.get("materialized_directory", "")).strip("/")
        id_prefix = str(family.get("id_prefix", "")).strip("-")
        family_variations = nominal_variations if family_id == "nominal" else [(0, 0, 1.0)]
        for polarized_member, unpolarized_member, scale in family_variations:
            central_variation = (
                polarized_member == 0 and unpolarized_member == 0
                and math.isclose(scale, 1.0)
            )
            variation_suffix = (
                "" if central_variation else
                f"-p{polarized_member:03d}-u{unpolarized_member:03d}-mu{_scale_token(scale)}"
            )
            for component in family["components"]:
                for helicity in family["helicities"]:
                    for order in family["orders"]:
                        option_name = str(family["event_options"][order])
                        event_splits = split_events(event_options[option_name], shards)
                        base_stem = stem_pattern.format(
                            component=component, helicity=helicity, order=order
                        )
                        stem = f"{base_stem}{variation_suffix}"
                        logical_parts = [
                            part for part in
                            (id_prefix, str(component), str(helicity), str(order)) if part
                        ]
                        logical_id = "-".join(logical_parts) + variation_suffix
                        for shard_index, events in enumerate(event_splits):
                            job_id = f"{logical_id}-s{shard_index + 1:03d}-of-{shards:03d}"
                            jobs.append(
                                {
                                    "id": job_id,
                                    "family": family_id,
                                    "family_label": str(family["label"]),
                                    "postprocess": str(family["postprocess"]),
                                    "component": str(component),
                                    "component_label": str(family["components"][component].get("label", component or "target")),
                                    "helicity": str(helicity),
                                    "order": str(order),
                                    "polarized_pdf_member": polarized_member,
                                    "unpolarized_pdf_member": unpolarized_member,
                                    "scale": scale,
                                    "stem": stem,
                                    "base_stem": base_stem,
                                    "card_source": f"{card_directory}/{base_stem}.in",
                                    "card_input": (
                                        f"{materialized_directory}/{stem}.in"
                                        if materialized_directory
                                        else f"{stem}.in"
                                    ),
                                    "shard": shard_index + 1,
                                    "shards": shards,
                                    "events": events,
                                    "seed": seed_base + slot,
                                    "initial_seed": seed_base + slot,
                                    "attempt": 0,
                                    "status": "planned",
                                    "output_yoda": f"yoda/{job_id}.yoda",
                                    "run_file": f"runs/{stem}.run",
                                }
                            )
                            slot += 1
    job_ids = [str(job["id"]) for job in jobs]
    if len(job_ids) != len(set(job_ids)):
        raise CampaignError("Comparison registry produces duplicate campaign job ids")
    return jobs


def _command_output(command: Sequence[str]) -> str:
    try:
        completed = subprocess.run(command, check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    except (OSError, subprocess.CalledProcessError) as exc:
        raise CampaignError(f"Preflight command failed: {' '.join(command)}: {exc}") from exc
    return completed.stdout.strip()


def _find_runtime_library(prefix: Path, relative_directory: str, pattern: str) -> Path:
    matches = sorted((prefix / relative_directory).glob(pattern))
    regular = [path.resolve() for path in matches if path.is_file()]
    if not regular:
        raise CampaignError(f"Active runtime {prefix} does not contain {relative_directory}/{pattern}")
    return regular[0]


def preflight_runtime(measurement: Mapping[str, Any]) -> dict[str, Any]:
    tools: dict[str, str] = {}
    for name in (
        "Herwig", "rivet", "rivet-build", "rivet-mkhtml", "rivet-config",
        "lhapdf", "lhapdf-config",
    ):
        resolved = shutil.which(name)
        if not resolved:
            raise CampaignError(
                f"Required executable {name!r} is not active. Load the polarized Herwig/Rivet environment first."
            )
        tools[name] = str(Path(resolved).resolve())

    herwig_prefix = Path(tools["Herwig"]).parent.parent.resolve()
    hwmedis = _find_runtime_library(herwig_prefix, "lib/Herwig", "HwMEDIS*.so*")
    hwmehadron = _find_runtime_library(
        herwig_prefix, "lib/Herwig", "HwMEHadron*.so*"
    )
    hwshower = _find_runtime_library(
        herwig_prefix, "lib/Herwig", "HwShower*.so*"
    )
    fixed_target = _find_runtime_library(herwig_prefix, "lib/ThePEG", "FixedTargetLuminosity*.so*")
    herwig_repository = herwig_prefix / "share" / "Herwig" / "HerwigDefaults.rpo"
    physics = measurement.get("physics", {})
    ensembles = measurement.get("pdf_ensembles", {})
    try:
        pdfs = [
            physics.get("unpolarized_pdf")
            or ensembles["unpolarized"]["set"],
            physics.get("polarized_pdf")
            or ensembles["polarized"]["set"],
        ]
    except (KeyError, TypeError) as exc:
        raise CampaignError(
            "Measurement does not define its polarized and unpolarized PDF sets"
        ) from exc
    for pdf in pdfs:
        _command_output([tools["lhapdf"], "show", str(pdf)])

    compiler = os.environ.get("CXX", "").strip()
    if compiler:
        compiler = shutil.which(compiler) or (compiler if Path(compiler).is_file() else "")
    if not compiler and sys.platform == "darwin":
        candidates = sorted(Path("/opt/homebrew/bin").glob("g++-[0-9]*"), reverse=True)
        compiler = str(candidates[0]) if candidates else ""

    try:
        file_provenance = provenance.runtime_record(
            repository=DISPOL_ROOT,
            tools=tools,
            herwig_prefix=herwig_prefix,
            artifact_paths={
                "Herwig": Path(tools["Herwig"]),
                "HerwigDefaults.rpo": herwig_repository,
                "HwMEDIS": hwmedis,
                "HwMEHadron": hwmehadron,
                "HwShower": hwshower,
                "FixedTargetLuminosity": fixed_target,
                "Rivet": Path(tools["rivet"]),
            },
            pdf_sets=[str(pdf) for pdf in pdfs],
            lhapdf_data_directory=Path(
                _command_output([tools["lhapdf-config"], "--datadir"])
            ),
        )
    except provenance.ProvenanceError as exc:
        raise CampaignError(str(exc)) from exc

    return {
        "checked_at": utc_now(),
        "tools": tools,
        "herwig_prefix": str(herwig_prefix),
        "herwig_version": _command_output([tools["Herwig"], "--version"]),
        "rivet_version": _command_output([tools["rivet"], "--version"]),
        "rivet_data_directory": _command_output([tools["rivet-config"], "--datadir"]),
        "hwmedis_library": str(hwmedis),
        "hwmehadron_library": str(hwmehadron),
        "hwshower_library": str(hwshower),
        "herwig_repository": str(herwig_repository.resolve()),
        "fixed_target_library": str(fixed_target),
        "rivet_plugin_compiler": compiler,
        "pdf_sets": pdfs,
        "provenance": file_provenance,
        "environment": {
            key: os.environ.get(key, "")
            for key in ("HERWIG_ENV", "RIVET_ANALYSIS_PATH", "RIVET_DATA_PATH", "DYLD_LIBRARY_PATH", "LD_LIBRARY_PATH")
        },
    }


def analysis_specs(measurement: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    """Return the primary Rivet analysis followed by optional companions.

    Descriptors without ``analysis.companions`` retain the historical
    single-analysis behaviour.  A companion is compiled into the same plugin
    and consumes the same generated event, but can have its own Rivet options
    and raw-object namespace.
    """

    primary = measurement["analysis"]
    companions = primary.get("companions", [])
    if not isinstance(companions, list) or any(
        not isinstance(spec, Mapping) for spec in companions
    ):
        raise CampaignError(
            f"{measurement['id']} analysis.companions must be a list of objects"
        )
    specs = [primary, *companions]
    names = [str(spec.get("name", "")) for spec in specs]
    if any(not name for name in names) or len(names) != len(set(names)):
        raise CampaignError(
            f"{measurement['id']} Rivet analysis names must be nonempty and unique"
        )
    for spec in specs:
        for key in ("name", "source", "info", "plot"):
            if not spec.get(key):
                raise CampaignError(
                    f"{measurement['id']} Rivet analysis {spec.get('name')!r} "
                    f"is missing {key}"
                )
    return specs


def analysis_environment(measurement: Mapping[str, Any], campaign_dir: Path, runtime: Mapping[str, Any]) -> dict[str, str]:
    environment = os.environ.copy()
    analysis_dirs = list(dict.fromkeys(
        str(resolve_dispol_path(spec["source"]).parent)
        for spec in analysis_specs(measurement)
    ))
    plugin_dir = campaign_dir / "build"
    old_analysis_path = environment.get("RIVET_ANALYSIS_PATH", "")
    old_data_path = environment.get("RIVET_DATA_PATH", "")
    environment["RIVET_ANALYSIS_PATH"] = os.pathsep.join(
        value for value in (str(plugin_dir), *analysis_dirs, old_analysis_path) if value
    )
    environment["RIVET_DATA_PATH"] = os.pathsep.join(
        value for value in (
            str(campaign_dir / "reference") if measurement["reference"].get("additional_snapshots") else "",
            *analysis_dirs, old_data_path,
        ) if value
    )
    if runtime.get("rivet_plugin_compiler"):
        environment["CXX"] = str(runtime["rivet_plugin_compiler"])
    return environment


def measurement_signature(
    measurement: Mapping[str, Any], *, plot_bytes: bytes | None = None
) -> str:
    registry_copy = copy.deepcopy(
        {key: value for key, value in measurement.items() if not key.startswith("_")}
    )
    campaign_config = registry_copy.get("campaign")
    if isinstance(campaign_config, dict):
        # Comparison families are optional extensions. Excluding them here
        # keeps pre-existing nominal manifests compatible; their own cards and
        # descriptor block are covered by comparison_signature().
        campaign_config.pop("comparison_profile", None)
    digest = hashlib.sha256(canonical_json_bytes(registry_copy))
    plot_path = resolve_dispol_path(measurement["analysis"]["plot"])
    files: list[Path] = []
    for spec in analysis_specs(measurement):
        files.extend(
            resolve_dispol_path(str(spec[key]))
            for key in ("source", "info", "plot")
        )
        files.extend(
            resolve_dispol_path(str(path))
            for path in spec.get("support_files", [])
        )
    files.append(resolve_dispol_path(measurement["reference"]["snapshot"]))
    additional_references = _additional_reference_paths(measurement)
    files.extend(additional_references)
    files.extend(dict.fromkeys(
        path for snapshot_path in additional_references
        for path in _reference_source_paths(load_json(snapshot_path))
    ))
    if measurement["reference"].get("raw_snapshot"):
        files.append(resolve_dispol_path(str(measurement["reference"]["raw_snapshot"])))
    card_dir = resolve_dispol_path(measurement["cards"]["directory"])
    files.extend(sorted(card_dir.glob("*.in")))
    for path in files:
        digest.update(str(path.relative_to(DISPOL_ROOT)).encode("utf-8"))
        digest.update(
            plot_bytes if plot_bytes is not None and path == plot_path
            else path.read_bytes()
        )
    return digest.hexdigest()


def comparison_signature(measurement: Mapping[str, Any]) -> str:
    profile = _comparison_profile(measurement)
    families = profile.get("families", {})
    if not isinstance(families, Mapping) or not families:
        raise CampaignError(f"{measurement['id']} does not define comparison families")
    digest = hashlib.sha256(canonical_json_bytes(profile))
    paths: set[Path] = set()
    for family in families.values():
        card_directory = resolve_dispol_path(str(family["card_directory"]))
        stem_pattern = str(family["stem_pattern"])
        components = family.get("components", _configured_components(measurement))
        for component in components:
            for helicity in family["helicities"]:
                for order in family["orders"]:
                    stem = stem_pattern.format(
                        component=component, helicity=helicity, order=order
                    )
                    path = (card_directory / f"{stem}.in").resolve()
                    if not path.is_file():
                        raise CampaignError(f"Missing comparison card {path}")
                    paths.add(path)
    for path in sorted(paths):
        digest.update(str(path.relative_to(DISPOL_ROOT)).encode("utf-8"))
        digest.update(path.read_bytes())
    return digest.hexdigest()


def _assert_manifest_signatures_current(
    manifest: Mapping[str, Any], measurement: Mapping[str, Any], *,
    allow_plot_metadata_refresh: bool = False,
) -> dict[str, Any] | None:
    """Refuse to reinterpret event products made with another definition."""

    configuration = manifest.get("configuration", {})
    recorded = configuration.get("measurement_signature")
    current = measurement_signature(measurement)
    refresh: dict[str, Any] | None = None
    if recorded != current:
        if allow_plot_metadata_refresh:
            refresh = authorize_plot_metadata_refresh(
                manifest,
                measurement,
                current_signature=current,
                signature_with_plot_bytes=lambda payload: measurement_signature(
                    measurement, plot_bytes=payload
                ),
            )
        else:
            raise CampaignError(
                f"{measurement['id']} campaign products were generated with a "
                "different analysis, reference, card, or metadata signature. "
                "Postprocessing and plotting are intentionally refused; use a "
                "new immutable tag and rerun the event campaign."
            )
    if bool(configuration.get("comparisons", False)):
        recorded_comparison = configuration.get("comparison_signature")
        current_comparison = comparison_signature(measurement)
        if recorded_comparison != current_comparison:
            raise CampaignError(
                f"{measurement['id']} comparison-family products have a "
                "stale signature. Use a new immutable tag and rerun them."
            )
    return refresh


def _git_file_at_commit(commit: str, relative_path: Path) -> bytes:
    if not re.fullmatch(r"[0-9a-f]{40}", commit):
        raise CampaignError(
            f"Cannot verify plot-only refresh against invalid source commit {commit!r}"
        )
    completed = subprocess.run(
        ["git", "-C", str(DISPOL_ROOT), "show", f"{commit}:{relative_path.as_posix()}"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if completed.returncode != 0:
        detail = completed.stderr.decode("utf-8", errors="replace").strip()
        raise CampaignError(
            f"Cannot read {relative_path} from generation commit {commit}: {detail}"
        )
    return completed.stdout


def authorize_plot_metadata_refresh(
    manifest: Mapping[str, Any],
    measurement: Mapping[str, Any],
    *,
    current_signature: str,
    signature_with_plot_bytes: Callable[[bytes], str],
) -> dict[str, Any]:
    """Prove that a stale campaign differs only in presentation metadata.

    The generation commit supplies the historical ``.plot`` bytes. Replacing
    only the current plot bytes with that historical payload must reproduce
    the immutable signature stored in the manifest. This keeps source,
    reference, cards, descriptor content, and postprocessed physics products
    under the original lock while permitting a recorded label/style refresh.
    """

    configuration = manifest.get("configuration", {})
    recorded = str(configuration.get("measurement_signature", ""))
    if manifest.get("status") != "complete":
        raise CampaignError(
            "Plot-metadata refresh is allowed only for a complete campaign"
        )
    jobs = manifest.get("jobs")
    if not isinstance(jobs, list) or not jobs or any(
        job.get("status") != "success" for job in jobs
    ):
        raise CampaignError(
            "Plot-metadata refresh requires a nonempty all-success job matrix"
        )
    source_control = (
        (manifest.get("runtime") or {}).get("provenance") or {}
    ).get("source_control") or {}
    generation_commit = str(source_control.get("commit", ""))
    plot_path = resolve_dispol_path(str(measurement["analysis"]["plot"]))
    relative_plot_path = plot_path.relative_to(DISPOL_ROOT)
    generation_plot_bytes = _git_file_at_commit(
        generation_commit, relative_plot_path
    )
    reconstructed = signature_with_plot_bytes(generation_plot_bytes)
    if reconstructed != recorded:
        raise CampaignError(
            f"{measurement['id']} differs from its generation signature in "
            "more than the Rivet .plot metadata; the presentation-only "
            "refresh is refused."
        )
    current_plot_bytes = plot_path.read_bytes()
    if current_plot_bytes == generation_plot_bytes:
        raise CampaignError(
            f"{measurement['id']} has unchanged .plot metadata, so its stale "
            "signature cannot be refreshed as presentation-only."
        )
    return {
        "mode": "presentation_only_plot_metadata_refresh",
        "generation_measurement_signature": recorded,
        "current_measurement_signature": current_signature,
        "generation_source_commit": generation_commit,
        "plot_path": relative_plot_path.as_posix(),
        "generation_plot_sha256": sha256_bytes(generation_plot_bytes),
        "current_plot_sha256": sha256_bytes(current_plot_bytes),
        "verified_at": utc_now(),
    }


def _resolved_campaign_options(args: argparse.Namespace, measurement: Mapping[str, Any]) -> dict[str, Any]:
    defaults = measurement["campaign"]
    smoke = bool(getattr(args, "smoke", False))
    profile = str(getattr(args, "profile", "central"))
    comparisons = bool(getattr(args, "comparisons", False)) or profile == "paper"
    smoke_events = int(defaults["smoke_events"])
    posnlo = smoke_events if smoke else int(getattr(args, "posnlo_events", None) or defaults["default_events"]["POSNLO"])
    negnlo = smoke_events if smoke else int(getattr(args, "negnlo_events", None) or defaults["default_events"]["NEGNLO"])
    lo = smoke_events if smoke else int(getattr(args, "lo_events", None) or posnlo)
    variation_points = _variation_points(
        profile,
        getattr(args, "polarized_pdf_members", None),
        getattr(args, "unpolarized_pdf_members", None),
        getattr(args, "scales", None),
    )
    return {
        "profile": profile,
        "jobs": int(getattr(args, "jobs", None) or defaults["default_jobs"]),
        "shards": int(getattr(args, "shards", None) or defaults["default_shards"]),
        "seed_base": int(getattr(args, "seed_base", None) or defaults["default_seed_base"]),
        "posnlo_events": posnlo,
        "negnlo_events": negnlo,
        "lo_events": lo,
        "smoke": smoke,
        "comparisons": comparisons,
        "variation_points": [list(point) for point in variation_points],
    }


def _plan_payload(measurement: Mapping[str, Any], tag: str, options: Mapping[str, Any]) -> dict[str, Any]:
    jobs = build_job_matrix(
        measurement,
        int(options["posnlo_events"]),
        int(options["negnlo_events"]),
        int(options["shards"]),
        int(options["seed_base"]),
        include_comparisons=bool(options["comparisons"]),
        lo_events=int(options["lo_events"]),
        variation_points=options.get("variation_points"),
    )
    logical_jobs = len({
        (
            str(job.get("family", "nominal")),
            str(job.get("component", "")),
            str(job["helicity"]),
            str(job["order"]),
            int(job.get("polarized_pdf_member", 0)),
            int(job.get("unpolarized_pdf_member", 0)),
            float(job.get("scale", 1.0)),
        )
        for job in jobs
    })
    return {
        "measurement": measurement["id"],
        "tag": tag,
        "campaign_directory": str(campaign_directory(str(measurement["id"]), tag)),
        "options": dict(options),
        "logical_jobs": logical_jobs,
        "shard_jobs": len(jobs),
        "jobs": jobs,
    }


def _copy_if_consistent(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        if destination.read_bytes() != source.read_bytes():
            raise CampaignError(f"Refusing to overwrite changed materialized input {destination}")
        return
    shutil.copy2(source, destination)


def _dis_variation_card_text(
    source: Path,
    common_source: Path,
    job: Mapping[str, Any],
) -> str:
    """Materialize a thin DIS card with coherent PDF and hard-scale axes."""
    text = source.read_text(encoding="utf-8")
    polarized_member = int(job.get("polarized_pdf_member", 0))
    unpolarized_member = int(job.get("unpolarized_pdf_member", 0))
    scale = float(job.get("scale", 1.0))
    base_stem = str(job.get("base_stem", job["stem"]))
    if (
        polarized_member == 0
        and unpolarized_member == 0
        and math.isclose(scale, 1.0)
        and str(job["stem"]) == base_stem
    ):
        return text

    common = common_source.read_text(encoding="utf-8")
    unpolarized_match = re.search(
        r"EPPolarizedExtractor:SecondPDF\s+(\S+)", common
    )
    polarized_match = re.search(
        r"EPPolarizedExtractor:SecondLongitudinalDifferencePDF\s+(\S+)", common
    )
    if unpolarized_match is None or polarized_match is None:
        raise CampaignError(
            f"Could not identify the polarized PDF objects in {common_source}"
        )
    overrides = [
        f"set {unpolarized_match.group(1)}:Member {unpolarized_member}",
        f"set /Herwig/Partons/HardLOPDF:Member {unpolarized_member}",
        f"set /Herwig/Partons/HardNLOPDF:Member {unpolarized_member}",
        f"set /Herwig/Partons/ShowerLOPDF:Member {unpolarized_member}",
        f"set /Herwig/Partons/ShowerNLOPDF:Member {unpolarized_member}",
        f"set {polarized_match.group(1)}:Member {polarized_member}",
    ]
    if not math.isclose(scale, 1.0):
        scale_objects = sorted(
            set(
                re.findall(
                    r"^set\s+(\S+):(?:MinimumScale|UseNativeDISWindowGeneration)\s+",
                    common,
                    re.MULTILINE,
                )
            )
        )
        if not scale_objects:
            raise CampaignError(f"No DIS scale objects were found in {common_source}")
        overrides.extend(
            f"set {object_path}:ScaleFactor {scale:.8g}"
            for object_path in scale_objects
        )
    saverun = f"saverun {job['stem']} EventGenerator"
    generated, replacements = re.subn(
        r"saverun\s+\S+\s+EventGenerator",
        "\n".join(overrides + [saverun]),
        text,
    )
    if replacements != 1:
        raise CampaignError(f"Expected one saverun command in {source}")
    return generated


def _run_logged(
    command: Sequence[str],
    cwd: Path,
    environment: Mapping[str, str],
    log_path: Path,
) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as log:
        log.write(f"command: {' '.join(command)}\n")
        log.write(f"cwd: {cwd}\n\n")
        log.flush()
        try:
            completed = subprocess.run(
                list(command), cwd=cwd, env=dict(environment), stdout=log, stderr=subprocess.STDOUT
            )
        except OSError as exc:
            raise CampaignError(f"Could not execute {command[0]}: {exc}") from exc
    if completed.returncode != 0:
        raise CampaignError(f"Command failed with status {completed.returncode}; see {log_path}")


def build_rivet_plugin(
    measurement: Mapping[str, Any],
    campaign_dir: Path,
    runtime: Mapping[str, Any],
) -> Path:
    build_dir = campaign_dir / "build"
    build_dir.mkdir(parents=True, exist_ok=True)
    plugin = build_dir / str(measurement["analysis"]["plugin"])
    specs = analysis_specs(measurement)
    sources = [resolve_dispol_path(str(spec["source"])) for spec in specs]
    include_dirs = list(dict.fromkeys(str(source.parent) for source in sources))
    environment = analysis_environment(measurement, campaign_dir, runtime)
    _run_logged(
        [
            runtime["tools"]["rivet-build"],
            str(plugin),
            *(str(source) for source in sources),
            *(f"-I{directory}" for directory in include_dirs),
        ],
        build_dir,
        environment,
        campaign_dir / "logs" / "rivet-build.log",
    )
    if not plugin.is_file() or plugin.stat().st_size == 0:
        raise CampaignError(f"rivet-build did not create {plugin}")
    for spec in specs:
        log_name = (
            "rivet-analysis-preflight.log"
            if len(specs) == 1
            else f"rivet-analysis-preflight-{spec['name']}.log"
        )
        _run_logged(
            [runtime["tools"]["rivet"], "--show-analysis", str(spec["name"])],
            build_dir,
            environment,
            campaign_dir / "logs" / log_name,
        )
    return plugin


def _manifest_configuration(
    measurement: Mapping[str, Any], tag: str, options: Mapping[str, Any], signature: str
) -> dict[str, Any]:
    configuration = {
        "measurement": measurement["id"],
        "tag": tag,
        "measurement_signature": signature,
        "jobs": int(options["jobs"]),
        "shards": int(options["shards"]),
        "seed_base": int(options["seed_base"]),
        "posnlo_events": int(options["posnlo_events"]),
        "negnlo_events": int(options["negnlo_events"]),
        "smoke": bool(options["smoke"]),
    }
    if bool(options.get("comparisons", False)):
        configuration.update(
            {
                "comparisons": True,
                "lo_events": int(options["lo_events"]),
                "comparison_signature": comparison_signature(measurement),
            }
        )
    profile = str(options.get("profile", "central"))
    variation_points = options.get("variation_points", [[0, 0, 1.0]])
    if profile != "central" or variation_points != [[0, 0, 1.0]]:
        configuration.update(
            {
                "profile": profile,
                "variation_points": variation_points,
                "pdf_uncertainty": (
                    "independent NNPDFpol2.0 and NNPDF4.0 replica variances, "
                    "combined in quadrature"
                ),
                "hard_scale_factors": [0.5, 1.0, 2.0],
            }
        )
    return configuration


def _assert_manifest_compatible(manifest: Mapping[str, Any], expected: Mapping[str, Any]) -> None:
    actual = manifest.get("configuration")
    if actual != expected:
        raise CampaignError(
            "An incompatible manifest already exists. Use a new tag rather than overwriting an existing campaign."
        )


def prepare_campaign(args: argparse.Namespace, measurement: Mapping[str, Any]) -> Path:
    options = _resolved_campaign_options(args, measurement)
    plan = _plan_payload(measurement, args.tag, options)
    campaign_dir = campaign_directory(str(measurement["id"]), args.tag)
    tracker_started = float(getattr(args, "_tracker_started_at", time.time()))
    if args.dry_run:
        print(json.dumps(plan, indent=2, sort_keys=True))
        return campaign_dir

    runtime = preflight_runtime(measurement)
    signature = measurement_signature(measurement)
    configuration = _manifest_configuration(measurement, args.tag, options, signature)
    manifest_path = campaign_dir / MANIFEST_NAME
    existing: dict[str, Any] | None = None
    if manifest_path.exists():
        existing = load_json(manifest_path)
        _assert_manifest_compatible(existing, configuration)

    for directory in ("build", "cards", "runs", "yoda", "logs", "work", "postprocess", "plots"):
        (campaign_dir / directory).mkdir(parents=True, exist_ok=True)
    write_campaign_monitor_files(
        campaign_dir,
        build_campaign_monitor_payload(
            measurement=measurement,
            tag=args.tag,
            phase="preparing",
            started_at=tracker_started,
            message="Building the Rivet plugin and materializing Herwig run files.",
        ),
    )

    snapshot = load_reference_snapshot(measurement)
    ensure_reference_yoda(measurement, campaign_dir, snapshot)

    source_card_dir = resolve_dispol_path(measurement["cards"]["directory"])
    materialized = campaign_dir / "cards"
    common_name = str(measurement["cards"]["common"])
    _copy_if_consistent(source_card_dir / common_name, materialized / common_name)
    stems = sorted({str(job["stem"]) for job in plan["jobs"]})
    cards_by_stem: dict[str, tuple[Path, Path, Mapping[str, Any]]] = {}
    for job in plan["jobs"]:
        stem = str(job["stem"])
        base_stem = str(job.get("base_stem", stem))
        default_source = f"{measurement['cards']['directory']}/{base_stem}.in"
        source = resolve_dispol_path(str(job.get("card_source", default_source)))
        destination = materialized / str(job.get("card_input", f"{stem}.in"))
        previous = cards_by_stem.setdefault(stem, (source, destination, job))
        if previous[0:2] != (source, destination):
            raise CampaignError(f"Logical card stem {stem} resolves to multiple sources")
    common_source = source_card_dir / common_name
    for _, (source, destination, job) in sorted(cards_by_stem.items()):
        generated = _dis_variation_card_text(source, common_source, job)
        if destination.exists():
            if destination.read_text(encoding="utf-8") != generated:
                raise CampaignError(
                    f"Refusing to overwrite changed materialized input {destination}"
                )
        else:
            atomic_write_text(destination, generated)

    plugin = build_rivet_plugin(measurement, campaign_dir, runtime)
    environment = analysis_environment(measurement, campaign_dir, runtime)
    manifest = existing or {
        "manifest_version": 1,
        "measurement": measurement["id"],
        "tag": args.tag,
        "created_at": utc_now(),
        "configuration": configuration,
        "jobs": plan["jobs"],
        "history": [],
    }
    manifest["status"] = "preparing"
    manifest["updated_at"] = utc_now()
    manifest["runtime"] = runtime
    manifest["plugin"] = str(plugin.relative_to(campaign_dir))
    manifest["plugin_provenance"] = provenance.file_record(plugin)
    manifest["reference_snapshot"] = {
        "path": measurement["reference"]["snapshot"],
        "sha256": sha256_file(resolve_dispol_path(measurement["reference"]["snapshot"])),
    }
    if measurement["reference"].get("additional_snapshots"):
        manifest["reference_snapshot"]["additional_snapshots"] = copy.deepcopy(
            measurement["reference"]["additional_snapshots"]
        )
    manifest["history"].append({"at": utc_now(), "action": "prepare"})
    atomic_write_json(manifest_path, manifest)

    card_inputs_by_stem = {
        str(job["stem"]): str(job.get("card_input", f"{job['stem']}.in"))
        for job in plan["jobs"]
    }
    for stem in stems:
        run_destination = campaign_dir / "runs" / f"{stem}.run"
        if run_destination.is_file() and run_destination.stat().st_size > 0:
            continue
        thin_name = card_inputs_by_stem[stem]
        _run_logged(
            [runtime["tools"]["Herwig"], "read", thin_name],
            materialized,
            environment,
            campaign_dir / "logs" / f"read-{stem}.log",
        )
        generated = materialized / f"{stem}.run"
        if not generated.is_file() or generated.stat().st_size == 0:
            raise CampaignError(f"Herwig read did not create {generated}")
        os.replace(generated, run_destination)

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
    manifest["updated_at"] = utc_now()
    manifest["history"].append({"at": utc_now(), "action": "prepared"})
    atomic_write_json(manifest_path, manifest)
    atomic_write_json(
        campaign_dir / "resolved-measurement.json",
        {key: value for key, value in measurement.items() if not key.startswith("_")},
    )
    write_campaign_monitor_files(
        campaign_dir,
        build_campaign_monitor_payload(
            measurement=measurement,
            tag=args.tag,
            phase="prepared",
            started_at=tracker_started,
            message=f"Prepared {len(plan['jobs'])} shard job(s).",
        ),
    )
    print(f"Prepared {measurement['id']} campaign at {campaign_dir}")
    return campaign_dir


def _nonempty(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0


def _reconcile_manifest_outputs(manifest: dict[str, Any], campaign_dir: Path) -> None:
    for job in manifest["jobs"]:
        output = campaign_dir / job["output_yoda"]
        if _nonempty(output):
            job["status"] = "success"
            job["output_size"] = output.stat().st_size
        elif job.get("status") in {"success", "running", "queued"}:
            job["status"] = "planned"


def pending_jobs(
    manifest: dict[str, Any],
    campaign_dir: Path,
    recover_failed: bool,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    _reconcile_manifest_outputs(manifest, campaign_dir)
    scheduled: list[dict[str, Any]] = []
    blocked_failed: list[dict[str, Any]] = []
    next_seed = max(int(job.get("seed", 0)) for job in manifest["jobs"]) + 1
    for job in manifest["jobs"]:
        status = job.get("status", "planned")
        if status == "success":
            continue
        if status == "failed":
            if not recover_failed:
                blocked_failed.append(job)
                continue
            old_seed = int(job["seed"])
            job["seed"] = next_seed
            next_seed += 1
            job["attempt"] = int(job.get("attempt", 0)) + 1
            job["recovered_from_seed"] = old_seed
        job["status"] = "queued"
        scheduled.append(job)
    return scheduled, blocked_failed


@dataclass
class JobResult:
    job_id: str
    success: bool
    returncode: int
    output_size: int
    message: str


@dataclass
class JobActivity:
    job: Mapping[str, Any]
    tag: str
    log_path: Path
    started_at: float


def _job_attempt_tag(
    job: Mapping[str, Any], measurement: Mapping[str, Any], campaign_tag: str
) -> str:
    """Return the compact tag passed to ``Herwig run -t``.

    Herwig constructs Rivet and EvtGen output names by appending this token to
    the already descriptive run-file stem.  Repeating the measurement,
    campaign, and complete job identity here can exceed the per-component
    filename limit (255 bytes on macOS) even when the surrounding path is
    valid.  The work directory and manifest already carry that full identity,
    so the Herwig-local token only needs to distinguish shard and attempt.
    """

    del measurement, campaign_tag
    return (
        f"s{int(job['shard']):03d}-of-{int(job['shards']):03d}-"
        f"a{int(job.get('attempt', 0)):02d}"
    )


def _job_log_path(job: Mapping[str, Any], campaign_dir: Path) -> Path:
    return campaign_dir / "logs" / (
        f"{job['id']}-a{int(job.get('attempt', 0)):02d}.log"
    )


def latest_progress_marker(log_path: Path) -> str:
    if not log_path.exists():
        return "-"
    try:
        with log_path.open("rb") as stream:
            stream.seek(0, os.SEEK_END)
            size = stream.tell()
            stream.seek(max(0, size - 65536))
            chunk = stream.read().decode("utf-8", errors="replace")
    except OSError:
        return "-"
    markers = list(PROGRESS_MARKER_RE.finditer(chunk.replace("\r", "\n")))
    if not markers:
        return "-"
    last = markers[-1]
    current = last.group("current")
    total = last.group("total") or last.group("total_alt")
    if current == "init":
        return "init"
    if not total:
        return current
    try:
        current_value = int(current)
        total_value = int(total)
    except ValueError:
        return f"{current}/{total}"
    if total_value <= 0:
        return f"{current_value}/{total_value}"
    return f"{current_value}/{total_value} ({100.0 * current_value / total_value:.1f}%)"


def _logical_group_counts(jobs: Sequence[Mapping[str, Any]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for job in jobs:
        label = str(job["stem"])
        counts[label] = counts.get(label, 0) + 1
    return counts


def build_herwig_progress_payload(
    all_jobs: Sequence[Mapping[str, Any]],
    pending: Sequence[Mapping[str, Any]],
    active: Sequence[JobActivity],
    max_listed: int,
    started_at: float,
) -> dict[str, Any]:
    completed = sum(1 for job in all_jobs if job.get("status") == "success")
    failed = sum(1 for job in all_jobs if job.get("status") == "failed")
    pending_counts = _logical_group_counts(pending)
    running_counts = _logical_group_counts([activity.job for activity in active])
    all_counts = _logical_group_counts(all_jobs)
    listing_limit = max(0, int(max_listed))

    logical_rows: list[list[str]] = []
    logical_names = sorted(set(pending_counts) | set(running_counts))
    for name in logical_names[:listing_limit]:
        logical_rows.append(
            [
                name,
                str(running_counts.get(name, 0)),
                str(pending_counts.get(name, 0)),
                str(all_counts.get(name, 0)),
            ]
        )

    active_rows: list[list[str]] = []
    for activity in active[:listing_limit]:
        job = activity.job
        active_rows.append(
            [
                str(job["stem"]),
                activity.tag,
                latest_progress_marker(activity.log_path),
                str(job["events"]),
                str(job["seed"]),
                f"{max(0.0, time.time() - activity.started_at):.0f}s",
            ]
        )

    return {
        "completed": completed,
        "running": len(active),
        "pending": len(pending),
        "failed": failed,
        "total": len(all_jobs),
        "elapsed_s": max(0.0, time.time() - started_at),
        "logical_rows": logical_rows,
        "active_rows": active_rows,
    }


def progress_lines(tag: str, herwig: Mapping[str, Any]) -> list[str]:
    lines = [
        f"Tag: {tag}",
        f"Elapsed: {float(herwig['elapsed_s']):.0f}s",
        f"Shards: completed {herwig['completed']}/{herwig['total']} | "
        f"running {herwig['running']} | pending {herwig['pending']} | "
        f"failed {herwig['failed']}",
    ]
    logical_rows = herwig.get("logical_rows")
    if isinstance(logical_rows, list) and logical_rows:
        lines.extend(["", "Logical runs with work remaining:"])
        lines.extend(
            render_table(
                ["Run", "Running", "Pending", "Total"],
                logical_rows,
                aligns=("l", "r", "r", "r"),
            )
        )
    active_rows = herwig.get("active_rows")
    if isinstance(active_rows, list) and active_rows:
        lines.extend(["", "Active shards:"])
        lines.extend(
            render_table(
                ["Run", "Tag", "Progress", "Events", "Seed", "Runtime"],
                active_rows,
                aligns=("l", "l", "r", "r", "r", "r"),
            )
        )
    return lines


def _run_campaign_job(
    job: Mapping[str, Any],
    measurement: Mapping[str, Any],
    campaign_dir: Path,
    runtime: Mapping[str, Any],
    campaign_tag: str,
) -> JobResult:
    job_id = str(job["id"])
    attempt = int(job.get("attempt", 0))
    run_file = campaign_dir / str(job["run_file"])
    output = campaign_dir / str(job["output_yoda"])
    work_dir = campaign_dir / "work" / job_id / f"attempt-{attempt:02d}"
    work_dir.mkdir(parents=True, exist_ok=True)
    tag = _job_attempt_tag(job, measurement, campaign_tag)
    output_stem = f"{job['stem']}-S{job['seed']}-{tag}"
    evtgen_sink = work_dir / f"{output_stem}-EvtGen.log"
    if not evtgen_sink.exists():
        try:
            evtgen_sink.symlink_to("/dev/null")
        except OSError:
            pass
    command = [
        runtime["tools"]["Herwig"],
        "run",
        str(run_file),
        "-N",
        str(job["events"]),
        "-s",
        str(job["seed"]),
        "-t",
        tag,
    ]
    environment = analysis_environment(measurement, campaign_dir, runtime)
    log_path = _job_log_path(job, campaign_dir)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as log:
        log.write(f"command: {' '.join(command)}\n")
        log.write(f"cwd: {work_dir}\n\n")
        log.flush()
        try:
            completed = subprocess.run(command, cwd=work_dir, env=environment, stdout=log, stderr=subprocess.STDOUT)
            returncode = completed.returncode
        except OSError as exc:
            return JobResult(job_id, False, 127, 0, str(exc))
    generated = work_dir / f"{output_stem}.yoda"
    if returncode != 0:
        return JobResult(job_id, False, returncode, 0, f"Herwig exited with {returncode}; see {log_path}")
    if not _nonempty(generated):
        return JobResult(job_id, False, returncode, 0, f"Missing or empty output {generated}")
    output.parent.mkdir(parents=True, exist_ok=True)
    os.replace(generated, output)
    return JobResult(job_id, True, returncode, output.stat().st_size, "ok")


def _run_tracked_campaign_job(
    job: Mapping[str, Any],
    measurement: Mapping[str, Any],
    campaign_dir: Path,
    runtime: Mapping[str, Any],
    campaign_tag: str,
    active: dict[str, JobActivity],
    active_lock: threading.Lock,
) -> JobResult:
    job_id = str(job["id"])
    activity = JobActivity(
        job=job,
        tag=_job_attempt_tag(job, measurement, campaign_tag),
        log_path=_job_log_path(job, campaign_dir),
        started_at=time.time(),
    )
    with active_lock:
        active[job_id] = activity
    try:
        return _run_campaign_job(job, measurement, campaign_dir, runtime, campaign_tag)
    finally:
        with active_lock:
            active.pop(job_id, None)


def run_campaign(args: argparse.Namespace, measurement: Mapping[str, Any]) -> Path:
    options = _resolved_campaign_options(args, measurement)
    campaign_dir = campaign_directory(str(measurement["id"]), args.tag)
    tracker_started = float(getattr(args, "_tracker_started_at", time.time()))
    manifest_path = campaign_dir / MANIFEST_NAME
    if args.dry_run and not manifest_path.exists():
        print(json.dumps(_plan_payload(measurement, args.tag, options), indent=2, sort_keys=True))
        return campaign_dir
    if not manifest_path.exists():
        prepare_campaign(args, measurement)
    manifest = load_json(manifest_path)
    expected = _manifest_configuration(measurement, args.tag, options, measurement_signature(measurement))
    _assert_manifest_compatible(manifest, expected)
    runtime = manifest.get("runtime") or preflight_runtime(measurement)
    scheduled, blocked = pending_jobs(manifest, campaign_dir, bool(args.recover_failed))
    if blocked:
        identifiers = ", ".join(job["id"] for job in blocked)
        raise CampaignError(
            f"Failed shards require fresh-seed recovery: {identifiers}. Re-run with --recover-failed."
        )
    if args.dry_run:
        print(json.dumps({"scheduled": scheduled, "already_complete": len(manifest["jobs"]) - len(scheduled)}, indent=2))
        return campaign_dir
    if not scheduled:
        manifest["status"] = "complete"
        manifest["updated_at"] = utc_now()
        manifest["history"].append(
            {"at": utc_now(), "action": "campaign-already-complete"}
        )
        atomic_write_json(manifest_path, manifest)
        herwig = build_herwig_progress_payload(
            manifest["jobs"], [], [], int(getattr(args, "max_listed", 12)), tracker_started
        )
        write_campaign_monitor_files(
            campaign_dir,
            build_campaign_monitor_payload(
                measurement=measurement,
                tag=args.tag,
                phase="campaign-complete",
                started_at=tracker_started,
                message="All non-empty shard outputs were already present.",
                herwig=herwig,
            ),
        )
        print(f"Campaign {measurement['id']}/{args.tag} is already complete")
        return campaign_dir

    manifest["status"] = "running"
    manifest["updated_at"] = utc_now()
    manifest["history"].append({"at": utc_now(), "action": "campaign", "scheduled": len(scheduled)})
    atomic_write_json(manifest_path, manifest)
    by_id = {job["id"]: job for job in manifest["jobs"]}
    failures: list[JobResult] = []
    max_workers = max(1, int(options["jobs"]))
    max_listed = max(0, int(getattr(args, "max_listed", 12)))
    progress_interval = float(getattr(args, "progress_interval", 5.0))
    interactive = sys.stdout.isatty()
    active: dict[str, JobActivity] = {}
    active_lock = threading.Lock()

    def tracker_snapshot(
        remaining: set[concurrent.futures.Future[JobResult]],
        futures: Mapping[concurrent.futures.Future[JobResult], Mapping[str, Any]],
        phase: str = "running-herwig",
        message: str | None = None,
    ) -> dict[str, Any]:
        with active_lock:
            active_jobs = list(active.values())
        order = {str(job["id"]): index for index, job in enumerate(manifest["jobs"])}
        active_jobs.sort(key=lambda item: order.get(str(item.job["id"]), len(order)))
        active_ids = {str(item.job["id"]) for item in active_jobs}
        pending = [
            futures[future]
            for future in remaining
            if str(futures[future]["id"]) not in active_ids
        ]
        pending.sort(key=lambda job: order.get(str(job["id"]), len(order)))
        herwig = build_herwig_progress_payload(
            manifest["jobs"], pending, active_jobs, max_listed, tracker_started
        )
        payload = build_campaign_monitor_payload(
            measurement=measurement,
            tag=args.tag,
            phase=phase,
            started_at=tracker_started,
            message=message,
            herwig=herwig,
        )
        write_campaign_monitor_files(campaign_dir, payload)
        return payload

    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                _run_tracked_campaign_job,
                job,
                measurement,
                campaign_dir,
                runtime,
                args.tag,
                active,
                active_lock,
            ): job
            for job in scheduled
        }
        remaining = set(futures)
        initial_payload = tracker_snapshot(remaining, futures)
        emit_progress(progress_lines(args.tag, initial_payload["herwig"]), interactive)
        last_progress = time.monotonic()
        while remaining:
            if progress_interval < 0.0:
                timeout: float | None = None
            else:
                refresh_period = max(1.0, progress_interval)
                until_refresh = refresh_period - (time.monotonic() - last_progress)
                timeout = max(0.05, until_refresh)
            done, _ = concurrent.futures.wait(
                remaining,
                timeout=timeout,
                return_when=concurrent.futures.FIRST_COMPLETED,
            )
            for future in done:
                job_spec = futures[future]
                try:
                    result = future.result()
                except Exception as exc:
                    result = JobResult(
                        str(job_spec["id"]), False, 1, 0, f"Campaign worker failed: {exc}"
                    )
                remaining.remove(future)
                job = by_id[result.job_id]
                job["finished_at"] = utc_now()
                job["returncode"] = result.returncode
                job["message"] = result.message
                job["output_size"] = result.output_size
                job["status"] = "success" if result.success else "failed"
                if not result.success:
                    failures.append(result)
                manifest["updated_at"] = utc_now()
                atomic_write_json(manifest_path, manifest)

            refresh_due = progress_interval >= 0.0 and (
                time.monotonic() - last_progress >= max(1.0, progress_interval)
            )
            if refresh_due and remaining:
                payload = tracker_snapshot(remaining, futures)
                emit_progress(progress_lines(args.tag, payload["herwig"]), interactive)
                last_progress = time.monotonic()

    manifest["status"] = "failed" if failures else "complete"
    manifest["updated_at"] = utc_now()
    manifest["history"].append(
        {"at": utc_now(), "action": "campaign-finished", "failures": len(failures)}
    )
    atomic_write_json(manifest_path, manifest)
    final_phase = "campaign-failed" if failures else "campaign-complete"
    final_message = (
        f"{len(failures)} shard(s) failed; use fresh-seed recovery after inspection."
        if failures
        else f"Completed {len(scheduled)} scheduled shard job(s)."
    )
    final_payload = tracker_snapshot(set(), {}, phase=final_phase, message=final_message)
    emit_progress(progress_lines(args.tag, final_payload["herwig"]), interactive)
    if failures:
        raise CampaignError(
            f"{len(failures)} shard(s) failed. Inspect logs and resume with --recover-failed for fresh seeds."
        )
    print(f"Completed {len(scheduled)} shard jobs in {campaign_dir}")
    return campaign_dir


# ---------------------------------------------------------------------------
# Normalized-bin postprocessing
# ---------------------------------------------------------------------------


@dataclass
class BinSeries:
    edges: list[float]
    values: list[float]
    variances: list[float]

    def __post_init__(self) -> None:
        if len(self.edges) != len(self.values) + 1 or len(self.values) != len(self.variances):
            raise CampaignError("Inconsistent binned-series dimensions")


def _same_edges(left: Sequence[float], right: Sequence[float]) -> bool:
    return len(left) == len(right) and all(
        math.isclose(a, b, rel_tol=0.0, abs_tol=1.0e-12) for a, b in zip(left, right)
    )


def _histogram_path(analysis: str, name: str) -> str:
    return f"/{analysis}/{name}"


def _rivet_instance_equivalent(left: str, right: str) -> bool:
    """Compare Rivet instances while tolerating numeric canonicalization.

    Rivet writes typed numeric options back to YODA paths (for example ``5``
    becomes ``5.0``), while the exact card instance remains the string stored
    in the immutable campaign manifest.  Option names and non-numeric values
    must still agree exactly.
    """

    def split(instance: str) -> tuple[str, dict[str, str]]:
        fields = instance.split(":")
        options: dict[str, str] = {}
        for field in fields[1:]:
            if "=" not in field:
                return fields[0], {"": instance}
            key, value = field.split("=", 1)
            options[key] = value
        return fields[0], options

    left_name, left_options = split(left)
    right_name, right_options = split(right)
    if left_name != right_name or left_options.keys() != right_options.keys():
        return False
    for key in left_options:
        left_value = left_options[key]
        right_value = right_options[key]
        if left_value == right_value:
            continue
        try:
            if math.isclose(
                float(left_value), float(right_value),
                rel_tol=0.0, abs_tol=1.0e-12,
            ):
                continue
        except ValueError:
            pass
        return False
    return True


def _resolve_histogram_path(
    objects: Mapping[str, Any], object_path: str,
) -> str:
    """Resolve one YODA path against Rivet's canonicalized option spelling."""

    if object_path in objects:
        return object_path
    stripped = object_path.removeprefix("/")
    if "/" not in stripped:
        return object_path
    requested_instance, object_name = stripped.split("/", 1)
    candidates: list[str] = []
    for candidate in objects:
        candidate_stripped = str(candidate).removeprefix("/")
        if "/" not in candidate_stripped:
            continue
        candidate_instance, candidate_object = candidate_stripped.split("/", 1)
        if (
            candidate_object == object_name
            and _rivet_instance_equivalent(
                requested_instance, candidate_instance
            )
        ):
            candidates.append(str(candidate))
    if len(candidates) == 1:
        return candidates[0]
    if len(candidates) > 1:
        raise CampaignError(
            f"Ambiguous canonical Rivet path for {object_path}: "
            + ", ".join(sorted(candidates))
        )
    return object_path


def _histogram_series_from_objects(
    objects: Mapping[str, Any], yoda_path: Path, object_path: str
) -> BinSeries:
    """Extract one normalized Estimate1D from an already-read YODA file."""

    resolved_path = (
        object_path
        if object_path in objects
        else _resolve_histogram_path(objects, object_path)
    )
    obj = objects.get(resolved_path)
    if obj is None:
        raise CampaignError(f"Missing {object_path} in {yoda_path}")
    bins = list(obj.bins()) if hasattr(obj, "bins") else []
    if not bins or not hasattr(bins[0], "val"):
        raise CampaignError(
            f"{resolved_path} in {yoda_path} is not the normalized Rivet Estimate1D; raw accumulators are not accepted"
        )
    edges = [float(bins[0].xMin())] + [float(bin_object.xMax()) for bin_object in bins]
    values = [float(bin_object.val()) for bin_object in bins]
    variances: list[float] = []
    for bin_object in bins:
        try:
            error = float(bin_object.errAvg("stats"))
        except Exception:
            error = 0.0
        variances.append(error * error if math.isfinite(error) else 0.0)
    return BinSeries(edges, values, variances)


def read_histogram_series_many(
    yoda_path: Path, object_paths: Mapping[str, str]
) -> dict[str, BinSeries]:
    """Read several normalized histograms while opening a YODA file once."""

    if not object_paths:
        return {}
    yoda = _import_yoda()
    try:
        objects = yoda.read(str(yoda_path))
    except Exception as exc:
        raise CampaignError(f"Could not read {yoda_path}: {exc}") from exc
    return {
        label: _histogram_series_from_objects(objects, yoda_path, object_path)
        for label, object_path in object_paths.items()
    }


def read_histogram_series(yoda_path: Path, object_path: str) -> BinSeries:
    return read_histogram_series_many(
        yoda_path, {"histogram": object_path}
    )["histogram"]


def combine_shard_series(series: Sequence[BinSeries], event_counts: Sequence[int]) -> BinSeries:
    if not series or len(series) != len(event_counts):
        raise CampaignError("Shard series and event counts must be non-empty and aligned")
    edges = series[0].edges
    if any(not _same_edges(edges, item.edges) for item in series[1:]):
        raise CampaignError("Cannot combine shard histograms with different bin edges")
    total = float(sum(event_counts))
    if total <= 0.0:
        raise CampaignError("Cannot normalize shards with zero generated events")
    fractions = [float(count) / total for count in event_counts]
    values = [
        sum(fraction * item.values[index] for fraction, item in zip(fractions, series))
        for index in range(len(edges) - 1)
    ]
    variances = [
        sum(fraction * fraction * item.variances[index] for fraction, item in zip(fractions, series))
        for index in range(len(edges) - 1)
    ]
    return BinSeries(list(edges), values, variances)


def add_independent_series(series: Sequence[BinSeries]) -> BinSeries:
    if not series:
        raise CampaignError("Need at least one series to add")
    edges = series[0].edges
    if any(not _same_edges(edges, item.edges) for item in series[1:]):
        raise CampaignError("Cannot add histograms with different bin edges")
    return BinSeries(
        list(edges),
        [sum(item.values[index] for item in series) for index in range(len(edges) - 1)],
        [sum(item.variances[index] for item in series) for index in range(len(edges) - 1)],
    )


def linear_combine_series(series: Mapping[str, BinSeries], coefficients: Mapping[str, float]) -> BinSeries:
    if set(series) != set(coefficients):
        raise CampaignError("Series and linear-combination coefficients do not have identical labels")
    first = next(iter(series.values()))
    if any(not _same_edges(first.edges, item.edges) for item in series.values()):
        raise CampaignError("Cannot linearly combine different bin edges")
    values: list[float] = []
    variances: list[float] = []
    for index in range(len(first.values)):
        values.append(sum(float(coefficients[label]) * series[label].values[index] for label in coefficients))
        variances.append(
            sum(float(coefficients[label]) ** 2 * series[label].variances[index] for label in coefficients)
        )
    return BinSeries(list(first.edges), values, variances)


def ratio_with_covariance(
    numerator: float,
    numerator_variance: float,
    denominator: float,
    denominator_variance: float,
    covariance: float,
    epsilon: float = 1.0e-15,
) -> tuple[float | None, float | None]:
    if not all(math.isfinite(value) for value in (numerator, numerator_variance, denominator, denominator_variance, covariance)):
        return None, None
    if abs(denominator) <= epsilon:
        return None, None
    value = numerator / denominator
    variance = (
        numerator_variance / (denominator * denominator)
        + numerator * numerator * denominator_variance / denominator**4
        - 2.0 * numerator * covariance / denominator**3
    )
    if variance < 0.0 and abs(variance) < 1.0e-12 * max(1.0, numerator_variance, denominator_variance):
        variance = 0.0
    if variance < 0.0 or not math.isfinite(value):
        return None, None
    return value, math.sqrt(variance)


def parity_residual(
    first: float,
    first_variance: float,
    second: float,
    second_variance: float,
    epsilon: float = 1.0e-15,
) -> tuple[float | None, float | None]:
    denominator = first + second
    if not all(math.isfinite(value) for value in (first, second, first_variance, second_variance)):
        return None, None
    if abs(denominator) <= epsilon:
        return None, None
    value = (first - second) / denominator
    derivative_first = 2.0 * second / (denominator * denominator)
    derivative_second = -2.0 * first / (denominator * denominator)
    variance = derivative_first**2 * first_variance + derivative_second**2 * second_variance
    return value, math.sqrt(max(0.0, variance))


def _logical_job_groups(
    manifest: Mapping[str, Any],
) -> dict[tuple[str, int, int, float, str, str, str], list[Mapping[str, Any]]]:
    groups: dict[tuple[str, int, int, float, str, str, str], list[Mapping[str, Any]]] = {}
    for job in manifest["jobs"]:
        key = (
            str(job.get("family", "nominal")),
            int(job.get("polarized_pdf_member", 0)),
            int(job.get("unpolarized_pdf_member", 0)),
            float(job.get("scale", 1.0)),
            str(job.get("component", "")),
            str(job["helicity"]),
            str(job["order"]),
        )
        groups.setdefault(key, []).append(job)
    for jobs in groups.values():
        jobs.sort(key=lambda job: int(job["shard"]))
    return groups


def _require_complete_matrix(
    manifest: dict[str, Any], measurement: Mapping[str, Any], campaign_dir: Path
) -> dict[tuple[str, int, int, float, str, str, str], list[Mapping[str, Any]]]:
    _reconcile_manifest_outputs(manifest, campaign_dir)
    groups = _logical_job_groups(manifest)
    include_comparisons = bool(manifest.get("configuration", {}).get("comparisons", False))
    family_specs = campaign_family_specs(measurement, include_comparisons)
    configured_variations = [
        (int(point[0]), int(point[1]), float(point[2]))
        for point in manifest.get("configuration", {}).get(
            "variation_points", [[0, 0, 1.0]]
        )
    ]
    required = {
        (family_id, polarized, unpolarized, scale,
         str(component), str(helicity), str(order))
        for family_id, family in family_specs.items()
        for polarized, unpolarized, scale in (
            configured_variations if family_id == "nominal" else [(0, 0, 1.0)]
        )
        for component in family["components"]
        for helicity in family["helicities"]
        for order in family["orders"]
    }
    missing_components = sorted(required - set(groups))
    if missing_components:
        raise CampaignError(f"Manifest is missing family/helicity/order components: {missing_components}")
    unexpected_components = sorted(set(groups) - required)
    if unexpected_components:
        raise CampaignError(
            f"Manifest has unexpected family/helicity/order components: {unexpected_components}"
        )
    incomplete: list[str] = []
    for key in sorted(required):
        jobs = groups[key]
        expected_shards = int(manifest["configuration"]["shards"])
        if len(jobs) != expected_shards:
            incomplete.append(
                f"{'/'.join(str(value) for value in key if value != '')} "
                f"({len(jobs)}/{expected_shards} shards)"
            )
            continue
        for job in jobs:
            if job.get("status") != "success" or not _nonempty(campaign_dir / str(job["output_yoda"])):
                incomplete.append(str(job["id"]))
    if incomplete:
        raise CampaignError("Refusing postprocessing; incomplete components: " + ", ".join(incomplete))
    return groups


def _load_component_series(
    groups: Mapping[
        tuple[str, int, int, float, str, str, str],
        Sequence[Mapping[str, Any]],
    ],
    campaign_dir: Path,
    analysis: str,
    object_name: str,
    family_id: str = "nominal",
    component: str = "",
    helicities: Sequence[str] = ("PP", "PM", "MP", "MM"),
    orders: Sequence[str] = ("POSNLO", "NEGNLO"),
    polarized_pdf_member: int = 0,
    unpolarized_pdf_member: int = 0,
    scale: float = 1.0,
    order_coefficients: Mapping[str, float] | None = None,
) -> dict[str, BinSeries]:
    by_helicity_order: dict[tuple[str, str], BinSeries] = {}
    for helicity in helicities:
        for order in orders:
            jobs = groups[(family_id, int(polarized_pdf_member),
                           int(unpolarized_pdf_member), float(scale),
                           str(component), str(helicity), str(order))]
            shard_series = [
                read_histogram_series(
                    campaign_dir / str(job["output_yoda"]),
                    _histogram_path(analysis, object_name),
                )
                for job in jobs
            ]
            by_helicity_order[(str(helicity), str(order))] = combine_shard_series(
                shard_series, [int(job["events"]) for job in jobs]
            )
    return {
        helicity: linear_combine_series(
            {str(order): by_helicity_order[(str(helicity), str(order))] for order in orders},
            {str(order): float((order_coefficients or {}).get(str(order), 1.0)) for order in orders},
        )
        for helicity in helicities
    }


def _covariance_array(
    covariance_series: Mapping[str, BinSeries],
    numerator_coefficients: Mapping[str, float],
    denominator_coefficients: Mapping[str, float],
) -> list[float]:
    return [
        sum(
            float(numerator_coefficients[label])
            * float(denominator_coefficients[label])
            * covariance_series[label].variances[index]
            for label in numerator_coefficients
        )
        for index in range(len(next(iter(covariance_series.values())).values))
    ]


def _estimate_from_values(
    yoda: Any,
    edges: Sequence[float],
    path: str,
    values: Sequence[float | None],
    errors: Sequence[float | None],
    annotations: Mapping[str, Any],
) -> Any:
    estimate = yoda.BinnedEstimate1D([float(edge) for edge in edges], path)
    for key, value in annotations.items():
        estimate.setAnnotation(str(key), value)
    for index, (value, error) in enumerate(zip(values, errors), start=1):
        if value is None or error is None or not math.isfinite(value) or not math.isfinite(error):
            try:
                estimate.maskBin(index)
            except Exception:
                pass
            continue
        bin_object = estimate.bin(index)
        bin_object.setVal(float(value))
        bin_object.setErr(-float(error), float(error), "stat")
    return estimate


def _series_estimate(
    yoda: Any,
    series: BinSeries,
    path: str,
    annotations: Mapping[str, Any],
    density: bool = False,
) -> Any:
    values: list[float] = []
    errors: list[float] = []
    for index, (value, variance) in enumerate(zip(series.values, series.variances)):
        width = series.edges[index + 1] - series.edges[index]
        factor = 1.0 / width if density else 1.0
        values.append(value * factor)
        errors.append(math.sqrt(max(0.0, variance)) * factor)
    return _estimate_from_values(yoda, series.edges, path, values, errors, annotations)


def _asymmetry_outputs(
    ordinary: Mapping[str, BinSeries],
    weighted: Mapping[str, BinSeries] | None,
    covariance_proxy: Mapping[str, BinSeries] | None,
    unpolarized_coefficients: Mapping[str, float],
    longitudinal_coefficients: Mapping[str, float],
    parity_ordinary: Mapping[str, BinSeries] | None = None,
) -> dict[str, Any]:
    parity_source = ordinary if parity_ordinary is None else parity_ordinary
    sigma_uu = linear_combine_series(ordinary, unpolarized_coefficients)
    sigma_ll = linear_combine_series(ordinary, longitudinal_coefficients)
    if (weighted is None) != (covariance_proxy is None):
        raise CampaignError("The inverse-D estimator requires both weighted and covariance histograms")
    sigma_ll_over_d = (
        linear_combine_series(weighted, longitudinal_coefficients)
        if weighted is not None else None
    )
    covariance_parallel = _covariance_array(ordinary, longitudinal_coefficients, unpolarized_coefficients)
    covariance_a1 = (
        _covariance_array(covariance_proxy, longitudinal_coefficients, unpolarized_coefficients)
        if covariance_proxy is not None else None
    )
    apar_values: list[float | None] = []
    apar_errors: list[float | None] = []
    a1_values: list[float | None] = []
    a1_errors: list[float | None] = []
    parity_pp_mm_values: list[float | None] = []
    parity_pp_mm_errors: list[float | None] = []
    parity_pm_mp_values: list[float | None] = []
    parity_pm_mp_errors: list[float | None] = []
    for index in range(len(sigma_uu.values)):
        value, error = ratio_with_covariance(
            sigma_ll.values[index], sigma_ll.variances[index],
            sigma_uu.values[index], sigma_uu.variances[index], covariance_parallel[index]
        )
        apar_values.append(value)
        apar_errors.append(error)
        if sigma_ll_over_d is not None and covariance_a1 is not None:
            value, error = ratio_with_covariance(
                sigma_ll_over_d.values[index], sigma_ll_over_d.variances[index],
                sigma_uu.values[index], sigma_uu.variances[index], covariance_a1[index]
            )
            a1_values.append(value)
            a1_errors.append(error)
        value, error = parity_residual(
            parity_source["PP"].values[index], parity_source["PP"].variances[index],
            parity_source["MM"].values[index], parity_source["MM"].variances[index]
        )
        parity_pp_mm_values.append(value)
        parity_pp_mm_errors.append(error)
        value, error = parity_residual(
            parity_source["PM"].values[index], parity_source["PM"].variances[index],
            parity_source["MP"].values[index], parity_source["MP"].variances[index]
        )
        parity_pm_mp_values.append(value)
        parity_pm_mp_errors.append(error)
    result = {
        "edges": sigma_uu.edges,
        "sigma_uu": sigma_uu,
        "sigma_ll": sigma_ll,
        "apar_values": apar_values,
        "apar_errors": apar_errors,
        "parity_pp_mm_values": parity_pp_mm_values,
        "parity_pp_mm_errors": parity_pp_mm_errors,
        "parity_pm_mp_values": parity_pm_mp_values,
        "parity_pm_mp_errors": parity_pm_mp_errors,
    }
    if sigma_ll_over_d is not None:
        result.update(sigma_ll_over_d=sigma_ll_over_d, a1_values=a1_values, a1_errors=a1_errors)
    return result


def _selection_supported_bins(output: Mapping[str, Any], count: int) -> set[int]:
    bins = set(output.get("supported_bins", range(1, count + 1)))
    if any(type(index) is not int or not 1 <= index <= count for index in bins):
        raise CampaignError("Selection supported_bins are outside the histogram bin range")
    return bins


def _mask_selection_results(result: dict[str, Any], output: Mapping[str, Any]) -> None:
    supported = _selection_supported_bins(output, len(result["edges"]) - 1)
    result["supported_bins"] = sorted(supported)
    result["axis"] = str(output.get("axis", "x"))
    for key in ("apar_values", "apar_errors", "a1_values", "a1_errors",
                "parity_pp_mm_values", "parity_pp_mm_errors", "parity_pm_mp_values", "parity_pm_mp_errors"):
        if key in result:
            result[key] = [value if index in supported else None
                           for index, value in enumerate(result[key], start=1)]


def _mask_unsupported_estimate(estimate: Any, supported: set[int], count: int) -> Any:
    for index in range(1, count + 1):
        if index not in supported:
            estimate.maskBin(index)
    return estimate


def _summary_rows(selection: str, result: Mapping[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    edges = result["edges"]
    axis = str(result.get("axis", "x"))
    supported_bins = set(result.get("supported_bins", range(1, len(edges))))
    for index in range(len(edges) - 1):
        supported = index + 1 in supported_bins
        row = {
                "selection": selection,
                "bin": index + 1,
                "axis": axis,
                "bin_low": edges[index],
                "bin_high": edges[index + 1],
                f"{axis}_low": edges[index],
                f"{axis}_high": edges[index + 1],
                "supported": supported,
                "sigma_uu_pb": result["sigma_uu"].values[index] if supported else None,
                "sigma_uu_stat_pb": math.sqrt(max(0.0, result["sigma_uu"].variances[index])) if supported else None,
                "sigma_ll_pb": result["sigma_ll"].values[index] if supported else None,
                "sigma_ll_stat_pb": math.sqrt(max(0.0, result["sigma_ll"].variances[index])) if supported else None,
                "a_parallel": result["apar_values"][index],
                "a_parallel_stat": result["apar_errors"][index],
                "parity_pp_mm": result["parity_pp_mm_values"][index],
                "parity_pp_mm_stat": result["parity_pp_mm_errors"][index],
                "parity_pm_mp": result["parity_pm_mp_values"][index],
                "parity_pm_mp_stat": result["parity_pm_mp_errors"][index],
            }
        if "a1_values" in result:
            row.update(a1=result["a1_values"][index], a1_stat=result["a1_errors"][index])
        rows.append(row)
    return rows


def _family_annotations(
    measurement: Mapping[str, Any], family_id: str, family: Mapping[str, Any]
) -> dict[str, str]:
    physics = measurement.get("physics", {})
    orders = [str(order) for order in family["orders"]]
    order_coefficients = _order_combination_coefficients(measurement, orders)
    order_formula = _order_combination_formula(order_coefficients)
    target_uu = _target_component_coefficients(measurement, family, "unpolarized")
    target_ll = _target_component_coefficients(measurement, family, "longitudinal")
    active_components = [
        component for component in target_uu
        if component and (target_uu[component] != 0.0 or target_ll[component] != 0.0)
    ]
    annotations = {
        "Generator": "HerwigPol full shower+hadronization",
        "CampaignFamily": family_id,
        "CampaignFamilyLabel": str(family["label"]),
        "PerturbativeOrder": str(family.get("perturbative_order", ",".join(orders))),
        "OrderCombination": order_formula,
        "OrderCoefficients": json.dumps(order_coefficients, sort_keys=True),
        "OrderInputConvention": str(physics.get(
            "order_input_convention", "coefficients multiply stored normalized order bins"
        )),
        "BeamPolarization": str(family.get("beam_polarization", "unspecified")),
        "RealEmissionSpinDensity": str(
            family.get("real_emission_spin_density", "unspecified")
        ),
        "ShowerSpinCorrelations": str(
            family.get("shower_spin_correlations", "unspecified")
        ),
        "A2G2Assumption": str(physics.get("a2_g2_assumption", "unspecified")),
        "DepolarizationModel": str(physics.get("depolarization_model", "R1990")),
        "TargetComponents": ",".join(active_components) or "single target",
        "TargetModel": str(physics.get("target_model", "single physical target")),
    }
    if "d_state_factor" in physics:
        annotations["DeuteronDStateFactor"] = str(physics["d_state_factor"])
    if set(orders) == {"POSNLO", "NEGNLO"}:
        annotations["NLOCombination"] = f"normalized {order_formula} bins"
    if str(family["postprocess"]) == "helicity_asymmetry":
        annotations["HelicityCombination"] = "independent PP,PM,MP,MM samples"
    else:
        annotations["HelicityCombination"] = "direct 00 unpolarized sample"
    return annotations


def _target_component_coefficients(
    measurement: Mapping[str, Any], family: Mapping[str, Any], channel: str,
    selection: str | None = None,
) -> dict[str, float]:
    components = [str(value) for value in family["components"]]
    configured = measurement["combination"].get("target_components")
    if selection is not None:
        configured = measurement["outputs"][selection].get("target_components", configured)
    if configured is None:
        if len(components) != 1:
            raise CampaignError("Multiple target components require explicit combination coefficients")
        return {components[0]: 1.0}
    _validate_target_component_maps(configured, components, "Target-component combination")
    coefficients = configured[channel]
    return {str(key): float(value) for key, value in coefficients.items()}


def _selection_annotations(
    measurement: Mapping[str, Any], selection: str, annotations: Mapping[str, str]
) -> dict[str, str]:
    """Annotate the target represented by this output rather than all generated beams."""
    output = measurement["outputs"][selection]
    result = dict(annotations)
    result["PlotAxis"] = "q2_mean" if output.get("axis", "x") == "q2" else "x_mean"
    if output.get("estimator", "a1") == "a_parallel":
        result.update(
            Estimator="ordinary sigma_LL / ordinary sigma_UU",
            DepolarizationModel="not used by the direct A_parallel estimator",
            A2G2Assumption="no A_parallel-to-A1 conversion; generator physics retained",
        )
    for key, annotation in (("projection", "Projection"), ("reference_status", "ReferenceStatus"),
                            ("apar_integration_variable", "AParallelIntegrationVariable"),
                            ("apar_integration_range", "AParallelIntegrationRange"),
                            ("integration_variable", "IntegrationVariable"),
                            ("integration_range", "IntegrationRange"),
                            ("x_range", "XRange"), ("supported_bins", "SupportedBins")):
        if key in output:
            result[annotation] = (json.dumps(output[key]) if isinstance(output[key], (list, dict))
                                  else str(output[key]))
    for key, annotation in (("target", "Target"), ("target_model", "TargetModel")):
        if key in output:
            result[annotation] = str(output[key])
    if "d_state_factor" in output:
        if output["d_state_factor"] is None:
            result.pop("DeuteronDStateFactor", None)
        else:
            result["DeuteronDStateFactor"] = str(output["d_state_factor"])
    target = output.get("target_components", measurement["combination"].get("target_components"))
    if isinstance(target, Mapping):
        result["TargetComponentCombination"] = json.dumps(target, sort_keys=True)
        result["TargetComponents"] = ",".join(
            component for component in target["unpolarized"]
            if any(float(target[channel][component]) != 0.0 for channel in ("unpolarized", "longitudinal"))
        ) or "single target"
    return result


def _apar_annotations(annotations: Mapping[str, Any]) -> dict[str, Any]:
    result = dict(annotations)
    suffix = {"proton": "p", "deuteron": "d", "P": "p", "D": "d"}.get(str(result.get("Target", "")), "")
    result.update(
        Observable=f"A_parallel{suffix}",
        Estimator="ordinary sigma_LL / ordinary sigma_UU",
        DepolarizationCorrection="none",
        DepolarizationModel="not used by the direct A_parallel estimator",
        A2G2Assumption="no analysis-level A2 neglect; generated cross-section model",
    )
    return result


def _component_helicity_label(component: str, helicity: str) -> str:
    return f"{component}:{helicity}" if component else helicity


def _load_asymmetry_family_products(
    groups: Mapping[
        tuple[str, int, int, float, str, str, str],
        Sequence[Mapping[str, Any]],
    ],
    campaign_dir: Path,
    measurement: Mapping[str, Any],
    family_id: str,
    family: Mapping[str, Any],
    polarized_pdf_member: int = 0,
    unpolarized_pdf_member: int = 0,
    scale: float = 1.0,
) -> tuple[dict[str, dict[str, Any]], dict[str, BinSeries]]:
    analysis = str(measurement["analysis"]["name"])
    components = tuple(str(value) for value in family["components"])
    helicities = tuple(str(value) for value in family["helicities"])
    orders = tuple(str(value) for value in family["orders"])
    order_coefficients = _order_combination_coefficients(measurement, orders)
    helicity_uu_coefficients = {
        key: float(value) for key, value in measurement["combination"]["unpolarized"].items()
    }
    helicity_ll_coefficients = {
        key: float(value) for key, value in measurement["combination"]["longitudinal"].items()
    }
    results: dict[str, dict[str, Any]] = {}
    for selection, names in measurement["raw_observables"].items():
        direct_apar = measurement["outputs"][selection].get("estimator", "a1") == "a_parallel"
        target_uu_coefficients = _target_component_coefficients(
            measurement, family, "unpolarized", str(selection)
        )
        target_ll_coefficients = _target_component_coefficients(
            measurement, family, "longitudinal", str(selection)
        )
        uu_coefficients = {
            _component_helicity_label(component, helicity):
            target_uu_coefficients[component] * helicity_uu_coefficients[helicity]
            for component in components for helicity in helicities
        }
        ll_coefficients = {
            _component_helicity_label(component, helicity):
            target_ll_coefficients[component] * helicity_ll_coefficients[helicity]
            for component in components for helicity in helicities
        }
        ordinary: dict[str, BinSeries] = {}
        weighted: dict[str, BinSeries] = {}
        covariance: dict[str, BinSeries] = {}
        ordinary_by_component: dict[str, dict[str, BinSeries]] = {}
        for component in components:
            loaded_ordinary = _load_component_series(
                groups, campaign_dir, analysis, str(names["ordinary"]),
                family_id, component, helicities, orders,
                polarized_pdf_member, unpolarized_pdf_member, scale, order_coefficients,
            )
            loaded_weighted = {} if direct_apar else _load_component_series(
                groups, campaign_dir, analysis, str(names["weighted"]),
                family_id, component, helicities, orders,
                polarized_pdf_member, unpolarized_pdf_member, scale, order_coefficients,
            )
            loaded_covariance = {} if direct_apar else _load_component_series(
                groups, campaign_dir, analysis, str(names["covariance"]),
                family_id, component, helicities, orders,
                polarized_pdf_member, unpolarized_pdf_member, scale, order_coefficients,
            )
            ordinary_by_component[component] = loaded_ordinary
            for helicity in helicities:
                label = _component_helicity_label(component, helicity)
                ordinary[label] = loaded_ordinary[helicity]
                if not direct_apar:
                    weighted[label] = loaded_weighted[helicity]
                    covariance[label] = loaded_covariance[helicity]
        parity_ordinary = {
            helicity: linear_combine_series(
                {component: ordinary_by_component[component][helicity] for component in components},
                target_uu_coefficients,
            )
            for helicity in helicities
        }
        results[str(selection)] = _asymmetry_outputs(
            ordinary, None if direct_apar else weighted, None if direct_apar else covariance,
            uu_coefficients, ll_coefficients, parity_ordinary
        )
        _mask_selection_results(results[str(selection)], measurement["outputs"][selection])

    diagnostics: dict[str, BinSeries] = {}
    diagnostic_target = _target_component_coefficients(measurement, family, "unpolarized")
    diagnostic_coefficients = {
        _component_helicity_label(component, helicity):
        diagnostic_target[component] * helicity_uu_coefficients[helicity]
        for component in components for helicity in helicities
    }
    for name in measurement.get("diagnostics", []):
        component_series: dict[str, BinSeries] = {}
        for component in components:
            loaded = _load_component_series(
                groups, campaign_dir, analysis, str(name), family_id,
                component, helicities, orders,
                polarized_pdf_member, unpolarized_pdf_member, scale, order_coefficients,
            )
            for helicity in helicities:
                component_series[_component_helicity_label(component, helicity)] = loaded[helicity]
        diagnostics[str(name)] = linear_combine_series(component_series, diagnostic_coefficients)
    return results, diagnostics


def _dis_replica_sigma(
    variation_results: Mapping[tuple[int, int, float], Mapping[str, Mapping[str, Any]]],
    keys: Sequence[tuple[int, int, float]],
    selection: str,
    value_key: str,
    central: Sequence[float | None],
) -> list[float | None] | None:
    available = [
        variation_results[key][selection][value_key]
        for key in keys
        if key in variation_results
    ]
    if len(available) < 2:
        return None
    sigma: list[float | None] = []
    for index, centre in enumerate(central):
        values = [row[index] for row in available if row[index] is not None]
        if centre is None or len(values) < 2:
            sigma.append(None)
            continue
        mean = sum(float(value) for value in values) / len(values)
        sigma.append(
            math.sqrt(
                sum((float(value) - mean) ** 2 for value in values)
                / (len(values) - 1)
            )
        )
    return sigma


def aggregate_dis_uncertainties(
    variation_results: Mapping[tuple[int, int, float], Mapping[str, Mapping[str, Any]]]
) -> dict[str, dict[str, dict[str, Any]]]:
    central_key = (0, 0, 1.0)
    if central_key not in variation_results:
        raise CampaignError("The DIS uncertainty matrix has no central prediction")
    central_results = variation_results[central_key]
    polarized_members = sorted(
        key for key in variation_results if key[0] > 0 and key[1] == 0
        and math.isclose(key[2], 1.0)
    )
    unpolarized_members = sorted(
        key for key in variation_results if key[0] == 0 and key[1] > 0
        and math.isclose(key[2], 1.0)
    )
    scale_keys = [
        key for key in ((0, 0, 0.5), central_key, (0, 0, 2.0))
        if key in variation_results
    ]
    output: dict[str, dict[str, dict[str, Any]]] = {}
    for selection, central_result in central_results.items():
        output[selection] = {}
        for observable, value_key, error_key in (
            ("a1", "a1_values", "a1_errors"),
            ("apar", "apar_values", "apar_errors"),
        ):
            if value_key not in central_result:
                continue
            central = central_result[value_key]
            polarized = _dis_replica_sigma(
                variation_results, polarized_members, selection, value_key, central
            )
            unpolarized = _dis_replica_sigma(
                variation_results, unpolarized_members, selection, value_key, central
            )
            pdf: list[float | None] | None = None
            if polarized is not None or unpolarized is not None:
                pdf = []
                for index in range(len(central)):
                    components = [
                        array[index]
                        for array in (polarized, unpolarized)
                        if array is not None and array[index] is not None
                    ]
                    pdf.append(
                        math.sqrt(sum(float(value) ** 2 for value in components))
                        if components else None
                    )
            scale_down: list[float | None] | None = None
            scale_up: list[float | None] | None = None
            if len(scale_keys) >= 2:
                scale_down, scale_up = [], []
                for index, centre in enumerate(central):
                    values = [
                        variation_results[key][selection][value_key][index]
                        for key in scale_keys
                    ]
                    finite = [float(value) for value in values if value is not None]
                    if centre is None or not finite:
                        scale_down.append(None)
                        scale_up.append(None)
                    else:
                        scale_down.append(min(finite) - float(centre))
                        scale_up.append(max(finite) - float(centre))
            output[selection][observable] = {
                "pdf_68": pdf,
                "polarized_pdf_68": polarized,
                "unpolarized_pdf_68": unpolarized,
                "scale_down": scale_down,
                "scale_up": scale_up,
                "monte_carlo": central_result[error_key],
            }
    return output


def _attach_dis_bands(
    estimate: Any,
    values: Sequence[float | None],
    bands: Mapping[str, Any] | None,
) -> None:
    if not bands:
        return
    for index, value in enumerate(values, start=1):
        if value is None:
            continue
        bin_object = estimate.bin(index)
        pdf = bands.get("pdf_68")
        if pdf is not None and pdf[index - 1] is not None:
            error = float(pdf[index - 1])
            bin_object.setErr(-error, error, "pdf68")
        lower = bands.get("scale_down")
        upper = bands.get("scale_up")
        if (
            lower is not None and upper is not None
            and lower[index - 1] is not None and upper[index - 1] is not None
        ):
            bin_object.setErr(
                float(lower[index - 1]), float(upper[index - 1]), "hard_scale"
            )


def _asymmetry_family_objects(
    yoda: Any,
    measurement: Mapping[str, Any],
    snapshot: Mapping[str, Any],
    results: Mapping[str, Mapping[str, Any]],
    diagnostics: Mapping[str, BinSeries],
    annotations: Mapping[str, str],
    uncertainty_bands: Mapping[str, Mapping[str, Mapping[str, Any]]] | None = None,
) -> list[Any]:
    analysis = str(measurement["analysis"]["name"])
    objects: list[Any] = []
    for selection, result in results.items():
        paths = measurement["outputs"][selection]
        selected_annotations = _selection_annotations(measurement, selection, annotations)
        selection_bands = (
            uncertainty_bands.get(selection, {}) if uncertainty_bands else {}
        )
        if "a1_values" in result:
            a1_estimate = _estimate_from_values(
                yoda,
                result["edges"],
                paths["a1"],
                result["a1_values"],
                result["a1_errors"],
                {
                    **selected_annotations,
                    "Observable": str(paths.get("observable", snapshot.get("observable", "A1"))),
                    "Selection": selection,
                },
            )
            _attach_dis_bands(a1_estimate, result["a1_values"], selection_bands.get("a1"))
            objects.append(a1_estimate)
        apar_estimate = _estimate_from_values(
            yoda,
            result["edges"],
            paths["apar"],
            result["apar_values"],
            result["apar_errors"],
            {**_apar_annotations(selected_annotations), "Selection": selection},
        )
        _attach_dis_bands(
            apar_estimate, result["apar_values"], selection_bands.get("apar")
        )
        objects.extend(
            [
                apar_estimate,
                _estimate_from_values(
                    yoda,
                    result["edges"],
                    paths["parity_pp_mm"],
                    result["parity_pp_mm_values"],
                    result["parity_pp_mm_errors"],
                    {
                        **selected_annotations,
                        "Observable": "(PP-MM)/(PP+MM)",
                        "Selection": selection,
                    },
                ),
                _estimate_from_values(
                    yoda,
                    result["edges"],
                    paths["parity_pm_mp"],
                    result["parity_pm_mp_values"],
                    result["parity_pm_mp_errors"],
                    {
                        **selected_annotations,
                        "Observable": "(PM-MP)/(PM+MP)",
                        "Selection": selection,
                    },
                ),
                _series_estimate(
                    yoda,
                    result["sigma_uu"],
                    f"/{analysis}/SigmaUU_{selection}",
                    {**selected_annotations, "Observable": "sigma_UU", "Selection": selection},
                ),
                _series_estimate(
                    yoda,
                    result["sigma_ll"],
                    f"/{analysis}/SigmaLL_{selection}",
                    {**selected_annotations, "Observable": "sigma_LL", "Selection": selection},
                ),
            ]
        )
        supported = set(result.get("supported_bins", range(1, len(result["edges"]))))
        for estimate in objects[-(6 if "a1_values" in result else 5):]:
            _mask_unsupported_estimate(estimate, supported, len(result["edges"]) - 1)
    for name, series in diagnostics.items():
        objects.append(
            _series_estimate(
                yoda,
                series,
                f"/{analysis}/UU_{name}",
                {**annotations, "Observable": name, "Combination": "sigma_UU"},
                density=True,
            )
        )
    return objects


def _load_direct_unpolarized_products(
    groups: Mapping[
        tuple[str, int, int, float, str, str, str],
        Sequence[Mapping[str, Any]],
    ],
    campaign_dir: Path,
    measurement: Mapping[str, Any],
    family_id: str,
    family: Mapping[str, Any],
    polarized_pdf_member: int = 0,
    unpolarized_pdf_member: int = 0,
    scale: float = 1.0,
) -> tuple[dict[str, BinSeries], dict[str, BinSeries]]:
    analysis = str(measurement["analysis"]["name"])
    components = tuple(str(value) for value in family["components"])
    target_coefficients = _target_component_coefficients(
        measurement, family, "unpolarized"
    )
    orders = tuple(str(value) for value in family["orders"])
    order_coefficients = _order_combination_coefficients(measurement, orders)
    results: dict[str, BinSeries] = {}
    for selection, names in measurement["raw_observables"].items():
        selected_target = _target_component_coefficients(
            measurement, family, "unpolarized", str(selection)
        )
        loaded_components = {
            component: _load_component_series(
                groups, campaign_dir, analysis, str(names["ordinary"]),
                family_id, component, ("00",), orders,
                polarized_pdf_member, unpolarized_pdf_member, scale, order_coefficients,
            )["00"]
            for component in components
        }
        results[str(selection)] = linear_combine_series(
            loaded_components, selected_target
        )
    diagnostics: dict[str, BinSeries] = {}
    for name in measurement.get("diagnostics", []):
        loaded_components = {
            component: _load_component_series(
                groups, campaign_dir, analysis, str(name),
                family_id, component, ("00",), orders,
                polarized_pdf_member, unpolarized_pdf_member, scale, order_coefficients,
            )["00"]
            for component in components
        }
        diagnostics[str(name)] = linear_combine_series(
            loaded_components, target_coefficients
        )
    return results, diagnostics


def _direct_unpolarized_objects(
    yoda: Any,
    measurement: Mapping[str, Any],
    direct_results: Mapping[str, BinSeries],
    diagnostics: Mapping[str, BinSeries],
    nominal_results: Mapping[str, Mapping[str, Any]],
    annotations: Mapping[str, str],
) -> tuple[list[Any], dict[str, dict[str, list[float | None]]]]:
    analysis = str(measurement["analysis"]["name"])
    objects: list[Any] = []
    closure: dict[str, dict[str, list[float | None]]] = {}
    for selection, direct in direct_results.items():
        selected_annotations = _selection_annotations(measurement, selection, annotations)
        supported = _selection_supported_bins(measurement["outputs"][selection], len(direct.values))
        nominal = nominal_results[selection]["sigma_uu"]
        values: list[float | None] = []
        errors: list[float | None] = []
        for index in range(len(direct.values)):
            value, error = parity_residual(
                nominal.values[index],
                nominal.variances[index],
                direct.values[index],
                direct.variances[index],
            )
            values.append(value if index + 1 in supported else None)
            errors.append(error if index + 1 in supported else None)
        closure[selection] = {"values": values, "errors": errors}
        objects.extend(
            [
                _series_estimate(
                    yoda,
                    direct,
                    f"/{analysis}/SigmaUU_{selection}",
                    {
                        **selected_annotations,
                        "Observable": "direct sigma_00",
                        "Selection": selection,
                    },
                ),
                _estimate_from_values(
                    yoda,
                    direct.edges,
                    f"/{analysis}/UnpolarizedClosure_{selection}",
                    values,
                    errors,
                    {
                        **selected_annotations,
                        "Observable": "(sigma_UU-sigma_00)/(sigma_UU+sigma_00)",
                        "Selection": selection,
                    },
                ),
            ]
        )
        for estimate in objects[-2:]:
            _mask_unsupported_estimate(estimate, supported, len(direct.values))
    for name, series in diagnostics.items():
        objects.append(
            _series_estimate(
                yoda,
                series,
                f"/{analysis}/UU_{name}",
                {**annotations, "Observable": name, "Combination": "direct sigma_00"},
                density=True,
            )
        )
    return objects, closure


def _direct_unpolarized_rows(
    direct_results: Mapping[str, BinSeries],
    closure: Mapping[str, Mapping[str, Sequence[float | None]]],
    measurement: Mapping[str, Any] | None = None,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for selection, series in direct_results.items():
        output = measurement["outputs"][selection] if measurement is not None else {}
        axis = str(output.get("axis", "x"))
        supported_bins = _selection_supported_bins(output, len(series.values))
        for index, value in enumerate(series.values):
            supported = index + 1 in supported_bins
            rows.append(
                {
                    "selection": selection,
                    "bin": index + 1,
                    "axis": axis,
                    "bin_low": series.edges[index],
                    "bin_high": series.edges[index + 1],
                    f"{axis}_low": series.edges[index],
                    f"{axis}_high": series.edges[index + 1],
                    "supported": supported,
                    "sigma_00_pb": value if supported else None,
                    "sigma_00_stat_pb": math.sqrt(max(0.0, series.variances[index])) if supported else None,
                    "uu_00_closure": closure[selection]["values"][index],
                    "uu_00_closure_stat": closure[selection]["errors"][index],
                }
            )
    return rows


def _write_summary_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        return
    with path.open("w", encoding="utf-8", newline="") as stream:
        # Mixed x/Q2 and A1/direct-A_parallel rows have different meaningful columns.
        fieldnames = list(dict.fromkeys(key for row in rows for key in row))
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def postprocess_campaign(args: argparse.Namespace, measurement: Mapping[str, Any]) -> Path:
    campaign_dir = campaign_directory(str(measurement["id"]), args.tag)
    tracker_started = float(getattr(args, "_tracker_started_at", time.time()))
    manifest_path = campaign_dir / MANIFEST_NAME
    if not manifest_path.exists():
        raise CampaignError(f"No prepared campaign at {campaign_dir}")
    manifest = load_json(manifest_path)
    _assert_manifest_signatures_current(manifest, measurement)
    if not args.dry_run:
        write_campaign_monitor_files(
            campaign_dir,
            build_campaign_monitor_payload(
                measurement=measurement,
                tag=args.tag,
                phase="postprocessing",
                started_at=tracker_started,
                message="Combining normalized generator components and helicity bins.",
                herwig=existing_herwig_monitor(campaign_dir),
            ),
        )
    groups = _require_complete_matrix(manifest, measurement, campaign_dir)
    snapshot = load_reference_snapshot(measurement)
    yoda = _import_yoda()
    include_comparisons = bool(manifest["configuration"].get("comparisons", False))
    family_specs = campaign_family_specs(measurement, include_comparisons)
    family_products: dict[str, dict[str, Any]] = {}

    nominal_spec = family_specs["nominal"]
    configured_variations = [
        (int(point[0]), int(point[1]), float(point[2]))
        for point in manifest.get("configuration", {}).get(
            "variation_points", [[0, 0, 1.0]]
        )
    ]
    nominal_variations: dict[
        tuple[int, int, float], dict[str, dict[str, Any]]
    ] = {}
    nominal_variation_diagnostics: dict[
        tuple[int, int, float], dict[str, BinSeries]
    ] = {}
    for polarized_member, unpolarized_member, scale in configured_variations:
        results, diagnostics = _load_asymmetry_family_products(
            groups,
            campaign_dir,
            measurement,
            "nominal",
            nominal_spec,
            polarized_member,
            unpolarized_member,
            scale,
        )
        key = (polarized_member, unpolarized_member, scale)
        nominal_variations[key] = results
        nominal_variation_diagnostics[key] = diagnostics
    central_variation = (0, 0, 1.0)
    if central_variation not in nominal_variations:
        raise CampaignError("The DIS postprocessor has no central nominal variation")
    nominal_results = nominal_variations[central_variation]
    nominal_diagnostics = nominal_variation_diagnostics[central_variation]
    uncertainty_bands = aggregate_dis_uncertainties(nominal_variations)
    family_products["nominal"] = {
        "results": nominal_results,
        "diagnostics": nominal_diagnostics,
        "objects": _asymmetry_family_objects(
            yoda,
            measurement,
            snapshot,
            nominal_results,
            nominal_diagnostics,
            _family_annotations(measurement, "nominal", nominal_spec),
            uncertainty_bands,
        ),
        "rows": [
            row
            for selection, result in nominal_results.items()
            for row in _summary_rows(selection, result)
        ],
    }

    variation_objects: list[Any] = []
    analysis_name = str(measurement["analysis"]["name"])
    for key, results in nominal_variations.items():
        if key == central_variation:
            continue
        polarized_member, unpolarized_member, scale = key
        token = (
            f"p{polarized_member:03d}-u{unpolarized_member:03d}-"
            f"mu{_scale_token(scale)}"
        )
        variation_measurement = copy.deepcopy(dict(measurement))
        variation_measurement["analysis"]["name"] = (
            f"{analysis_name}/VARIATIONS/{token}"
        )
        for selection, paths in variation_measurement["outputs"].items():
            for observable in ("a1", "apar", "parity_pp_mm", "parity_pm_mp"):
                if observable not in paths:
                    continue
                paths[observable] = (
                    f"/{analysis_name}/VARIATIONS/{token}/"
                    f"{observable}_{selection}"
                )
        annotations = {
            **_family_annotations(measurement, "nominal", nominal_spec),
            "PolarizedPDFMember": str(polarized_member),
            "UnpolarizedPDFMember": str(unpolarized_member),
            "HardScaleFactor": str(scale),
        }
        variation_objects.extend(
            _asymmetry_family_objects(
                yoda,
                variation_measurement,
                snapshot,
                results,
                {},
                annotations,
            )
        )

    for family_id, family in family_specs.items():
        if family_id == "nominal":
            continue
        annotations = _family_annotations(measurement, family_id, family)
        if family["postprocess"] == "helicity_asymmetry":
            results, diagnostics = _load_asymmetry_family_products(
                groups, campaign_dir, measurement, family_id, family
            )
            family_products[family_id] = {
                "results": results,
                "diagnostics": diagnostics,
                "objects": _asymmetry_family_objects(
                    yoda, measurement, snapshot, results, diagnostics, annotations
                ),
                "rows": [
                    row
                    for selection, result in results.items()
                    for row in _summary_rows(selection, result)
                ],
            }
            continue
        direct_results, diagnostics = _load_direct_unpolarized_products(
            groups, campaign_dir, measurement, family_id, family
        )
        objects, closure = _direct_unpolarized_objects(
            yoda,
            measurement,
            direct_results,
            diagnostics,
            nominal_results,
            annotations,
        )
        family_products[family_id] = {
            "direct_results": direct_results,
            "diagnostics": diagnostics,
            "closure": closure,
            "objects": objects,
            "rows": _direct_unpolarized_rows(direct_results, closure, measurement),
        }

    postprocess_dir = campaign_dir / "postprocess"
    postprocess_dir.mkdir(parents=True, exist_ok=True)
    yoda_outputs = {
        family_id: postprocess_dir
        / (
            f"{measurement['id']}-analyzed.yoda"
            if family_id == "nominal"
            else f"{measurement['id']}-{family_id}.yoda"
        )
        for family_id in family_specs
    }
    variation_output = postprocess_dir / f"{measurement['id']}-variations.yoda"
    yoda_output = yoda_outputs["nominal"]
    if args.dry_run:
        for family_id, output in yoda_outputs.items():
            print(f"Would write {len(family_products[family_id]['objects'])} objects to {output}")
        if variation_objects:
            print(f"Would write {len(variation_objects)} objects to {variation_output}")
        return yoda_output
    for family_id, output in yoda_outputs.items():
        _write_yoda_objects(yoda, family_products[family_id]["objects"], output)
        if not _nonempty(output):
            raise CampaignError(f"Postprocessing produced empty output {output}")
    if variation_objects:
        _write_yoda_objects(yoda, variation_objects, variation_output)
        if not _nonempty(variation_output):
            raise CampaignError(
                f"Postprocessing produced empty variation output {variation_output}"
            )

    rows = family_products["nominal"]["rows"]
    summary = {
        "measurement": measurement["id"],
        "tag": args.tag,
        "created_at": utc_now(),
        "assumptions": list(measurement.get("physics", {}).get("assumptions", [])),
        "sign_convention": {
            "sigma_UU": "(PP+PM+MP+MM)/4",
            "sigma_LL": "(PP+MM-PM-MP)/4",
        },
        "target_combination": measurement["combination"].get("target_components", {}),
        "order_combination": {
            "coefficients": _order_combination_coefficients(measurement, nominal_spec["orders"]),
            "formula": _order_combination_formula(
                _order_combination_coefficients(measurement, nominal_spec["orders"])
            ),
            "input_convention": measurement.get("physics", {}).get(
                "order_input_convention", "coefficients multiply stored normalized order bins"
            ),
        },
        "output_target_combinations": {
            selection: {
                channel: _target_component_coefficients(
                    measurement, nominal_spec, channel, selection
                )
                for channel in ("unpolarized", "longitudinal")
            }
            for selection in measurement["outputs"]
        },
        "reference": snapshot["provenance"],
        "uncertainties": uncertainty_bands,
        "variation_points": [list(point) for point in configured_variations],
        "masked_bins": {
            selection: [
                index + 1
                for index, value in enumerate(result.get("a1_values", result["apar_values"]))
                if value is None
            ]
            for selection, result in nominal_results.items()
        },
        "bins": rows,
        "diagnostics": {
            name: {
                "edges": series.edges,
                "values": series.values,
                "errors": [math.sqrt(max(0.0, value)) for value in series.variances],
            }
            for name, series in nominal_diagnostics.items()
        },
        "prediction_families": [
            {
                "id": family_id,
                "label": family["label"],
                "postprocess": family["postprocess"],
                "yoda": str(yoda_outputs[family_id].relative_to(campaign_dir)),
            }
            for family_id, family in family_specs.items()
        ],
    }
    atomic_write_json(postprocess_dir / "summary.json", summary)
    _write_summary_csv(postprocess_dir / "summary.csv", rows)

    comparison_summary: dict[str, Any] = {
        "measurement": measurement["id"],
        "tag": args.tag,
        "created_at": utc_now(),
        "families": {},
    }
    for family_id, family in family_specs.items():
        if family_id == "nominal":
            continue
        product = family_products[family_id]
        family_summary: dict[str, Any] = {
            "label": family["label"],
            "postprocess": family["postprocess"],
            "perturbative_order": family.get("perturbative_order"),
            "real_emission_spin_density": family.get("real_emission_spin_density"),
            "shower_spin_correlations": family.get("shower_spin_correlations"),
            "beam_polarization": family.get("beam_polarization"),
            "order_coefficients": _order_combination_coefficients(measurement, family["orders"]),
            "yoda": str(yoda_outputs[family_id].relative_to(campaign_dir)),
            "bins": product["rows"],
            "diagnostics": {
                name: {
                    "edges": series.edges,
                    "values": series.values,
                    "errors": [math.sqrt(max(0.0, value)) for value in series.variances],
                }
                for name, series in product["diagnostics"].items()
            },
        }
        if "results" in product:
            family_summary["masked_bins"] = {
                selection: [
                    index + 1
                    for index, value in enumerate(result.get("a1_values", result["apar_values"]))
                    if value is None
                ]
                for selection, result in product["results"].items()
            }
        if "closure" in product:
            family_summary["unpolarized_closure"] = product["closure"]
        comparison_summary["families"][family_id] = family_summary
        _write_summary_csv(postprocess_dir / f"summary-{family_id}.csv", product["rows"])
    comparison_summary_path: Path | None = None
    if include_comparisons:
        comparison_summary_path = postprocess_dir / "comparison-summary.json"
        atomic_write_json(comparison_summary_path, comparison_summary)

    predictions = [
        {
            "family": family_id,
            "label": str(family["label"]),
            "postprocess": str(family["postprocess"]),
            "yoda": str(yoda_outputs[family_id].relative_to(campaign_dir)),
        }
        for family_id, family in family_specs.items()
    ]
    manifest["postprocess"] = {
        "created_at": utc_now(),
        "yoda": str(yoda_output.relative_to(campaign_dir)),
        "predictions": predictions,
        "summary_json": "postprocess/summary.json",
        "summary_csv": "postprocess/summary.csv",
    }
    if variation_objects:
        manifest["postprocess"]["variations_yoda"] = str(
            variation_output.relative_to(campaign_dir)
        )
    if comparison_summary_path is not None:
        manifest["postprocess"]["comparison_summary_json"] = str(
            comparison_summary_path.relative_to(campaign_dir)
        )
    manifest["updated_at"] = utc_now()
    manifest["history"].append({"at": utc_now(), "action": "postprocess"})
    atomic_write_json(manifest_path, manifest)
    write_campaign_monitor_files(
        campaign_dir,
        build_campaign_monitor_payload(
            measurement=measurement,
            tag=args.tag,
            phase="postprocessed",
            started_at=tracker_started,
            message=f"Wrote {len(yoda_outputs)} normalized prediction family output(s).",
            herwig=existing_herwig_monitor(campaign_dir),
        ),
    )
    print(
        f"Wrote normalized-bin {measurement['id']} observables for "
        f"{len(yoda_outputs)} prediction family output(s) below {postprocess_dir}"
    )
    return yoda_output


def write_plot_indexes(
    output_dir: Path,
    measurement: Mapping[str, Any],
    plot_scripts: Sequence[Path],
) -> Path:
    """Wrap sequentially generated Rivet plots in a small static HTML report.

    ``rivet-mkhtml --dry-run`` is used deliberately because its normal parallel
    driver opens a local multiprocessing socket, which is unavailable in some
    batch and sandbox environments.  The generated plotting scripts are still
    Rivet's; this function supplies only the indexes that the dry-run mode omits.
    """

    def plot_cards(directory: Path, scripts: Sequence[Path]) -> str:
        cards: list[str] = []
        for script in sorted(scripts):
            stem = script.stem
            png = directory / f"{stem}.png"
            pdf = directory / f"{stem}.pdf"
            if not _nonempty(png):
                raise CampaignError(f"Generated Rivet plot is missing or empty: {png}")
            escaped_stem = html.escape(stem)
            links = [f'<a href="{html.escape(png.name)}">PNG</a>']
            if _nonempty(pdf):
                links.append(f'<a href="{html.escape(pdf.name)}">PDF</a>')
            cards.append(
                '<section class="plot">'
                f"<h2>{escaped_stem}</h2>"
                f'<a href="{html.escape(png.name)}">'
                f'<img src="{html.escape(png.name)}" alt="{escaped_stem}"></a>'
                f"<p>{' &middot; '.join(links)}</p>"
                "</section>"
            )
        return "\n".join(cards)

    def document(title: str, body: str) -> str:
        escaped_title = html.escape(title)
        return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{escaped_title}</title>
  <style>
    body {{ font-family: sans-serif; margin: 2rem auto; max-width: 1100px; padding: 0 1rem; }}
    .plot {{ border-top: 1px solid #ccc; margin-top: 2rem; padding-top: 1rem; }}
    img {{ height: auto; max-width: 100%; }}
  </style>
</head>
<body>
  <h1>{escaped_title}</h1>
  <p>Plots generated by Rivet from the postprocessed physical-helicity campaign.</p>
  {body}
</body>
</html>
"""

    grouped: dict[Path, list[Path]] = {}
    for script in plot_scripts:
        try:
            relative_directory = script.parent.relative_to(output_dir)
        except ValueError as exc:
            raise CampaignError(f"Plot script lies outside output directory: {script}") from exc
        grouped.setdefault(relative_directory, []).append(script)

    title = str(measurement.get("title", measurement["id"]))
    root_sections: list[str] = []
    for relative_directory, scripts in sorted(grouped.items(), key=lambda item: str(item[0])):
        directory = output_dir / relative_directory
        directory.mkdir(parents=True, exist_ok=True)
        cards = plot_cards(directory, scripts)
        if relative_directory == Path("."):
            root_sections.append(cards)
            continue
        child_title = f"{title} — {relative_directory.as_posix()}"
        child_index = directory / "index.html"
        child_index.write_text(document(child_title, cards), encoding="utf-8")
        link = html.escape((relative_directory / "index.html").as_posix())
        label = html.escape(relative_directory.as_posix())
        root_sections.append(f'<p><a href="{link}">{label}</a> ({len(scripts)} plots)</p>')

    index = output_dir / "index.html"
    index.write_text(document(title, "\n".join(root_sections)), encoding="utf-8")
    return index


def plot_script_has_finite_y(script: Path) -> bool:
    """Whether a generated Rivet plot has finite ordinates and axis limits.

    Fully masked theory diagnostics are expected in low-statistics smoke runs.
    Rivet correctly writes ``nan`` for them, but its automatic limits can also
    become non-finite when only part of a low-statistics series is masked.
    Matplotlib cannot render either case. Such panels are omitted from the
    HTML rather than being converted to misleading zeros or arbitrary ranges.
    """
    source = script.read_text(encoding="utf-8")
    nonfinite_token = re.compile(
        r"(?<![A-Za-z0-9_])(?:nan|[+-]?inf)(?![A-Za-z0-9_])",
        re.IGNORECASE,
    )
    for line in source.splitlines():
        if line.startswith(("xLims =", "yLims =")) and nonfinite_token.search(line):
            return False
    data_script = script.with_name(f"{script.stem}__data.py")
    if not data_script.is_file():
        return True
    namespace = runpy.run_path(str(data_script))
    yvalues = namespace.get("yvals")
    if not isinstance(yvalues, Mapping):
        return True
    return any(
        math.isfinite(float(value))
        for values in yvalues.values()
        for value in values
    )


def apply_published_reference_coordinates(
    script: Path, points: Sequence[Mapping[str, Any]] | None,
) -> None:
    """Place the reference markers at published means, keeping physical bin edges.

    YODA's plot generator forces all estimates to their bin centres. Reference
    means are display coordinates; the theory still represents bin integrals.
    """
    if not points:
        return
    if any(not all(key in point and math.isfinite(float(point[key])) for key in ("plot_x", "value"))
           for point in points):
        raise CampaignError("Published reference coordinates require finite means and values")
    source = script.read_text(encoding="utf-8")
    marker = "# curve from input yoda files in main panel"
    if marker not in source:
        raise CampaignError(f"Cannot locate the Rivet reference plotting point in {script}")
    coordinates = [float(point["plot_x"]) for point in points]
    values = [float(point["value"]) for point in points]
    adjustment = [
        "# Use published reference means, preserving theory bin edges and integrals",
        "_published_ref_key = next(iter(dataf['yvals']))",
        f"_published_ref_x = np.asarray({coordinates!r}, dtype=float)",
        f"_published_ref_y = np.asarray({values!r}, dtype=float)",
        "_published_ref_actual_y = np.asarray(dataf['yvals'][_published_ref_key], dtype=float)",
        "if (_published_ref_actual_y.shape != _published_ref_y.shape or",
        "        not np.allclose(_published_ref_actual_y, _published_ref_y, rtol=1e-8, atol=1e-12)):",
        "    raise RuntimeError('Published reference points do not match the first plotted data curve')",
        "_published_ref_edges = np.asarray(dataf['xedges'][_published_ref_key], dtype=float)",
        "if (_published_ref_edges.size != _published_ref_x.size + 1 or",
        "        np.any(_published_ref_x < _published_ref_edges[:-1]) or",
        "        np.any(_published_ref_x > _published_ref_edges[1:])):",
        "    raise RuntimeError('Published reference means lie outside their physical bins')",
        "dataf['xpoints'][_published_ref_key] = _published_ref_x.tolist()",
        "dataf['xerrs'][_published_ref_key] = [",
        "    (_published_ref_x - _published_ref_edges[:-1]).tolist(),",
        "    (_published_ref_edges[1:] - _published_ref_x).tolist()]",
        "if 'ref_xerrs' in dataf:",
        "    dataf['ref_xerrs'] = dataf['xerrs'][_published_ref_key]",
    ]
    source = source.replace(marker, "\n".join(adjustment) + "\n\n" + marker, 1)
    script.write_text(source, encoding="utf-8")


def add_experimental_error_overlay(
    script: Path,
    points: Sequence[Mapping[str, Any]] | None,
    show_statistical: bool = False,
) -> None:
    """Optionally add statistical bars and retain low-Q2 display markers.

    The Rivet reference curve already carries the total experimental
    uncertainty.  Statistical-only bars are therefore an opt-in display
    component rather than a second default error bar on the data.
    """
    if not points:
        return
    required = ("plot_x", "value", "stat") if show_statistical else ("plot_x", "value")
    normalized = [
        point for point in points
        if all(key in point for key in required)
        and all(math.isfinite(float(point[key])) for key in required)
    ]
    if not normalized:
        return
    source = script.read_text(encoding="utf-8")
    marker = "# set plot metadata as defined above"
    if marker not in source:
        raise CampaignError(f"Cannot locate the Rivet plot insertion point in {script}")
    x_values = [float(point["plot_x"]) for point in normalized]
    y_values = [float(point["value"]) for point in normalized]
    low_q2 = [
        point for point in normalized
        if bool(point.get("nonperturbative_extrapolation", False))
    ]
    if not show_statistical and not low_q2:
        return
    overlay = ["# Experimental display components added by the campaign runner"]
    if show_statistical:
        stat_values = [float(point["stat"]) for point in normalized]
        lower_legend = re.search(r"\bloc\s*=\s*['\"]lower", source) is not None
        stat_note_y = 0.97 if lower_legend else (0.84 if low_q2 else 0.025)
        stat_note_x = 0.98 if lower_legend else 0.02
        stat_note_horizontal = "right" if lower_legend else "left"
        stat_note_vertical = "top" if lower_legend else "bottom"
        overlay.extend(
            [
                f"_data_x = {x_values!r}",
                f"_data_y = {y_values!r}",
                f"_data_stat = {stat_values!r}",
                "ax.errorbar(_data_x, _data_y, yerr=_data_stat, fmt='none',",
                "            ecolor='black', elinewidth=1.7, capsize=2.5,",
                "            capthick=1.2, zorder=10)",
                f"ax.text({stat_note_x}, {stat_note_y}, 'inner bars: statistical; outer bars: total',",
                f"        transform=ax.transAxes, fontsize=7, ha={stat_note_horizontal!r}, va={stat_note_vertical!r})",
            ]
        )
    if low_q2:
        overlay.extend(
            [
                f"_lowq_x = {[float(point['plot_x']) for point in low_q2]!r}",
                f"_lowq_y = {[float(point['value']) for point in low_q2]!r}",
                "ax.scatter(_lowq_x, _lowq_y, marker='s', s=38,",
                "           facecolors='none', edgecolors='#AA3377',",
                "           linewidths=1.3, zorder=11)",
                "ax.text(0.02, 0.89, r'open squares: $\\langle Q^2\\rangle<1$ GeV$^2$; display only',",
                "        transform=ax.transAxes, fontsize=7, ha='left', va='bottom',",
                "        color='#AA3377')",
            ]
        )
    source = source.replace(marker, "\n".join(overlay) + "\n\n" + marker, 1)
    script.write_text(source, encoding="utf-8")


def ensure_plot_canvas_draw(script: Path) -> None:
    """Use complete MathText fonts and initialize layout before each format save."""
    source = script.read_text(encoding="utf-8")
    # Rivet's no-TeX fallback selects DejaVu Sans but leaves the style's custom
    # Palatino math fonts active. Missing font variants can hide math axis labels.
    source = re.sub(
        r"(?m)^([ \t]*)(plt\.rcParams\[['\"]text\.usetex['\"]\]\s*=\s*False)\s*$",
        r"\1\2\n\1plt.rcParams['mathtext.fontset'] = 'dejavusans'",
        source,
    )
    source, count = re.subn(
        r"(?m)^([ \t]*)plt\.savefig\(",
        r"\1fig.canvas.draw()\n\1plt.savefig(",
        source,
    )
    if count == 0:
        raise CampaignError(f"Cannot locate the Rivet figure save calls in {script}")
    script.write_text(source, encoding="utf-8")


def add_theory_uncertainty_overlay(
    script: Path,
    bands: Mapping[str, Any] | None,
    theory_label: str,
) -> None:
    """Render MC, PDF, and hard-scale uncertainties as distinct components.

    YODA preserves named uncertainty components, but Rivet's standard plotting
    path reduces them to a single quadrature error.  This deterministic layer
    keeps those components separate in the generated paper plots.  The caller
    supplies the nominal prediction's generated-series label explicitly since
    reference plots conventionally place ``Data`` first in ``dataf['yvals']``.
    """
    if not bands:
        return

    def component(name: str) -> list[float | None] | None:
        raw = bands.get(name)
        if not isinstance(raw, Sequence) or isinstance(raw, (str, bytes)):
            return None
        values: list[float | None] = []
        for value in raw:
            if value is None:
                values.append(None)
                continue
            numeric = float(value)
            values.append(numeric if math.isfinite(numeric) else None)
        return values if any(value is not None for value in values) else None

    monte_carlo = component("monte_carlo")
    pdf = component("pdf_68")
    scale_down = component("scale_down")
    scale_up = component("scale_up")
    if monte_carlo is None and pdf is None and (scale_down is None or scale_up is None):
        return

    def array_literal(values: Sequence[float | None] | None) -> str:
        if values is None:
            return "None"
        return "[" + ", ".join(
            "nan" if value is None else repr(value) for value in values
        ) + "]"

    source = script.read_text(encoding="utf-8")
    curve_marker = "# curve from input yoda files in main panel"
    legend_marker = "legend_items = list(legend_handles.values())"
    if curve_marker not in source or legend_marker not in source:
        raise CampaignError(f"Cannot locate Rivet uncertainty insertion points in {script}")

    pre_curve = [
        "# Named theory-uncertainty components added by the campaign runner",
        f"_theory_label = {str(theory_label)!r}",
        "if _theory_label not in dataf['yvals']:",
        "    raise KeyError(f'Nominal theory series {_theory_label!r} is absent from generated plot data')",
        "_theory_y = np.asarray(dataf['yvals'][_theory_label], dtype=float)",
        "_theory_x = np.asarray(dataf['xedges'][_theory_label], dtype=float)",
        f"_theory_mc = {array_literal(monte_carlo)}",
        f"_theory_pdf = {array_literal(pdf)}",
        f"_theory_scale_dn = {array_literal(scale_down)}",
        f"_theory_scale_up = {array_literal(scale_up)}",
        "_theory_mc_handle = None",
        "_theory_pdf_handle = None",
        "_theory_scale_handle = None",
        "styles[_theory_label]['yerrorbars'] = 0",
        "def _theory_step(values):",
        "    return np.insert(np.asarray(values, dtype=float), 0, float(values[0]))",
        "if (_theory_scale_dn is not None and _theory_scale_up is not None",
        "        and len(_theory_scale_dn) == len(_theory_y)",
        "        and len(_theory_scale_up) == len(_theory_y)):",
        "    _scale_lo = _theory_y + np.asarray(_theory_scale_dn, dtype=float)",
        "    _scale_hi = _theory_y + np.asarray(_theory_scale_up, dtype=float)",
        "    _theory_scale_handle = ax.fill_between(",
        "        _theory_x, _theory_step(_scale_lo), _theory_step(_scale_hi),",
        "        step='pre', color='#EE7733', alpha=0.22, linewidth=0, zorder=1)",
        "if _theory_pdf is not None and len(_theory_pdf) == len(_theory_y):",
        "    _pdf_err = np.asarray(_theory_pdf, dtype=float)",
        "    _theory_pdf_handle = ax.fill_between(",
        "        _theory_x, _theory_step(_theory_y-_pdf_err),",
        "        _theory_step(_theory_y+_pdf_err), step='pre',",
        "        color='#0077BB', alpha=0.25, linewidth=0, zorder=2)",
        "if _theory_mc is not None and len(_theory_mc) == len(_theory_y):",
        "    _mc_err = np.asarray(_theory_mc, dtype=float)",
        "    _theory_mc_handle = ax.fill_between(",
        "        _theory_x, _theory_step(_theory_y-_mc_err),",
        "        _theory_step(_theory_y+_mc_err), step='pre',",
        "        color=styles[_theory_label]['color'], alpha=0.24,",
        "        linewidth=0, zorder=3)",
    ]
    source = source.replace(
        curve_marker, "\n".join(pre_curve) + "\n\n" + curve_marker, 1
    )

    pre_legend = [
        "# Keep the uncertainty components distinct in the legend",
        "_theory_titles = dict(zip(dataf.get('add_legend_handle', []), labels['legend'][0]))",
        "labels['legend'][0][:] = [_theory_titles.get(key, key) for key in legend_handles]",
        "if _theory_mc_handle is not None:",
        "    legend_handles['__mc_stat'] = _theory_mc_handle",
        "    labels['legend'][0].append('MC statistical')",
        "if _theory_pdf_handle is not None:",
        "    legend_handles['__pdf68'] = _theory_pdf_handle",
        "    labels['legend'][0].append('PDF 68%')",
        "if _theory_scale_handle is not None:",
        "    legend_handles['__hard_scale'] = _theory_scale_handle",
        "    labels['legend'][0].append(r'hard scale $\\mu/\\mu_0=0.5,2$')",
    ]
    source = source.replace(
        legend_marker, "\n".join(pre_legend) + "\n\n" + legend_marker, 1
    )
    script.write_text(source, encoding="utf-8")


def _fixed_reference_overlay_points(
    measurement: Mapping[str, Any],
    snapshot: Mapping[str, Any],
    plot_stem: str,
) -> list[dict[str, Any]] | None:
    datasets = snapshot.get("datasets")
    if isinstance(datasets, list):
        for dataset in datasets:
            if Path(str(dataset.get("rivet_path", ""))).name != plot_stem:
                continue
            axis = str(dataset.get("plot_axis", "q2_mean"))
            if axis not in {"x_mean", "q2_mean"}:
                raise CampaignError(f"Unsupported multi-dataset reference plot axis {axis!r}")
            return [
                {**dict(point), "plot_x": float(point[axis])}
                for point in dataset.get("points", [])
            ]
        return None
    reference_path = str(snapshot.get("rivet_path", ""))
    if Path(reference_path).name != plot_stem:
        return None
    return [
        {**dict(point), "plot_x": float(point["x_mean"])}
        for point in snapshot.get("points", [])
    ]


def _fixed_theory_uncertainty_bands(
    measurement: Mapping[str, Any],
    summary: Mapping[str, Any],
    plot_stem: str,
) -> Mapping[str, Any] | None:
    uncertainties = summary.get("uncertainties", {})
    for selection, outputs in measurement.get("outputs", {}).items():
        if not isinstance(outputs, Mapping):
            continue
        for observable in ("a1", "apar"):
            if Path(str(outputs.get(observable, ""))).name != plot_stem:
                continue
            selected = uncertainties.get(selection, {})
            if isinstance(selected, Mapping):
                bands = selected.get(observable)
                return bands if isinstance(bands, Mapping) else None
    return None


def plot_campaign(args: argparse.Namespace, measurement: Mapping[str, Any]) -> Path:
    campaign_dir = campaign_directory(str(measurement["id"]), args.tag)
    tracker_started = float(getattr(args, "_tracker_started_at", time.time()))
    manifest_path = campaign_dir / MANIFEST_NAME
    if not manifest_path.exists():
        raise CampaignError(f"No prepared campaign at {campaign_dir}")
    manifest = load_json(manifest_path)
    plot_metadata_refresh = _assert_manifest_signatures_current(
        manifest,
        measurement,
        allow_plot_metadata_refresh=bool(
            getattr(args, "allow_plot_metadata_refresh", False)
        ),
    )
    postprocess = manifest.get("postprocess", {})
    all_prediction_entries = postprocess.get("predictions")
    if not isinstance(all_prediction_entries, list) or not all_prediction_entries:
        all_prediction_entries = [
            {
                "family": "nominal",
                "label": "NLO+PS (polarized; full spin)",
                "yoda": postprocess.get(
                    "yoda", f"postprocess/{measurement['id']}-analyzed.yoda"
                ),
            }
        ]
    plot_comparisons = bool(getattr(args, "plot_comparisons", False))
    if plot_comparisons:
        prediction_entries = all_prediction_entries
    else:
        prediction_entries = [
            entry for entry in all_prediction_entries
            if str(entry.get("family", "")) == "nominal"
        ]
        if not prediction_entries:
            prediction_entries = all_prediction_entries[:1]
    predictions: list[tuple[Path, str]] = []
    nominal_theory_label: str | None = None
    for entry in prediction_entries:
        prediction = campaign_dir / str(entry["yoda"])
        if not _nonempty(prediction):
            raise CampaignError(f"Missing postprocessed YODA file {prediction}")
        predictions.append(
            (prediction, str(entry.get("label", entry.get("family", "HerwigPol"))))
        )
        if str(entry.get("family", "")) == "nominal":
            nominal_theory_label = prediction.name
    if nominal_theory_label is None:
        nominal_theory_label = predictions[0][0].name
    runtime = manifest.get("runtime") or preflight_runtime(measurement)
    output_dir = campaign_dir / "plots"
    safe_wrapper = DISPOL_ROOT / "scripts" / "rivet_mkhtml_safe.py"
    command = [
        sys.executable,
        str(safe_wrapper),
        runtime["tools"]["rivet-mkhtml"],
        "--dry-run",
        "-o",
        str(output_dir),
    ]
    plotting = measurement.get("plotting", {})
    if plot_comparisons and len(prediction_entries) > 1:
        plotting = _comparison_profile(measurement).get("plotting", plotting)
    for pattern in plotting.get("unmatch", []):
        command.extend(["-M", str(pattern)])
    for prediction, label in predictions:
        command.append(f"{prediction}:Title={label}")
    if args.dry_run:
        print(" ".join(command))
        return output_dir
    write_campaign_monitor_files(
        campaign_dir,
        build_campaign_monitor_payload(
            measurement=measurement,
            tag=args.tag,
            phase="plotting",
            started_at=tracker_started,
            message="Generating Rivet HTML and plot images.",
            herwig=existing_herwig_monitor(campaign_dir),
        ),
    )
    ensure_reference_yoda(measurement, campaign_dir)
    environment = analysis_environment(measurement, campaign_dir, runtime)
    mpl_cache = campaign_dir / "work" / "matplotlib-cache"
    mpl_cache.mkdir(parents=True, exist_ok=True)
    environment["MPLCONFIGDIR"] = str(mpl_cache)
    environment["DISPOL_FORCE_NO_ERROR_BANDS"] = "1"
    _run_logged(command, campaign_dir, environment, campaign_dir / "logs" / "rivet-mkhtml-generate.log")
    plot_scripts = sorted(
        path for path in output_dir.rglob("*.py") if not path.name.endswith("__data.py")
    )
    if not plot_scripts:
        raise CampaignError(f"rivet-mkhtml did not generate plot scripts below {output_dir}")
    script_log = campaign_dir / "logs" / "rivet-plot-scripts.log"
    rendered_scripts: list[Path] = []
    snapshot = load_reference_snapshot(measurement)
    summary_path = campaign_dir / "postprocess" / "summary.json"
    summary = load_json(summary_path) if summary_path.is_file() else {}
    with script_log.open("w", encoding="utf-8") as log:
        for script in plot_scripts:
            log.write(f"script: {script}\n")
            log.flush()
            reference_points = _fixed_reference_overlay_points(measurement, snapshot, script.stem)
            apply_published_reference_coordinates(script, reference_points)
            add_theory_uncertainty_overlay(
                script,
                _fixed_theory_uncertainty_bands(measurement, summary, script.stem),
                nominal_theory_label,
            )
            add_experimental_error_overlay(
                script,
                reference_points,
                show_statistical=bool(
                    getattr(args, "plot_data_components", False)
                ),
            )
            ensure_plot_canvas_draw(script)
            if not plot_script_has_finite_y(script):
                log.write(
                    "skipped: non-renderable finite-data/axis-limit state\n"
                )
                continue
            completed = subprocess.run(
                [sys.executable, str(script)],
                cwd=script.parent,
                env=environment,
                stdout=log,
                stderr=subprocess.STDOUT,
            )
            if completed.returncode != 0:
                raise CampaignError(
                    f"Generated Rivet plot script failed with status {completed.returncode}: {script}; see {script_log}"
                )
            rendered_scripts.append(script)
    if not any(output_dir.rglob("*.png")):
        raise CampaignError(f"Generated Rivet scripts produced no PNG plots below {output_dir}")
    index = write_plot_indexes(output_dir, measurement, rendered_scripts)
    if not _nonempty(index):
        raise CampaignError(f"Could not create Rivet plot index {index}")
    manifest["plots"] = {
        "created_at": utc_now(),
        "index": str(index.relative_to(campaign_dir)),
        "plot_metadata_sha256": sha256_file(
            resolve_dispol_path(str(measurement["analysis"]["plot"]))
        ),
    }
    if plot_metadata_refresh is not None:
        manifest["plots"]["presentation_only_refresh"] = plot_metadata_refresh
    manifest["updated_at"] = utc_now()
    manifest["history"].append(
        {
            "at": utc_now(),
            "action": (
                "plot-metadata-refresh"
                if plot_metadata_refresh is not None else "plot"
            ),
        }
    )
    atomic_write_json(manifest_path, manifest)
    write_campaign_monitor_files(
        campaign_dir,
        build_campaign_monitor_payload(
            measurement=measurement,
            tag=args.tag,
            phase="complete",
            started_at=tracker_started,
            message=f"Rivet HTML is available at {index.relative_to(campaign_dir)}.",
            herwig=existing_herwig_monitor(campaign_dir),
        ),
    )
    print(f"Wrote Rivet HTML to {index}")
    return output_dir


# ---------------------------------------------------------------------------
# Command-line interface
# ---------------------------------------------------------------------------


def _add_measurement(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--measurement", required=True, help="Measurement registry id")


def _add_tag(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--tag", required=True, help="Immutable campaign tag")


def _add_campaign_options(parser: argparse.ArgumentParser) -> None:
    _add_measurement(parser)
    _add_tag(parser)
    parser.add_argument("--jobs", type=int, help="Concurrent Herwig shard jobs")
    parser.add_argument("--shards", type=int, help="Shards per helicity/order logical job")
    parser.add_argument("--seed-base", type=int, help="First deterministic random seed")
    parser.add_argument("--posnlo-events", type=int, help="POSNLO events per physical helicity")
    parser.add_argument("--negnlo-events", type=int, help="NEGNLO events per physical helicity")
    parser.add_argument(
        "--lo-events",
        type=int,
        help="LO events per physical helicity when --comparisons is enabled; defaults to POSNLO events",
    )
    parser.add_argument(
        "--comparisons",
        action="store_true",
        help="Add registry-defined unpolarized, LO, and spin-treatment comparison families",
    )
    parser.add_argument(
        "--profile",
        choices=("central", "paper"),
        default="central",
        help=(
            "central runs the nominal prediction; paper also enables the "
            "comparison families, PDF replicas, and hard-scale envelope"
        ),
    )
    parser.add_argument(
        "--polarized-pdf-members",
        help="central, all, or comma-separated NNPDFpol2.0 member numbers",
    )
    parser.add_argument(
        "--unpolarized-pdf-members",
        help="central, all, or comma-separated NNPDF4.0 member numbers",
    )
    parser.add_argument(
        "--scales",
        help="central, all, or comma-separated hard-scale factors from 0.5,1,2",
    )
    parser.add_argument(
        "--progress-interval",
        type=float,
        default=5.0,
        help="Seconds between progress refreshes; a negative value disables live refreshes",
    )
    parser.add_argument(
        "--max-listed",
        type=int,
        default=12,
        help="Maximum number of logical runs and active shards in the live display",
    )
    parser.add_argument("--smoke", action="store_true", help="Use the registry smoke event count")
    parser.add_argument("--dry-run", action="store_true", help="Print the resolved plan without executing")
    parser.add_argument(
        "--recover-failed",
        "--rerun-failed-random-seed",
        dest="recover_failed",
        action="store_true",
        help="Retry failed shards with fresh, non-reused seeds",
    )


def _add_plot_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--plot-comparisons",
        action="store_true",
        help=(
            "Overlay postprocessed comparison-family predictions; by default "
            "only the nominal NLO+PS (polarized; full spin) prediction is plotted"
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


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    subparsers.add_parser("list", help="List discovered measurement descriptors")

    fetch = subparsers.add_parser("fetch-data", help="Refresh and checksum-validate reference data")
    _add_measurement(fetch)

    prepare = subparsers.add_parser("prepare", help="Preflight, build the plugin, and create .run files")
    _add_campaign_options(prepare)

    campaign = subparsers.add_parser("campaign", help="Run or safely resume the physical-helicity job matrix")
    _add_campaign_options(campaign)

    postprocess = subparsers.add_parser("postprocess", help="Construct A_parallel, A1, and diagnostics")
    _add_measurement(postprocess)
    _add_tag(postprocess)
    postprocess.add_argument("--dry-run", action="store_true")

    plot = subparsers.add_parser("plot", help="Generate Rivet HTML from postprocessed outputs")
    _add_measurement(plot)
    _add_tag(plot)
    plot.add_argument("--dry-run", action="store_true")
    plot.add_argument(
        "--allow-plot-metadata-refresh",
        action="store_true",
        help=(
            "Replot a complete campaign only when its immutable generation "
            "signature can be reconstructed by restoring the historical "
            "Rivet .plot file; no event or postprocess products are changed"
        ),
    )
    _add_plot_options(plot)

    full = subparsers.add_parser("full", help="Prepare, run, postprocess, and plot")
    _add_campaign_options(full)
    _add_plot_options(full)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = make_parser().parse_args(argv)
    args._tracker_started_at = time.time()
    try:
        if args.command == "list":
            registry = discover_registry()
            for identifier, measurement in registry.items():
                print(f"{identifier}\t{measurement.get('title', '')}")
            return 0

        measurement = get_measurement(args.measurement)
        if args.command == "fetch-data":
            archive = fetch_reference_data(measurement)
            print(f"Verified {measurement['id']} source data and refreshed reference YODA; cached archive: {archive}")
        elif args.command == "prepare":
            prepare_campaign(args, measurement)
        elif args.command == "campaign":
            run_campaign(args, measurement)
        elif args.command == "postprocess":
            postprocess_campaign(args, measurement)
        elif args.command == "plot":
            plot_campaign(args, measurement)
        elif args.command == "full":
            prepare_campaign(args, measurement)
            if not args.dry_run:
                run_campaign(args, measurement)
                postprocess_campaign(args, measurement)
                plot_campaign(args, measurement)
        else:
            raise CampaignError(f"Unsupported command {args.command}")
    except CampaignError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
