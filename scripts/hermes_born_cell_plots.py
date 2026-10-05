"""Paper-sized, direct HERMES Born-cell comparisons and checked plot caches.

The published cell values are copied unchanged.  The nominal theory ratios
already contain the physical-helicity and target combinations performed by
postprocessing; this module applies no fitted weights or data corrections.
"""
from __future__ import annotations

import csv
import hashlib
import html
import json
import math
import os
from pathlib import Path
from typing import Any, Mapping


ROOT = Path(__file__).resolve().parents[1]
MEASUREMENT = "HERMES_2007_I726689"
REFERENCE = ROOT / "data/experimental" / MEASUREMENT / "born-apar-reference.json"
ESTIMATOR = "Published Born cells; direct helicity LL/UU"
CACHE_MANIFEST = "cache-manifest.json"
SETTINGS = {
    "columns": 5, "rows": 4, "panels": 19, "figure_mm": [180, 160],
    "font_pt": 8, "png_dpi": 300, "q2_limits": [0.18, 20.0],
    "row_y_padding": 0.10, "theory_family": "nominal",
    "data_outer_error": "sqrt(stat^2 + published_systematic^2)",
    "theory_error": "MC statistical only", "theory_color": "#d55e00",
}
FIGURE_STEMS = {target: f"AParallel_{target}_BornCells_5x4" for target in ("P", "D")}
REQUIRED_OUTPUTS = {"born-cells.json", "born-cells.csv", "index.html"} | {
    f"{stem}.{suffix}" for stem in FIGURE_STEMS.values() for suffix in ("pdf", "png")
}


class BornCellError(RuntimeError):
    """A direct-cell comparison cannot satisfy its input/output contract."""


def _finite(value: Any, label: str) -> float:
    if isinstance(value, bool):
        raise BornCellError(f"Invalid {label}")
    try:
        number = float(value)
    except (ValueError, TypeError) as exc:
        raise BornCellError(f"Invalid {label}") from exc
    if not math.isfinite(number):
        raise BornCellError(f"Non-finite {label}")
    return number


def _integer(value: Any, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise BornCellError(f"Invalid {label}: expected an integer")
    return value


def _boolean(value: Any, label: str) -> bool:
    if not isinstance(value, bool):
        raise BornCellError(f"Invalid {label}: expected a support boolean")
    return value


def _same(left: float, right: float) -> bool:
    return math.isclose(left, right, rel_tol=0, abs_tol=1.e-12)


def _reference_panels(reference: Mapping[str, Any], target: str) -> list[dict[str, Any]]:
    datasets = reference.get("datasets")
    if not isinstance(datasets, list):
        raise BornCellError("Born reference needs its published datasets")
    selected = {dataset.get("selection"): dataset for dataset in datasets
                if dataset.get("target") == target}
    if len(selected) != 19 or sum(dataset.get("target") == target for dataset in datasets) != 19:
        raise BornCellError(f"Target {target} must have exactly 19 fixed-x panels")
    panels = []
    for index in range(19):
        selection = f"Born{target}_X{index + 1:02d}"
        if selection not in selected:
            raise BornCellError(f"Born reference is missing {selection}")
        dataset = selected[selection]
        panel = {"selection": selection, "panel_index": index,
                 "row_index": index // 5, "column_index": index % 5}
        for key in ("x_low", "x_high", "display_x_low", "display_x_high"):
            panel[key] = _finite(dataset.get(key), f"{selection} {key}")
        if not (0 < panel["x_low"] < panel["x_high"] <= 1):
            raise BornCellError(f"Invalid x boundaries for {selection}")
        if not (0 < panel["display_x_low"] < panel["display_x_high"] <= 1):
            raise BornCellError(f"Invalid display x interval for {selection}")
        if index and not _same(panels[-1]["x_high"], panel["x_low"]):
            raise BornCellError("Published fixed-x intervals must be contiguous")
        points, edges = dataset.get("points"), dataset.get("bin_edges")
        if not isinstance(points, list) or not points or not isinstance(edges, list):
            raise BornCellError(f"Missing Born cells/physical Q2 edges for {selection}")
        if len(edges) != len(points) + 1:
            raise BornCellError(f"Born cell/edge count differs for {selection}")
        edges = [_finite(edge, f"{selection} Q2 edge") for edge in edges]
        cells = []
        for bin_index, point in enumerate(points, 1):
            cell = dict(point)
            cell.update(target=target, selection=selection,
                        panel_index=index, row_index=index // 5,
                        column_index=index % 5, **{key: panel[key] for key in
                            ("x_low", "x_high", "display_x_low", "display_x_high")})
            if _integer(cell.get("bin"), f"{selection} bin") != bin_index:
                raise BornCellError(f"Cell bin order differs for {selection}")
            _integer(cell.get("global_bin"), f"{selection} global bin")
            for key in ("value", "stat", "systematic_combined", "x_mean", "q2_mean",
                        "q2_low", "q2_high"):
                cell[key] = _finite(cell.get(key), f"{selection} bin {bin_index} {key}")
            if "y_mean" in cell:
                cell["y_mean"] = _finite(cell["y_mean"], f"{selection} y_mean")
            if not (0 < cell["q2_low"] < cell["q2_high"] <= 20):
                raise BornCellError(f"Invalid physical Q2 edges for {selection}")
            if not (_same(cell["q2_low"], edges[bin_index - 1])
                    and _same(cell["q2_high"], edges[bin_index])):
                raise BornCellError(f"Physical Q2 edges disagree for {selection}")
            if not (cell["q2_low"] <= cell["q2_mean"] <= cell["q2_high"]
                    and cell["x_low"] <= cell["x_mean"] <= cell["x_high"]):
                raise BornCellError(f"Published mean lies outside its cell for {selection}")
            if cell["stat"] < 0 or cell["systematic_combined"] < 0:
                raise BornCellError(f"Negative published uncertainty for {selection}")
            supported = _boolean(cell.get("supported_prediction"), f"{selection} support")
            if supported != (cell["q2_low"] >= 1.0) or (not supported and cell["q2_high"] > 1.0):
                raise BornCellError(f"Support mask disagrees with Q2=1 for {selection}")
            if cell["q2_low"] < SETTINGS["q2_limits"][0]:
                raise BornCellError(f"Cell lies outside the published Q2 window for {selection}")
            cell["total_error"] = math.hypot(cell["stat"], cell["systematic_combined"])
            cells.append(cell)
        panel["cells"] = cells
        panels.append(panel)
    cells = [cell for panel in panels for cell in panel["cells"]]
    if [cell["global_bin"] for cell in cells] != list(range(1, 46)):
        raise BornCellError(f"Target {target} must contain the 45 ordered published cells")
    if sum(cell["supported_prediction"] for cell in cells) != 37:
        raise BornCellError(f"Target {target} must have 37 supported and eight Q2<1 cells")
    return panels


def build_snapshot(summary: dict[str, Any], reference: dict[str, Any]) -> dict[str, Any]:
    """Join nominal normalized MC ratios to all native published Born cells."""
    if summary.get("measurement") != MEASUREMENT or reference.get("measurement") != MEASUREMENT:
        raise BornCellError(f"Direct Born-cell figures require {MEASUREMENT} inputs")
    bins = summary.get("bins")
    if not isinstance(bins, list):
        raise BornCellError("Theory summary needs Born cells; run postprocess for a Born-cell campaign")
    theory = {}
    for row in bins:
        if not str(row.get("selection", "")).startswith("Born"):
            continue
        identity = (row.get("selection"), _integer(row.get("bin"), "theory cell bin"))
        if identity in theory:
            raise BornCellError(f"Repeated theory cell {identity}")
        theory[identity] = row
    if not theory:
        raise BornCellError("Theory summary has no Born cells; older campaigns cannot provide these figures")
    targets = {}
    used = set()
    for target in ("P", "D"):
        panels = _reference_panels(reference, target)
        cells = [cell for panel in panels for cell in panel["cells"]]
        for cell in cells:
            identity = (cell["selection"], cell["bin"])
            row = theory.get(identity)
            if row is None:
                raise BornCellError(f"Theory summary is missing Born cell {identity}")
            used.add(identity)
            if _boolean(row.get("supported"), f"{identity} theory support") != cell["supported_prediction"]:
                raise BornCellError(f"Theory/reference support mismatch for {identity}")
            for key, source in (("bin_low", "q2_low"), ("bin_high", "q2_high")):
                if not _same(_finite(row.get(key), f"{identity} {key}"), cell[source]):
                    raise BornCellError(f"Theory/reference Q2 boundaries differ for {identity}")
            value, error = row.get("a_parallel"), row.get("a_parallel_stat")
            if (value is None) != (error is None):
                raise BornCellError(f"Inconsistent null theory ratio/error for {identity}")
            if value is not None:
                value = _finite(value, f"{identity} theory A_parallel")
                error = _finite(error, f"{identity} theory statistical error")
                if error < 0:
                    raise BornCellError(f"Negative theory statistical error for {identity}")
                if not cell["supported_prediction"]:
                    raise BornCellError(f"Unsupported Q2<1 theory cell must be null for {identity}")
            cell.update(mc_supported=value is not None, mc_a_parallel=value, mc_stat=error,
                        mc_missing_reason=("Q2 below generator support" if not cell["supported_prediction"]
                                           else "Missing nominal MC ratio" if value is None else None))
        targets[target] = {"label": "Proton" if target == "P" else "Deuteron",
                           "cells": cells, "panels": panels}
    if used != set(theory):
        raise BornCellError("Theory summary has unexpected Born-cell selections or bins")
    for proton, deuteron in zip(targets["P"]["panels"], targets["D"]["panels"]):
        if any(not _same(proton[key], deuteron[key]) for key in
               ("x_low", "x_high", "display_x_low", "display_x_high")):
            raise BornCellError("Proton/deuteron fixed-x intervals differ")
        if [(c["q2_low"], c["q2_high"]) for c in proton["cells"]] != [
                (c["q2_low"], c["q2_high"]) for c in deuteron["cells"]]:
            raise BornCellError("Proton/deuteron physical Q2 cells differ")
    row_y_limits = []
    for row_index in range(4):
        limits = [0.0]
        for block in targets.values():
            for cell in block["cells"]:
                if cell["row_index"] != row_index:
                    continue
                limits.extend((cell["value"] - cell["total_error"],
                               cell["value"] + cell["total_error"]))
                if cell["mc_supported"]:
                    limits.extend((cell["mc_a_parallel"] - cell["mc_stat"],
                                   cell["mc_a_parallel"] + cell["mc_stat"]))
        low, high = min(limits), max(limits)
        padding = SETTINGS["row_y_padding"] * max(high - low, 0.05)
        row_y_limits.append([low - padding, high + padding])
    return {
        "schema_version": 1, "measurement": MEASUREMENT,
        "observable": "Published Born A_parallel in native (x, Q2) cells",
        "estimator": ESTIMATOR, "theory_family": "nominal", "tag": summary.get("tag"),
        "layout": dict(SETTINGS), "row_y_limits": row_y_limits, "targets": targets,
        "reference_provenance": reference.get("provenance", {}),
        "theory_assumptions": summary.get("assumptions", []),
        "interpretation": [
            "All 45 published Born cells per target are retained unchanged at their published mean Q2.",
            "Data inner errors are statistical; outer errors add published systematics in quadrature, including normalization once.",
            "Nominal Herwig direct helicity LL/UU ratios and MC statistical errors use physical Q2 cell edges.",
            "Q2<1 cells and missing MC ratios have no theory curve; no zero or interpolation is substituted.",
            "No GD11 weights, depolarization factor or additional data normalization/nuclear correction is applied.",
        ],
    }


def _mc_groups(cells: list[dict[str, Any]]) -> list[list[dict[str, Any]]]:
    """Keep steps contiguous only where every physical source cell is present."""
    groups: list[list[dict[str, Any]]] = []
    current: list[dict[str, Any]] = []
    for cell in cells:
        if not cell["mc_supported"]:
            if current:
                groups.append(current)
                current = []
            continue
        if current and not _same(current[-1]["q2_high"], cell["q2_low"]):
            groups.append(current)
            current = []
        current.append(cell)
    if current:
        groups.append(current)
    return groups


def _panel(axis: Any, panel: dict[str, Any], y_limits: list[float]) -> None:
    import numpy as np
    from matplotlib.ticker import FixedLocator, FuncFormatter, LogLocator, MaxNLocator, NullFormatter

    cells = panel["cells"]
    axis.set_xscale("log")
    axis.set_xlim(*SETTINGS["q2_limits"])
    axis.set_ylim(*y_limits)
    axis.set_title(rf"${panel['display_x_low']:g}<x<{panel['display_x_high']:g}$", pad=4, fontsize=8)
    axis.axvspan(SETTINGS["q2_limits"][0], 1.0, facecolor="0.93", edgecolor="none", zorder=0)
    axis.axhline(0, color="0.6", linewidth=0.5, zorder=1)
    axis.xaxis.set_major_locator(FixedLocator([0.2, 1.0, 10.0]))
    axis.xaxis.set_major_formatter(FuncFormatter(lambda value, position: f"{value:g}"))
    axis.xaxis.set_minor_locator(LogLocator(base=10, subs=(2., 3., 4., 5., 6., 7., 8., 9.)))
    axis.xaxis.set_minor_formatter(NullFormatter())
    axis.yaxis.set_major_locator(MaxNLocator(nbins=3, min_n_ticks=3))
    axis.tick_params(axis="both", which="major", labelsize=8, length=3, width=0.6, pad=2)
    axis.tick_params(axis="x", which="minor", length=1.5, width=0.5)
    axis.tick_params(labelbottom=panel["row_index"] == 3,
                     labelleft=panel["column_index"] == 0)
    for spine in axis.spines.values():
        spine.set_linewidth(0.6)
    for group in _mc_groups(cells):
        edges = [group[0]["q2_low"], *(cell["q2_high"] for cell in group)]
        values = np.asarray([cell["mc_a_parallel"] for cell in group])
        errors = np.asarray([cell["mc_stat"] for cell in group])
        axis.stairs(values + errors, edges, baseline=values - errors, fill=True,
                    color=SETTINGS["theory_color"], alpha=0.22, linewidth=0, zorder=2)
        axis.stairs(values, edges, baseline=None, color=SETTINGS["theory_color"],
                    linewidth=1.15, zorder=3)
    x = [cell["q2_mean"] for cell in cells]
    values = [cell["value"] for cell in cells]
    axis.errorbar(x, values, yerr=[cell["total_error"] for cell in cells], fmt="none",
                  color="black", elinewidth=0.6, capsize=2.0, capthick=0.6, zorder=4)
    axis.errorbar(x, values, yerr=[cell["stat"] for cell in cells], fmt="o",
                  color="black", markersize=2.5, elinewidth=1.1, capsize=0, zorder=5)
    if any(cell["supported_prediction"] and not cell["mc_supported"] for cell in cells):
        axis.text(0.97, 0.96, "MC missing", transform=axis.transAxes, va="top", ha="right",
                  color=SETTINGS["theory_color"], fontsize=8)


def _write_numbers(snapshot: dict[str, Any], output: Path) -> list[Path]:
    json_path = output / "born-cells.json"
    json_path.write_text(json.dumps(snapshot, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    csv_path = output / "born-cells.csv"
    columns = ["target", "selection", "panel_index", "row_index", "column_index", "global_bin",
               "bin", "x_low", "x_high", "x_mean", "q2_low", "q2_high", "q2_mean", "value",
               "stat", "systematic_combined", "total_error", "supported_prediction", "mc_supported",
               "mc_a_parallel", "mc_stat", "mc_missing_reason"]
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for target in ("P", "D"):
            writer.writerows(snapshot["targets"][target]["cells"])
    return [json_path, csv_path]


def render(snapshot: dict[str, Any], output: Path) -> list[Path]:
    """Write the two matching 180x160 mm atlases, numbers and portable index."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    paths = _write_numbers(snapshot, output)
    with plt.rc_context({"text.usetex": False, "font.family": "DejaVu Sans",
                         "font.size": SETTINGS["font_pt"], "pdf.fonttype": 42,
                         "savefig.dpi": SETTINGS["png_dpi"]}):
        for target in ("P", "D"):
            figure, axes = plt.subplots(4, 5, figsize=tuple(mm / 25.4 for mm in SETTINGS["figure_mm"]))
            figure.subplots_adjust(left=0.10, right=0.985, bottom=0.095, top=0.90,
                                   hspace=0.35, wspace=0.14)
            for panel in snapshot["targets"][target]["panels"]:
                _panel(axes.flat[panel["panel_index"]], panel,
                       snapshot["row_y_limits"][panel["row_index"]])
            legend = axes.flat[19]
            legend.set_axis_off()
            legend.legend(handles=[
                Line2D([], [], color="black", marker="o", markersize=3, linewidth=0,
                       label="HERMES Born\ninner: stat.\nouter: total"),
                (Patch(facecolor=SETTINGS["theory_color"], alpha=0.22,
                       edgecolor=SETTINGS["theory_color"], linewidth=1.15)),
                Patch(facecolor="0.93", edgecolor="none"),
            ], labels=["HERMES Born\ninner: stat.\nouter: total",
                       "Herwig nominal\nband: MC stat.", "$Q^2<1$ GeV$^2$\ndata only"],
                loc="center", frameon=False, fontsize=8, handlelength=1.0,
                handletextpad=0.5, labelspacing=1.1, borderaxespad=0)
            figure.suptitle(rf"HERMES Born $A_\parallel$: {snapshot['targets'][target]['label']}",
                            y=0.978, fontsize=10)
            figure.text(0.54, 0.935, "Published fixed-x cells; direct helicity LL/UU",
                        ha="center", fontsize=8)
            figure.supxlabel(r"$Q^2$ [GeV$^2$]", x=0.54, y=0.027, fontsize=10)
            figure.supylabel(r"$A_\parallel$", x=0.018, y=0.51, fontsize=10)
            figure.canvas.draw()
            for suffix in ("pdf", "png"):
                path = output / f"{FIGURE_STEMS[target]}.{suffix}"
                # Do not use bbox_inches='tight': the paper's physical size is fixed.
                figure.savefig(path)
                paths.append(path)
            plt.close(figure)
    sections = []
    for target in ("P", "D"):
        block = snapshot["targets"][target]
        low = sum(not cell["supported_prediction"] for cell in block["cells"])
        missing = sum(cell["supported_prediction"] and not cell["mc_supported"] for cell in block["cells"])
        stem = FIGURE_STEMS[target]
        sections.append(f'<section><h2>{block["label"]}</h2><p>45 published cells; {low} below Q²=1 '
                        f'GeV² with data only; {missing} additional cells with an unavailable MC ratio.</p>'
                        f'<p><a href="{stem}.pdf">Paper PDF (180 × 160 mm)</a> · '
                        f'<a href="{stem}.png">PNG</a></p><a href="{stem}.pdf">'
                        f'<img src="{stem}.png" alt="{block["label"]} Born asymmetry in 19 fixed-x panels"></a></section>')
    tag = html.escape(str(snapshot.get("tag") or "unspecified"))
    page = """<!doctype html><html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>HERMES direct Born-cell asymmetries</title>
<style>body{font:16px system-ui,sans-serif;max-width:1100px;margin:2rem auto;padding:0 1rem;color:#222}
p{line-height:1.5}img{width:100%;height:auto}section{margin:2rem 0}a{color:#075b9c}</style>
<h1>HERMES direct Born-cell asymmetries</h1>
<p>All 45 published Born A_parallel cells for each target, in 19 fixed Bjorken-x panels.
Proton and deuteron figures use matching row scales. Data markers use the published mean Q²;
Herwig steps and transparent bands use the physical cell edges.</p>
<p>Black inner errors are statistical. Outer errors are the quadrature sum of statistical and
published systematic errors; the published normalization is already included. Orange curves show
the nominal direct helicity LL/UU estimator and MC statistical uncertainty only.</p>
<p>The gray Q²&lt;1 GeV² region has no generator support. Missing MC ratios are masked while all data
are retained. No continuous interpolation, fitted GD11 weighting, depolarization conversion or
extra data normalization/deuteron correction is applied.</p>
""" + f"<p>Campaign: {tag}.</p>" + """
<p><a href="born-cells.csv">Direct-cell CSV</a> · <a href="born-cells.json">Direct-cell JSON and provenance</a></p>
""" + "\n".join(sections) + "</html>\n"
    index = output / "index.html"
    index.write_text(page, encoding="utf-8")
    paths.append(index)
    return paths


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _input_hashes(summary_path: Path) -> dict[str, str]:
    return {"postprocess/summary.json": _digest(summary_path),
            REFERENCE.relative_to(ROOT).as_posix(): _digest(REFERENCE),
            "scripts/hermes_born_cell_plots.py": _digest(Path(__file__).resolve())}


def _valid_cache(directory: Path, key: str) -> bool:
    """A reusable revision must have every expected, recorded output intact."""
    try:
        manifest_path = directory / CACHE_MANIFEST
        if directory.is_symlink() or manifest_path.is_symlink():
            return False
        record = json.loads(manifest_path.read_text(encoding="utf-8"))
        outputs = record["outputs"]
        if (record["cache_key"] != key or not isinstance(outputs, dict)
                or set(outputs) != REQUIRED_OUTPUTS):
            return False
        provenance = {"inputs": record["inputs"], "settings": record["settings"]}
        if hashlib.sha256(json.dumps(provenance, sort_keys=True).encode("utf-8")).hexdigest() != key:
            return False
        if {path.name for path in directory.iterdir()} != REQUIRED_OUTPUTS | {CACHE_MANIFEST}:
            return False
        for name, digest in outputs.items():
            path = directory / name
            if path.is_symlink() or not path.is_file() or path.stat().st_size == 0 or _digest(path) != digest:
                return False
        return True
    except (OSError, ValueError, KeyError, TypeError):
        return False


def _gallery_info(campaign_dir: Path, directory: Path, key: str, *, created: bool) -> dict[str, Any]:
    def relative(path: Path) -> str:
        return Path(os.path.relpath(path, campaign_dir)).as_posix()
    return {
        "directory": relative(directory), "index": relative(directory / "index.html"),
        "cache_key": key, "created": created, "estimator": ESTIMATOR,
        "theory_family": "nominal",
        "figures": [{"target": target, "title": f"{label} Born A_parallel: 19 fixed-x panels",
                     "pdf": relative(directory / f"{FIGURE_STEMS[target]}.pdf"),
                     "png": relative(directory / f"{FIGURE_STEMS[target]}.png")}
                    for target, label in (("P", "Proton"), ("D", "Deuteron"))],
    }


def ensure_campaign_plots(campaign_dir: Path, output_dir: Path) -> dict[str, Any]:
    """Publish or reuse complete direct-cell plots inside the requested gallery.

    Failed, damaged and older revisions are retained.  A manifest is written
    last, only after every expected artifact and unchanged input is verified.
    """
    campaign_dir, output_dir = Path(campaign_dir).resolve(), Path(output_dir).resolve()
    summary_path = campaign_dir / "postprocess/summary.json"
    if not summary_path.is_file():
        raise BornCellError("Direct HERMES Born-cell plots need postprocess/summary.json; run postprocess first")
    inputs = _input_hashes(summary_path)
    provenance = {"inputs": inputs, "settings": json.loads(json.dumps(SETTINGS))}
    key = hashlib.sha256(json.dumps(provenance, sort_keys=True).encode("utf-8")).hexdigest()
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    reference = json.loads(REFERENCE.read_text(encoding="utf-8"))
    snapshot = build_snapshot(summary, reference)
    snapshot["provenance"] = {"cache_key": key, **provenance}
    cache_root = output_dir / MEASUREMENT / "born-cell-panels"
    for directory in sorted(cache_root.glob(f"{key}*")):
        if directory.is_dir() and _valid_cache(directory, key):
            if inputs != _input_hashes(summary_path):
                raise BornCellError("Born-cell inputs changed during cache validation; retry the plot stage")
            return _gallery_info(campaign_dir, directory, key, created=False)
    cache_root.mkdir(parents=True, exist_ok=True)
    attempt = 0
    while True:
        directory = cache_root / (key if attempt == 0 else f"{key}-{attempt:03d}")
        try:
            directory.mkdir()
            break
        except FileExistsError:
            attempt += 1
    paths = render(snapshot, directory)
    if {path.relative_to(directory).as_posix() for path in paths} != REQUIRED_OUTPUTS:
        raise BornCellError(f"Incomplete direct-cell HERMES gallery at {directory}")
    outputs = {}
    for name in sorted(REQUIRED_OUTPUTS):
        path = directory / name
        if path.is_symlink() or not path.is_file() or path.stat().st_size == 0:
            raise BornCellError(f"Missing or empty direct-cell output: {path}")
        outputs[name] = _digest(path)
    if inputs != _input_hashes(summary_path):
        raise BornCellError("Born-cell inputs changed during rendering; retry the plot stage")
    record = {"schema_version": 1, "cache_key": key, "estimator": ESTIMATOR,
              **provenance, "outputs": outputs}
    temporary = directory / f".{CACHE_MANIFEST}.tmp"
    temporary.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    temporary.replace(directory / CACHE_MANIFEST)
    if not _valid_cache(directory, key):
        raise BornCellError(f"Direct-cell output checksum verification failed at {directory}")
    return _gallery_info(campaign_dir, directory, key, created=True)
