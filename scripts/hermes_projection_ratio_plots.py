"""Additive A1 and matched GD11 projection figures with MC/data panels.

The reconstructed comparisons deliberately use the GD11-weighted MC cell
projection, not the separate ratio of event-level integrated cross sections.
Existing absolute-only renderers and campaign products are not modified.
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Mapping, Sequence

MEASUREMENT = "HERMES_2007_I726689"
MC_COLOR = "#d55e00"
# Preserve the original d14 plot metadata; ratio panels add no new physics.
A1_PRESENTATION = {
    target: {
        "rivet_id": f"d14-x01-y0{index}",
        "title": rf"HERMES {label} virtual-photon asymmetry ($Q^2>1\,\mathrm{{GeV}}^2$)",
        "ylabel": rf"$A_1^{{{target.lower()}}}$ ($\eta A_2$ neglected in prediction)",
        "absolute_ylim": (-.15, 1.50),
        "theory_label": "NLO+PS (polarized; full spin)",
    }
    for index, (target, label) in enumerate((("P", "proton"), ("D", "deuteron")), 1)
}
COMPARISON_IDS = ("A1_P_Q2GT1", "A1_D_Q2GT1") + tuple(
    f"AParallel_{target}_vs_{axis}_{selection}_reconstructed"
    for target in ("P", "D") for selection in ("Q2GT1", "Q2GT4") for axis in ("x", "q2"))
REQUIRED_OUTPUTS = {f"{name}{variant}.{extension}" for name in COMPARISON_IDS
                    for variant in ("_ratio", "_ratio_fullrange") for extension in ("pdf", "png")}


def _finite(value: Any) -> bool:
    return isinstance(value, (int, float)) and math.isfinite(value)


def _row(*, low: float, high: float, marker: float | None, data: float | None,
         stat: float | None, total: float | None, mc: float | None,
         mc_stat: float | None, supported: bool, **extra: Any) -> dict[str, Any]:
    from hermes_ratio_common import ratio_record

    if not all(_finite(value) for value in (low, high)) or not 0 < low < high:
        raise ValueError("Comparison bins require finite positive ordered boundaries")
    marker = math.sqrt(low * high) if marker is None else marker
    if not _finite(marker) or not low <= marker <= high:
        raise ValueError("Comparison marker must lie within its physical bin")
    row = dict(low=low, high=high, marker=marker, data=data, stat=stat,
               total=total, mc=mc, mc_stat=mc_stat, supported=bool(supported), **extra)
    row.update(ratio_record(data, stat, total, mc, mc_stat, supported=bool(supported)))
    return row


def _a1_total_error(point: Mapping[str, Any]) -> float:
    """Match the original Rivet A1 errors before rounding the combined systematic.

    Older reference fixtures without component fields can use the published
    rounded combined value. A partially specified component mapping is invalid.
    """
    components = point.get("systematics")
    if components is None:
        return math.hypot(point["stat"], point["systematic_combined"])
    required = ("experimental", "parameterization", "evolution")
    if any(name not in components for name in required):
        raise ValueError("A1 systematic components must include experimental, parameterization and evolution")
    return math.sqrt(point["stat"]**2 + sum(components[name]**2 for name in required))


def _check_identities(summary: Mapping[str, Any], reference: Mapping[str, Any],
                      integrated: Mapping[str, Any]) -> None:
    for document in (summary, reference, integrated):
        if document.get("measurement", MEASUREMENT) != MEASUREMENT:
            raise ValueError("Ratio comparisons require HERMES_2007_I726689 inputs")
    if (summary.get("tag") and integrated.get("campaign_tag")
            and summary["tag"] != integrated["campaign_tag"]):
        raise ValueError("Summary and reconstructed projections belong to different campaigns")
    acceptance = integrated.get("fit", {}).get("cuts", {}).get("acceptance")
    if acceptance and summary.get("acceptance") and acceptance != summary["acceptance"]:
        raise ValueError("Summary and reconstructed projections have different acceptances")


def build_comparisons(summary: Mapping[str, Any], referenceA1: Mapping[str, Any],
                      integrated: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Return two published A1 and eight matched reconstructed comparisons.

    Each standardized row retains the original bin/marker and separate data
    and MC errors. Experimental uncertainty is displayed around unity in the
    ratio panel; only the MC statistical error accompanies MC/data itself.
    """
    _check_identities(summary, referenceA1, integrated)
    datasets = {dataset["target"]: dataset for dataset in referenceA1["datasets"]}
    if len(datasets) != len(referenceA1["datasets"]) or set(datasets) != {"P", "D"}:
        raise ValueError("Expected one proton and one deuteron A1 reference")
    comparisons: list[dict[str, Any]] = []
    ratio_notes = [
        "MC/data errors contain MC statistics only; gray bands at unity show the separate data statistical and total errors.",
        "Nonzero data compatible with zero retain their ratios and are flagged; an exactly zero denominator has no defined ratio.",
        "Focused ratio range is [-1, 3]; full-range companions retain every finite ratio and uncertainty extent.",
    ]
    for target, label, selection in (("P", "Proton", "Q2GT1"),
                                      ("D", "Deuteron", "D_Q2GT1")):
        reference = datasets[target]
        edges, points = reference["bin_edges"], reference["points"]
        native = [row for row in summary["bins"] if row["selection"] == selection]
        by_bin = {row["bin"]: row for row in native}
        if len(by_bin) != len(native) or set(by_bin) != {point["bin"] for point in points}:
            raise ValueError(f"A1 {target} prediction and reference bins do not match")
        if len(edges) != len(points) + 1:
            raise ValueError(f"A1 {target} reference has inconsistent bin edges")
        rows = []
        for point in sorted(points, key=lambda row: row["bin"]):
            number = point["bin"]
            predicted = by_bin[number]
            low, high = edges[number - 1:number + 1]
            if not (math.isclose(predicted["bin_low"], low, rel_tol=0., abs_tol=1e-10)
                    and math.isclose(predicted["bin_high"], high, rel_tol=0., abs_tol=1e-10)):
                raise ValueError(f"A1 {target} bin {number} boundaries do not match")
            rows.append(_row(low=low, high=high, marker=point["x_mean"],
                             data=point["value"], stat=point["stat"],
                             total=_a1_total_error(point),
                             mc=predicted.get("a1"), mc_stat=predicted.get("a1_stat"),
                             supported=predicted.get("supported", True),
                             bin=number, target=target, selection=selection))
        comparisons.append(dict(
            id=f"A1_{target}_Q2GT1", **A1_PRESENTATION[target],
            observable="A1", axis="x", kind="a1", rows=rows,
            notes=["Published Q2-averaged A1 points from HERMES Table XXII; no Q2>4 A1 data comparison is constructed.",
                   "Herwig is the inverse-D proxy and omits eta*A2; its event-weighted averaging differs from the experimental prescription.",
                   "The original Rivet title, target label and prediction legend are retained; the common absolute range is extended from [-0.15, 1.45] to [-0.15, 1.50] to show the complete final proton data error.",
                   "Data total errors combine statistics and the three published systematic components in quadrature, matching the original Rivet panels; normalization is already included. References without components use the rounded combined systematic as a fallback.",
                   *ratio_notes]))
    for target, label in (("P", "Proton"), ("D", "Deuteron")):
        for selection in ("Q2GT1", "Q2GT4"):
            for axis in ("x", "q2"):
                source = sorted((row for row in integrated["targets"][target]["bins"]
                                 if row["projection"] == axis and row["selection"] == selection),
                                key=lambda row: row["bin"])
                if not source or len({row["bin"] for row in source}) != len(source):
                    raise ValueError(f"Missing or duplicate reconstructed bins: {target}/{selection}/{axis}")
                rows = [_row(low=row["bin_low"], high=row["bin_high"], marker=row["marker"],
                             data=row["a_parallel"], stat=row["stat"],
                             total=row["stat_plus_systematic_upper_bound"],
                             mc=row.get("mc_a_parallel"),
                             mc_stat=row.get("mc_stat_independent_cell_approx"),
                             supported=row["supported"], bin=row["bin"], target=target,
                             selection=selection, source_id=row["id"],
                             systematic_upper_bound=row["systematic_upper_bound"],
                             mc_supported=row.get("mc_supported", row.get("mc_a_parallel") is not None))
                        for row in source]
                threshold = 1 if selection == "Q2GT1" else 4
                comparisons.append(dict(
                    id=f"AParallel_{target}_vs_{axis}_{selection}_reconstructed",
                    title=label + rf": $Q^2>{threshold}$ GeV$^2$",
                    observable="A_parallel", axis=axis, kind="reconstructed", rows=rows,
                    notes=["Model-assisted Born A_parallel projection: identical independent GD11 weights multiply the published and Herwig cell asymmetries.",
                           "The MC values are mc_a_parallel from the reconstruction; event-level integrated A_parallel is a different estimator and is not overlaid.",
                           "Experimental statistical covariance is propagated by the reconstruction. The outer error combines statistics with a conservative systematic correlation bound, not a recovered systematic covariance.",
                           "Errors are conditional on central fit weights; split source cells assume constant asymmetry. MC errors use the independent-cell approximation.",
                           *ratio_notes]))
    return comparisons


def _draw_absolute(axis: Any, comparison: Mapping[str, Any]) -> None:
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    rows = comparison["rows"]
    data = [row for row in rows if row["supported"] and _finite(row["data"])]
    for row in rows:
        if not row["supported"] or not _finite(row["mc"]):
            continue
        low, high, value, error = row["low"], row["high"], row["mc"], row["mc_stat"]
        if _finite(error):
            axis.fill_between([low, high], [value-error]*2, [value+error]*2,
                              color=MC_COLOR, alpha=.2, linewidth=0, zorder=1)
        axis.hlines(value, low, high, color=MC_COLOR, linewidth=1.05, zorder=2)
    for row in data:
        x, y = row["marker"], row["data"]
        if _finite(row["total"]):
            axis.errorbar(x, y, yerr=row["total"], fmt="none", color=".65",
                          elinewidth=2.2, capsize=0, zorder=3)
        if _finite(row["stat"]):
            axis.errorbar(x, y, yerr=row["stat"], fmt="o", markersize=2.8,
                          xerr=[[x-row["low"]], [row["high"]-x]], color="black",
                          linewidth=.65, capsize=1.4, zorder=4)
    reconstructed = comparison["kind"] == "reconstructed"
    if reconstructed:
        handles = [
            Line2D([], [], color="black", marker="o", linestyle="none", markersize=3,
                   label="HERMES (GD11 weighted)"),
            Line2D([], [], color=".65", linewidth=2.2, label="Stat. + syst. bound"),
            Patch(facecolor=MC_COLOR, alpha=.35, edgecolor=MC_COLOR,
                  label="Herwig (same weights)"),
        ]
    else:
        handles = [
            Line2D([], [], color="black", marker="o", linestyle="none", markersize=3,
                   label="Data"),
            Line2D([], [], color=".65", linewidth=2.2, label="Data stat. + syst."),
            Line2D([], [], color=MC_COLOR, linewidth=1.05,
                   label=comparison["theory_label"]),
            Patch(facecolor=MC_COLOR, alpha=.2, edgecolor="none", label="MC statistical"),
        ]
    axis.legend(handles=handles, loc="best" if reconstructed else "upper left",
                frameon=False, fontsize=7. if reconstructed else 8.,
                handlelength=1.5, labelspacing=.25)
    axis.set_ylabel(r"$A_\parallel$" if reconstructed else comparison["ylabel"])
    axis.set_title(comparison["title"], fontsize=8.5 if reconstructed else 9., pad=5)
    axis.axhline(0., color=".7", linewidth=.5, zorder=0)
    axis.grid(axis="y", color=".91", linewidth=.45)
    axis.margins(y=.15)
    if not reconstructed:
        axis.set_ylim(comparison["absolute_ylim"])


def render(comparisons: Sequence[Mapping[str, Any]], output: Path) -> list[Path]:
    """Write focus/full-range PDF and 300-dpi PNG versions of each comparison."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter
    from hermes_ratio_common import draw_ratio, ratio_limits

    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    style = {"font.family": "DejaVu Sans", "font.size": 8.5, "text.usetex": False,
             "axes.labelsize": 8.5, "xtick.labelsize": 8., "ytick.labelsize": 8.,
             "axes.linewidth": .65, "pdf.fonttype": 42, "ps.fonttype": 42,
             "savefig.dpi": 300, "xtick.direction": "in", "ytick.direction": "in"}
    with plt.rc_context(style):
        for comparison in comparisons:
            rows = comparison["rows"]
            if not rows:
                raise ValueError("Cannot render a comparison without bins")
            # The original A1 caveat and full prediction label need more space
            # than the concise model-assisted projection captions.
            width_mm, height_mm = (140, 125) if comparison["kind"] == "a1" else (90, 98)
            for focus in (True, False):
                figure, (absolute, ratio) = plt.subplots(
                    2, 1, sharex=True, figsize=(width_mm/25.4, height_mm/25.4),
                    gridspec_kw={"height_ratios": (2.4, 1.2)}, layout="constrained")
                try:
                    _draw_absolute(absolute, comparison)
                    draw_ratio(ratio, rows, ylim=ratio_limits(rows, focus=focus))
                    ratio.set_ylabel("MC / data")
                    ratio.set_xlabel("$x$" if comparison["axis"] == "x" else r"$Q^2$ [GeV$^2$]")
                    for axis in (absolute, ratio):
                        axis.set_xscale("log")
                        axis.set_xlim(min(row["low"] for row in rows), max(row["high"] for row in rows))
                        axis.tick_params(which="both", top=True, right=True)
                        axis.xaxis.set_major_locator(LogLocator(base=10, subs=(1., 3.), numticks=6))
                        axis.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
                        axis.xaxis.set_minor_locator(LogLocator(base=10, subs=(2., 4., 5., 6., 7., 8., 9.)))
                        axis.xaxis.set_minor_formatter(NullFormatter())
                    suffix = "_ratio" if focus else "_ratio_fullrange"
                    for extension in ("pdf", "png"):
                        path = output / f"{comparison['id']}{suffix}.{extension}"
                        notes = "; ".join(comparison["notes"])
                        metadata = ({"Title": comparison["title"], "Subject": notes,
                                     "Creator": "HerwigPolarizedPheno HERMES ratio comparison"}
                                    if extension == "pdf" else
                                    {"Title": comparison["title"], "Description": notes})
                        figure.savefig(path, dpi=300, metadata=metadata)
                        paths.append(path)
                finally:
                    plt.close(figure)
    return paths
