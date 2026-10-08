"""Additive MC/data figures for the published HERMES Born cells.

The existing absolute atlas remains unchanged. This module reuses its checked
snapshot and renders additional fixed-x comparisons and ratio-only atlases.
"""
from __future__ import annotations

import copy
from pathlib import Path
from typing import Any, Mapping

import hermes_born_cell_plots as born
import hermes_ratio_common as ratio

SETTINGS = {
    "figure_mm": [180, 160], "single_figure_mm": [100, 112],
    "rows": 4, "columns": 5, "font_pt": 8, "png_dpi": 300,
    "focus_ratio_limits": [-1., 3.], "q2_limits": [0.18, 20.],
    "full_range_scales": "independent panel scales, including all supported uncertainties",
    "theory_color": born.SETTINGS["theory_color"],
}
FIGURE_STEMS = {target: f"AParallel_{target}_BornCells_Ratio_5x4" for target in ("P", "D")}
REQUIRED_OUTPUTS = {
    f"{stem}{extra}.{suffix}" for stem in FIGURE_STEMS.values()
    for extra in ("", "_FullRange") for suffix in ("pdf", "png")
} | {
    f"AParallel_{target}_X{index:02d}_WithRatio.{suffix}"
    for target in ("P", "D") for index in range(1, 20) for suffix in ("pdf", "png")
}


def adapt_rows(cells: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Preserve each measured cell and join ratio diagnostics without reweighting."""
    rows = []
    for cell in cells:
        row = {
            "low": cell["q2_low"], "high": cell["q2_high"], "marker": cell["q2_mean"],
            "data": cell["value"], "stat": cell["stat"], "total": cell["total_error"],
            "mc": cell["mc_a_parallel"], "mc_stat": cell["mc_stat"],
            "supported": cell["supported_prediction"], "global_bin": cell["global_bin"],
            "selection": cell["selection"], "bin": cell["bin"],
        }
        row.update(ratio.ratio_record(row["data"], row["stat"], row["total"],
                                      row["mc"], row["mc_stat"], supported=row["supported"]))
        rows.append(row)
    return rows


def build_snapshot(summary: dict[str, Any], reference: dict[str, Any]) -> dict[str, Any]:
    """Build the existing checked absolute snapshot plus additive ratio metadata."""
    return with_ratio_metadata(born.build_snapshot(summary, reference))


def with_ratio_metadata(snapshot: Mapping[str, Any]) -> dict[str, Any]:
    result = copy.deepcopy(snapshot)
    result["ratio_layout"] = copy.deepcopy(SETTINGS)
    result["ratio_interpretation"] = [
        "MC/data is retained for every nonzero measured central value with a supported nominal prediction.",
        "Orange ratio uncertainties are MC statistical only; gray bands about one show relative data statistical and total uncertainties.",
        "Open markers flag measured central values compatible with zero within their total uncertainty; these ratios are unstable diagnostics.",
        "Focused figures use -1 to 3 with explicitly labelled arrows for cropped central values; full-range atlases show all supported uncertainties on independent panel scales.",
        "Published mean Q2 positions, native cell boundaries, low-Q2 masks, and target normalization remain unchanged.",
    ]
    for block in result["targets"].values():
        for panel in block["panels"]:
            panel["ratio_rows"] = adapt_rows(panel["cells"])
            panel["full_range_ratio_limits"] = list(ratio.ratio_limits(panel["ratio_rows"], focus=False))
    return result


def _ratio_axis(axis: Any, panel: dict[str, Any], *, full_range: bool = False,
                single: bool = False) -> None:
    from matplotlib.ticker import FixedLocator, FuncFormatter, LogLocator, MaxNLocator, NullFormatter

    axis.set_xscale("log")
    axis.set_xlim(*SETTINGS["q2_limits"])
    axis.axvspan(SETTINGS["q2_limits"][0], 1., facecolor="0.93", edgecolor="none", zorder=0)
    rows = panel["ratio_rows"]
    limits = panel["full_range_ratio_limits"] if full_range else SETTINGS["focus_ratio_limits"]
    ratio.draw_ratio(axis, rows, ylim=limits)
    axis.xaxis.set_major_locator(FixedLocator([.2, 1., 10.]))
    axis.xaxis.set_major_formatter(FuncFormatter(lambda value, position: f"{value:g}"))
    axis.xaxis.set_minor_locator(LogLocator(base=10, subs=(2., 3., 4., 5., 6., 7., 8., 9.)))
    axis.xaxis.set_minor_formatter(NullFormatter())
    axis.yaxis.set_major_locator(MaxNLocator(nbins=3, min_n_ticks=3) if full_range
                                 else FixedLocator([-1., 0., 1., 2., 3.]))
    axis.tick_params(axis="both", which="major", labelsize=7 if full_range else 8,
                     length=3, width=.6, pad=2)
    axis.tick_params(axis="x", which="minor", length=1.5, width=.5)
    axis.tick_params(labelbottom=single or panel["row_index"] == 3,
                     labelleft=single or full_range or panel["column_index"] == 0)
    for spine in axis.spines.values():
        spine.set_linewidth(.6)
    if not single:
        axis.set_title(rf"${panel['display_x_low']:g}<x<{panel['display_x_high']:g}$",
                       pad=4, fontsize=8)
    if not any(row["ratio"] is not None for row in rows):
        reason = "No MC support" if not any(row["supported"] for row in rows) else "Ratio undefined"
        axis.text(.5, .5, reason, transform=axis.transAxes, ha="center", va="center", fontsize=7,
                  color="0.4", bbox={"facecolor": "white", "edgecolor": "none", "alpha": .8, "pad": 1})


def _save(figure: Any, output: Path, stem: str) -> list[Path]:
    figure.canvas.draw()
    paths = []
    for suffix in ("pdf", "png"):
        path = output/f"{stem}.{suffix}"
        figure.savefig(path)
        paths.append(path)
    return paths


def _legend(axis: Any) -> None:
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    axis.set_axis_off()
    axis.legend(handles=[
        Patch(facecolor=SETTINGS["theory_color"], alpha=.22,
              edgecolor=SETTINGS["theory_color"], label="Herwig / HERMES\nband: MC stat."),
        Patch(facecolor="0.75", edgecolor="none", label="Data stat. / |data|"),
        Patch(facecolor="0.9", edgecolor="none", label="Data total / |data|"),
        Line2D([], [], color=SETTINGS["theory_color"], marker="o", markerfacecolor="white",
               linewidth=0, markersize=4, label="Data compatible\nwith zero"),
    ], loc="center", frameon=False, fontsize=7, handlelength=.9,
        handletextpad=.5, labelspacing=.8, borderaxespad=0)


def _render_atlas(snapshot: dict[str, Any], target: str, output: Path,
                  *, full_range: bool = False) -> list[Path]:
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(4, 5, figsize=tuple(mm/25.4 for mm in SETTINGS["figure_mm"]))
    figure.subplots_adjust(left=.095, right=.985, bottom=.115, top=.895,
                           hspace=.40, wspace=.42 if full_range else .16)
    for panel in snapshot["targets"][target]["panels"]:
        _ratio_axis(axes.flat[panel["panel_index"]], panel, full_range=full_range)
    _legend(axes.flat[19])
    figure.suptitle(rf"HERMES Born $A_\parallel$: {snapshot['targets'][target]['label']} MC/data",
                   y=.978, fontsize=10)
    subtitle = ("Full range; independent panel scales" if full_range
                else "Focused range; arrows label off-scale central values")
    figure.text(.54, .937, subtitle, ha="center", fontsize=8)
    figure.supxlabel(r"$Q^2$ [GeV$^2$]", x=.54, y=.042, fontsize=10)
    figure.supylabel("Herwig / HERMES", x=.018, y=.51, fontsize=9)
    figure.text(.54, .014, r"Gray $Q^2<1$ GeV$^2$: no prediction. Data bands are relative to unity.",
                ha="center", fontsize=7)
    stem = FIGURE_STEMS[target]+("_FullRange" if full_range else "")
    try:
        return _save(figure, output, stem)
    finally:
        plt.close(figure)


def _render_single(snapshot: dict[str, Any], target: str, panel: dict[str, Any],
                   output: Path) -> list[Path]:
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    figure, (upper, lower) = plt.subplots(2, 1, sharex=True,
        figsize=tuple(mm/25.4 for mm in SETTINGS["single_figure_mm"]),
        gridspec_kw={"height_ratios": [2., 1.]})
    figure.subplots_adjust(left=.18, right=.97, top=.855, bottom=.225, hspace=.07)
    born._panel(upper, panel, snapshot["row_y_limits"][panel["row_index"]])
    upper.tick_params(labelleft=True, labelbottom=False)
    upper.set_ylabel(rf"$A_\parallel^{{{'p' if target == 'P' else 'd'}}}$", fontsize=9)
    _ratio_axis(lower, panel, single=True)
    lower.set_ylabel("MC/data", fontsize=8)
    lower.set_xlabel(r"$Q^2$ [GeV$^2$]", fontsize=9)
    figure.suptitle(rf"HERMES Born $A_\parallel$: {snapshot['targets'][target]['label']}",
                   y=.98, fontsize=10)
    figure.legend(handles=[
        Line2D([], [], color="black", marker="o", linewidth=0, markersize=3,
               label="HERMES: stat., total"),
        Patch(facecolor=SETTINGS["theory_color"], alpha=.22,
              edgecolor=SETTINGS["theory_color"], label="Herwig: MC stat."),
    ], loc="upper center", bbox_to_anchor=(.57, .943), frameon=False, fontsize=6.7,
        ncol=2, handlelength=1.1, columnspacing=.8, labelspacing=.3)
    figure.text(.57, .060, "Ratio: MC stat. only; gray bands: data stat./total.\n"
                "Open markers: data compatible with zero.\n"
                "Arrows label off-scale values; full-range atlas available.",
                ha="center", va="center", fontsize=6.5, linespacing=1.5)
    stem = f"AParallel_{target}_X{panel['panel_index']+1:02d}_WithRatio"
    try:
        return _save(figure, output, stem)
    finally:
        plt.close(figure)


def render(snapshot: dict[str, Any], output: Path) -> list[Path]:
    """Write 38 new absolute/ratio figures and four ratio atlases, PDF and PNG."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    # Accept either the existing checked Born snapshot or this module's richer
    # snapshot. In both cases no published value or primary file is mutated.
    snapshot = with_ratio_metadata(snapshot)
    paths = []
    with plt.rc_context({"text.usetex": False, "font.family": "DejaVu Sans", "font.size": 8,
                         "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 300,
                         "savefig.bbox": None}):
        for target in ("P", "D"):
            for panel in snapshot["targets"][target]["panels"]:
                paths.extend(_render_single(snapshot, target, panel, output))
            for full_range in (False, True):
                paths.extend(_render_atlas(snapshot, target, output, full_range=full_range))
    if {path.name for path in paths} != REQUIRED_OUTPUTS:
        raise born.BornCellError("Incomplete additive Born ratio figure set")
    return paths
