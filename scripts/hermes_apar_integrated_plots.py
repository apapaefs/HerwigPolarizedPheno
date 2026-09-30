"""Render the explicitly model-assisted HERMES Born-cell projections."""

from __future__ import annotations

import html
import math
from pathlib import Path
from typing import Any


def _panel(axis: Any, rows: list[dict[str, Any]], target: str, projection: str,
           selection: str, sensitivity: bool = False) -> None:
    from matplotlib.patches import Rectangle
    from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter

    valid = [row for row in rows if row["supported"]]
    label = "proton" if target == "P" else "deuteron"
    threshold = 1 if selection == "Q2GT1" else 4
    axis.set_title(f"{label.capitalize()}, $Q^2 > {threshold}$ GeV$^2$")
    axis.set_xlabel("$x$" if projection == "x" else "$Q^2$ [GeV$^2$]")
    axis.set_ylabel(r"$A_\parallel$")
    axis.set_xscale("log")
    axis.xaxis.set_major_locator(LogLocator(base=10, subs=(1., 3.), numticks=5))
    axis.xaxis.set_major_formatter(FuncFormatter(lambda value, position: f"{value:g}"))
    axis.xaxis.set_minor_locator(LogLocator(base=10, subs=(2., 4., 5., 6., 7., 8., 9.)))
    axis.xaxis.set_minor_formatter(NullFormatter())
    axis.axhline(0, color="0.7", linewidth=0.7, zorder=0)
    axis.grid(axis="y", alpha=0.15)
    if rows:
        axis.set_xlim(min(row["bin_low"] for row in rows),
                      max(row["bin_high"] for row in rows))
    if not valid:
        axis.text(0.5, 0.5, "No supported bins", transform=axis.transAxes,
                  ha="center")
        return
    x = [row["marker"] for row in valid]
    values = [row["a_parallel"] for row in valid]
    xerrors = [[row["marker"] - row["bin_low"] for row in valid],
               [row["bin_high"] - row["marker"] for row in valid]]
    outer = [row["stat_plus_systematic_upper_bound"] for row in valid]
    axis.errorbar(x, values, yerr=outer, fmt="none", ecolor="0.55",
                  elinewidth=2.7, alpha=0.65,
                  label="Statistical + systematic bound", zorder=3)
    axis.errorbar(x, values, xerr=xerrors, yerr=[row["stat"] for row in valid],
                  fmt="o", color="black", markersize=3.5, capsize=2,
                  linewidth=0.9, label="HERMES Born cells, GD11 weighted", zorder=4)
    first = True
    for row in valid:
        value = row.get("mc_a_parallel")
        if value is None or not math.isfinite(value):
            continue
        error = row.get("mc_stat_independent_cell_approx")
        axis.hlines(value, row["bin_low"], row["bin_high"], color="#d55e00",
                    linewidth=1.4, label="Herwig cells, same GD11 weights" if first else None,
                    zorder=2)
        if error is not None and math.isfinite(error) and error > 0:
            axis.add_patch(Rectangle((row["bin_low"], value - error),
                                     row["bin_high"] - row["bin_low"], 2 * error,
                                     facecolor="#d55e00", alpha=0.18, edgecolor="none",
                                     zorder=1))
        first = False
    if sensitivity:
        data_x, data_y, mc_x, mc_y = [], [], [], []
        for row in valid:
            control = row.get("weight_controls", {}).get("W2GT4", {})
            if control.get("data_shift") is not None:
                data_x.append(row["marker"])
                data_y.append(row["a_parallel"] + control["data_shift"])
            if control.get("mc_shift") is not None and row.get("mc_a_parallel") is not None:
                mc_x.append(row["marker"])
                mc_y.append(row["mc_a_parallel"] + control["mc_shift"])
        axis.plot(data_x, data_y, "o--", markersize=4, markerfacecolor="none",
                  color="#0072b2", linewidth=0.8, label="Data: $W^2>4$ weights (sensitivity)")
        axis.plot(mc_x, mc_y, "+--", color="#009e73", linewidth=0.8,
                  label="Herwig: $W^2>4$ weights (sensitivity)")
    axis.legend(loc="best", fontsize=7.2, frameon=False)


def render(snapshot: dict[str, Any], output: Path) -> list[Path]:
    """Write a self-contained gallery, eight panels and two overview figures."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    output = Path(output)
    paths: list[Path] = []
    panel_names: dict[tuple[str, str, str], str] = {}
    with plt.rc_context({"text.usetex": False, "font.family": "DejaVu Sans",
                         "font.size": 10, "savefig.dpi": 180}):
        for selection in ("Q2GT1", "Q2GT4"):
            for target in ("P", "D"):
                for projection in ("x", "q2"):
                    rows = [row for row in snapshot["targets"][target]["bins"]
                            if row["projection"] == projection and row["selection"] == selection]
                    name = f"AParallel_{target}_vs_{projection}_{selection}_reconstructed"
                    panel_names[target, projection, selection] = name
                    figure, axis = plt.subplots(figsize=(6.6, 4.9), layout="constrained")
                    _panel(axis, rows, target, projection, selection)
                    figure.suptitle("Reconstructed Born longitudinal asymmetry", fontsize=12)
                    assumption = "Cell-constant asymmetry where output cuts split source cells."
                    if projection == "x" and selection == "Q2GT1":
                        assumption = "Whole source cells; independent unpolarized fit weights."
                    fraction = max((row.get("fraction_of_fitted_uu_below_W2_4", 0)
                                    for row in rows if row["supported"]), default=0)
                    figure.supxlabel(assumption + "\nErrors conditional on central fit; MC cell covariance approximated."
                                      + f"\nLow-$W^2$ fit extrapolation: up to {100*fraction:.1f}% of bin weight; see sensitivity plots.",
                                      fontsize=7)
                    figure.canvas.draw()
                    for suffix in ("png", "pdf"):
                        path = output / f"{name}.{suffix}"
                        figure.savefig(path)
                        paths.append(path)
                    plt.close(figure)
            figure, axes = plt.subplots(2, 2, figsize=(12, 9), layout="constrained")
            for row_index, projection in enumerate(("x", "q2")):
                for column_index, target in enumerate(("P", "D")):
                    rows = [row for row in snapshot["targets"][target]["bins"]
                            if row["projection"] == projection and row["selection"] == selection]
                    _panel(axes[row_index, column_index], rows, target, projection, selection)
            figure.suptitle("HERMES Born data and Herwig: independent GD11-weighted cell projections")
            figure.supxlabel("Gray errors: statistical plus conservative systematic bound. "
                              "Errors conditional on the central fit.\n"
                              "Split source cells use a constant asymmetry; no extra D or deuteron factor is applied to data.",
                              fontsize=8)
            figure.canvas.draw()
            for suffix in ("png", "pdf"):
                path = output / f"integrated-overview-{selection}.{suffix}"
                figure.savefig(path)
                paths.append(path)
            plt.close(figure)
            figure, axes = plt.subplots(2, 2, figsize=(12, 9.5), layout="constrained")
            for row_index, projection in enumerate(("x", "q2")):
                for column_index, target in enumerate(("P", "D")):
                    rows = [row for row in snapshot["targets"][target]["bins"]
                            if row["projection"] == projection and row["selection"] == selection]
                    _panel(axes[row_index, column_index], rows, target, projection, selection,
                           sensitivity=True)
            figure.suptitle("Sensitivity to removing extrapolated low-$W^2$ contributions from the weights")
            figure.supxlabel("Dashed curves change weights only; source-cell asymmetries remain fixed. "
                              "These shifts are not one-sigma errors.\n"
                              "They do not reconstruct experimental or Herwig events with a new cut.", fontsize=8)
            figure.canvas.draw()
            for suffix in ("png", "pdf"):
                path = output / f"fit-domain-sensitivity-{selection}.{suffix}"
                figure.savefig(path)
                paths.append(path)
            plt.close(figure)
    sections = []
    for selection in ("Q2GT1", "Q2GT4"):
        threshold = 1 if selection == "Q2GT1" else 4
        parts = [f"<h2>Q² &gt; {threshold} GeV²</h2>",
                 f'<p><a href="integrated-overview-{selection}.pdf">Overview PDF</a> · '
                 f'<a href="fit-domain-sensitivity-{selection}.pdf">Fit-domain sensitivity PDF</a></p>']
        for projection in ("x", "q2"):
            for target in ("P", "D"):
                name = panel_names[target, projection, selection]
                label = ("Proton" if target == "P" else "Deuteron") + " versus " + ("x" if projection == "x" else "Q²")
                parts.append(f'<article><h3>{html.escape(label)}</h3><a href="{name}.pdf">'
                             f'<img src="{name}.png" alt="{html.escape(label)}"></a></article>')
        sections.append("\n".join(parts))
    domain_rows = []
    for target in ("P", "D"):
        for projection in ("x", "q2"):
            rows = [row for row in snapshot["targets"][target]["bins"]
                    if row["projection"] == projection and row["selection"] == "Q2GT1" and row["supported"]]
            fraction = max((row.get("fraction_of_fitted_uu_below_W2_4", 0) for row in rows), default=0)
            shift = max((abs(row.get("weight_controls", {}).get("W2GT4", {}).get("data_shift") or 0)
                         for row in rows), default=0)
            domain_rows.append(f"<li>{'Proton' if target == 'P' else 'Deuteron'} versus "
                               f"{'x' if projection == 'x' else 'Q²'}: up to {100*fraction:.1f}% "
                               f"low-W² weight; largest absolute data shift {shift:.4f}.</li>")
    domain = ("<h2>Fit-domain check</h2><p>GD11 is fitted for W² &gt; 4 GeV²; the analysis extends to W² &gt; 3.24 GeV². "
              "The central weights therefore extrapolate the fit. The separate sensitivity plots remove that contribution "
              "from the weights while keeping the same measured cells. This check is especially important at high x.</p><ul>"
              + "".join(domain_rows) + "</ul>")
    page = """<!doctype html><html lang="en"><meta charset="utf-8">
<title>Reconstructed HERMES integrated Born asymmetries</title>
<style>body{font:16px system-ui,sans-serif;max-width:1200px;margin:2rem auto;padding:0 1rem;color:#222}
article{display:inline-block;width:49%;vertical-align:top}img{width:100%}p{line-height:1.5}
@media(max-width:700px){article{width:100%}}</style>
<h1>Reconstructed HERMES integrated Born asymmetries</h1>
<p>Published proton and deuteron Born cells projected using independent GD11 unpolarized fits.
These are model-assisted reconstructions. The Herwig comparison uses the same weights on its cell results.</p>
<p>Statistical errors use the complete experimental covariance. Gray errors add a conservative
bound for unknown systematic correlations; the published normalization is already included.
Errors are conditional on the central fit. The x-integrated Q² projection and Q²&gt;4 x controls
assume a constant asymmetry within split source cells. Fit/acceptance sensitivity shifts are recorded
separately in the numerical output and are not one-sigma uncertainties.</p>
<p><a href="integrated.csv">Numerical CSV</a> · <a href="integrated.json">Weights, covariance and provenance (JSON)</a></p>
""" + domain + "\n".join(sections) + "</html>\n"
    index = output / "index.html"
    index.write_text(page, encoding="utf-8")
    paths.append(index)
    return paths
