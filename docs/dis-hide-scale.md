# Plot DIS results without the hard-scale envelope

`--hide-scale-uncertainty` is a plotting-only option of
`scripts/run_experimental_campaign.py`, the fixed-target DIS route in
`scripts/run_phenomenology_campaign.py`, and the E142, E154, and HERMES 1997 wrappers.
It removes only the named `hard_scale` error from a derived plotting input,
retaining MC statistical errors, any PDF errors, central values, masks, and
experimental uncertainties. Both main and residual panels use this input.
The separate scale overlay is also suppressed; MC and PDF overlays are retained.

This matches the uncertainty content of a central-scale campaign with otherwise
identical PDF/family settings. It uses the central events of the scale study,
so its central values need not numerically match an independent production run.
It does not merge events from the different scale points.

From the repository root on the production host, using the activated Herwig
environment (replace `/path/to/work` and the tag with your existing campaign):

```bash
python3 scripts/run_E142_1996_I424108.py \
  --output-root /path/to/work \
  plot --tag e142-scale-study \
  --hide-scale-uncertainty --allow-plot-metadata-refresh
```

For E154 substitute `run_E154_1997_I443170.py` and `e154-scale-study`;
for HERMES substitute `run_HERMES_1997_I440904.py` and `hermes1997-scale-study`.
`--dry-run` prints the plan without writing anything.

Cleaner plots go into `plots-central/`, with separate logs and a `plots_central`
manifest record. The derived YODA is stored under `work/plot-central/`.
Existing `plots/`, analyzed YODA, summaries, and generated events are preserved.
Repeating the option updates its own output directory. Omitting it retains the
existing full-uncertainty plotting behavior. It is also accepted by `full`, but
use `plot` to avoid requesting generation/postprocessing. Q2-window selections
and optional comparison families are unaffected.

These display controls apply to fixed-target DIS. The generic campaign CLI
rejects them for other process kinds rather than silently ignoring them.

## Plot one scale prediction with its own MC uncertainty

Use `--plot-scale 0.5`, `--plot-scale 1`, or `--plot-scale 2` to display the
central-PDF prediction at that scale with **its own MC statistical uncertainty**.
This reads the saved analyzed/variations YODA; no regeneration or reprocessing
is needed for an already postprocessed scale study. In particular, the half-
and double-scale predictions do not inherit the central-scale MC errors.

The main plot shows the selected prediction and its MC statistical band;
the residual panel also uses that prediction and its statistical errors.
The label identifies `mu/mu0`. There is no scale envelope or PDF band in
this mode. `--plot-scale 1` therefore differs from `--hide-scale-uncertainty`
if the campaign contains PDF uncertainties: the latter retains the PDF band.
These options cannot be combined. `--plot-comparisons` is also incompatible
because comparison families are not generated at all three scales.

The output directories are `plots-mu0p5/`, `plots-mu1/`, and `plots-mu2/`,
with separate logs and manifest records. Derived YODA inputs go under
`work/plot-mu0p5/`, etc. The original full plots, central-only plots,
postprocessing, and event samples are preserved. Repeating a selection
updates its own output directory.

For example, to generate half- and double-scale plots for existing E154 and
HERMES 1997 scale studies:

```bash
for scale in 0.5 2; do
  python3 scripts/run_E154_1997_I443170.py \
    --output-root /path/to/work \
    plot --tag e154-scale-study --plot-scale "$scale" || break
  python3 scripts/run_HERMES_1997_I440904.py \
    --output-root /path/to/work \
    plot --tag hermes1997-scale-study --plot-scale "$scale" || break
done
```

Add `--dry-run` to preview without writing files. Use `1` as an additional
loop value to create a third set of statistical-only plots. The E142 wrapper
accepts the same option, selecting both helium and neutron outputs at the
requested scale. The option is also available on `full`; use `plot` for
existing campaigns to avoid requesting generation/postprocessing.

Only the configured asymmetry and parity outputs are selected. Accepted-
kinematics diagnostics are not copied from the central sample and relabeled
as another scale. Masked/invalid bins in a varied-scale prediction remain
masked; no nominal-scale values are substituted. A missing scale sample or
missing variation object is an error, rather than a fallback to the central
prediction. For old E142 campaigns, the existing output-schema compatibility
check still applies; any required schema reprocessing must precede plotting.

HERMES 2007 supplemental galleries read nominal/full-uncertainty summaries.
The alternate display modes therefore omit those galleries, print an explicit
notice, and record the omission in their manifest entry. The default `plots/`
view retains the reconstructed Born, cell, acceptance-control, and data-ratio
galleries. Summary-only output selections are skipped when selecting YODA
objects for a single scale.

Run the focused tests in the activated environment:

```bash
python3 -m unittest discover -s scripts/tests -p test_dis_hide_scale_plot_nbaldona.py -v
python3 -m unittest discover -s scripts/tests -p test_dis_plot_scale_nbaldona.py -v
```
