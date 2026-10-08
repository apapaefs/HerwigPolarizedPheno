# Additional HERMES MC/data figures

The HERMES `plot` and `full` stages automatically add a separate ratio gallery
under `plots/HERMES_2007_I726689/data-ratios/<input-checksum>/`. Existing
absolute plots and their estimators remain available. The additional figures
are presentation products: no events, unfolding, PDF evaluation, depolarization
conversion or target correction is repeated.

The gallery contains 62 figures, each as a vector PDF with embedded fonts and
a 300-dpi PNG:

- two published A1 comparisons, proton and deuteron, with a lower ratio panel;
- eight GD11 reconstructed Born comparisons, both targets versus x and Q2,
  with Q2>1 and Q2>4 selections, each with a lower ratio panel;
- full-range companions for those ten comparisons;
- 38 separate fixed-x Born comparisons, each with absolute and ratio panels;
- separate proton/deuteron 5-column by 4-row ratio atlases, plus full-range
  companions. These retain 19 increasing-x panels and a legend tile.

The Born atlases are 180 by 160 mm. Their common logarithmic Q2 axis runs
from 0.18 to 20 GeV2; cell markers stay at the published means, while steps
and bands retain the physical cell boundaries. The eight cells per target
below Q2=1 retain their support mask. No prediction or ratio is invented.

The A1 versions identify the original Rivet plots `d14-x01-y01` (proton)
and `d14-x01-y02` (deuteron). Their visible labels retain the original HERMES
virtual-photon title, Q2 selection, target superscript, eta*A2 caveat, full
NLO+PS spin-treatment label and MC statistical legend. The absolute axes use
a common y range of -0.15 to 1.50, extending the original upper limit of
1.45 to retain the entire last proton error bar. Both figures are embedded
directly in the complete
gallery, with PDF, PNG, full-range and original-plot links; the original
absolute figures remain intact.

## Ratio and uncertainty convention

For a prediction M and experimental central value d, the plotted ratio is
M/d. The orange uncertainty is sigma_MC/abs(d), holding the observed denominator
fixed. Dark and light gray bands about one display sigma_data_stat/abs(d)
and sigma_data_outer/abs(d), separately. These displays are not a Gaussian
confidence interval for a quotient with an uncertain denominator, and the
experimental and MC errors are not combined into one ratio error.

The direct Born outer experimental errors combine statistical and published
combined systematic errors in quadrature. The A1 figures reproduce the original
Rivet errors by combining the statistical, experimental, parameterization and
evolution components in quadrature, avoiding the rounding in the separately
printed combined systematic column. Normalization is already included. In the
model-assisted reconstructed plots, the outer error remains
the statistical-plus-conservative-systematic bound. MC errors there retain the
existing independent-cell approximation. The covariance matrices used in the
reconstruction are unchanged; the marginal bands do not imply independent bins
or supply a new goodness-of-fit test.

Every nonzero experimental denominator is retained with its sign. An open
orange symbol marks a measured outer interval containing zero. An exact-zero
denominator is explicitly undefined, never replaced by an epsilon. Missing MC
and absent geometrical support are distinguished in the numerical output.

Focused ratio axes use -1 to 3. Triangles mark any clipped experimental band
or MC interval; a clipped MC central value carries its numerical label and
retains the open-symbol warning where applicable. Full-range companion
figures include all MC and experimental uncertainty extents. Full-range Born
atlases use independent panel scales, stated in the figure. This matters in
v5: deuteron cells 16 and 40 have ratios about 11.22 and 23.93 because their
measured asymmetries are close to zero. No such cell is silently discarded.

## Matching the observable

Direct Born ratios use the existing native-cell LL/UU prediction against the
published Born data. A1 ratios use the existing inverse-D proxy against the
published A1 values, only for the Q2>1 comparison; the proxy's eta*A2 and
averaging limitations remain. No new Q2>4 A1 measurement is inferred.

Reconstructed comparisons divide the GD11-weighted MC cell projection by the
GD11-weighted experimental cell projection from the same `integrated.json`.
They do not divide event-integrated MC asymmetries by reconstructed data.
The full source-summary hash and Born-reference hash must match before rendering.
GD11 extrapolation, split-cell assumptions and systematic-correlation bounds
remain those documented in [the reconstruction method](hermes-apar-integrated-data.md).

## Add figures to an already completed campaign

On Odysseus, load the existing `herwig/pol` environment, enter
`/home/apapaefs/Projects/HerwigPolarizedPheno`, and run:

```bash
python3 scripts/hermes_data_ratio_gallery.py \
  --campaign campaigns/experimental/HERMES_2007_I726689/hermes_pd_born_3m_20260930_v5 \
  --publish
```

This command reads existing numerical products, creates a new checked revision,
adds a bounded section to both gallery indexes and records it in the manifest.
It does not rerender existing figures. Old index text is backed up in the
campaign's `work/ratio-gallery-index-backups/`. Identical complete revisions
are reused; incomplete or corrupted revisions are retained and rebuilt in a
new sibling. Every output is checksum-verified. The ratio subgallery is
self-contained, so the analysis directory remains portable.

`ratios.json` and `ratios.csv` retain the actual data, MC, ratio, separate
uncertainties, support state and denominator flags. PNG/PDF metadata and the
HTML explanation carry the plotting conventions; these uncertainty definitions
should also be stated in a publication caption.

## V5 validation — 7 October 2026

All 402 native tests passed on Odysseus without skips. The completed v5
gallery contains 62 new PDFs and 62 new PNGs. The preservation audit checked
631 prior files: all 629 figure/numerical products are byte-identical, and
the two HTML indexes differ only by the bounded additional section. All
new and existing gallery links resolve; the production physics signature
is unchanged. Representative A1, reconstructed and focused/full-range Born
figures were visually inspected, including the small-denominator cases.
