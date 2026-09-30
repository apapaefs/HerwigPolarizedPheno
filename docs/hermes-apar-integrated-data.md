# Reconstructed HERMES Born longitudinal projections

This calculation constructs proton and deuteron experimental projections of
the published Born `A_parallel` cells using independent unpolarized cross-section
fits. These are derived, model-assisted results, rather than additional HERMES
published measurements. The existing event-level Herwig projections remain
available separately.

## Estimator

For an output region R and published cell i, integrate the unpolarized Born
cross section from the HERMES GD11-P or GD11-D fit over their intersection:

```
U[R,i] = integral_(R intersection cell_i) d sigma_UU^GD11
W[R,i] = U[R,i] / sum_j U[R,j]
A_data[R] = sum_i W[R,i] A_data[i]
A_Herwig[R] = sum_i W[R,i] A_Herwig[i]
```

The same fixed, positive weights act on data and the Herwig published-cell
predictions. They are independent of polarized PDFs and of Herwig's
unpolarized cross-section predictions. This comparison tests the same
compressed-cell estimator on both sides. It is distinct from Herwig's direct
event-level ratio of integrals, whose weights come from its own cross sections.

The inputs are the unchanged 45 Born cells per target and their statistical
covariance from [HERMES 2007](https://www.hepdata.net/record/ins726689).
The independent unpolarized fit is GD11-P/D from
[HERMES 2011](https://arxiv.org/abs/1103.5704). Its parameterization and numerical
source provenance are pinned in `data/experimental/HERMES_2007_I726689/`.
R1998 is used solely to convert the unpolarized structure function into its
lepton scattering cross section. No depolarization factor is applied to
`A_parallel`, and no additional deuteron D-state or published-data correction
is applied to the experimental cells.

## Integration regions and assumptions

Use the physical Born-cell x edges, not the slightly rounded historical
one-dimensional edges. The Q2-integrated projection has 15 x bins over
0.02124 < x < 0.9 and 1 < Q2 < 20 GeV2. It combines whole published cells and
needs no assumption about the asymmetry inside an individual cell.

The x-integrated projection uses Q2 edges 1, 1.5, 2, 3, 4, 6, 8, 12, 20 GeV2.
Those edges split some published cells. Assigning a fraction of a cell's
unpolarized cross section its published asymmetry assumes that asymmetry is
constant inside that source cell. The same assumption acts on the Herwig
cell results. This resolution limitation is not represented by a new
statistical error. The Q2 > 4 x projection also splits cells and has the
same limitation.

The integration uses the actual analysis cuts: a 27.6 GeV lepton beam,
0.10 < y < 0.91, W2 > 3.24 GeV2 and 0.04 < theta < 0.22 radians.
The eight cells below Q2 = 1 GeV2 carry zero weight. No asymmetry is inferred
outside the published cell union.

GD11's stated fit domain is W2 > 4 GeV2. Its use in the accepted region
3.24 < W2 < 4 is an extrapolation, recorded in the derived output. Removing
that part from the weights is an acceptance sensitivity check, rather than
a measurement of a new cut on the unavailable underlying experimental events.
In the completed v4 sample's reconstruction, this extrapolated region supplies
about 45% of the fitted weight in the final x bin. The fit-domain control is
therefore shown separately and should accompany interpretation of the high-x
points.

## Uncertainties

The full statistical covariance propagates as `C_projected = W C_cells W^T`.
Cells shared by multiple projection bins induce correlations between those
bins. The output retains the complete projection operator and covariance;
the plotted point errors alone do not describe their correlations.

Only pointwise experimental systematic uncertainties are public. The default
systematic error is the conservative correlation bound
`sum_i abs(W[R,i]) * s_i`. A diagonal-correlation alternative is reported as
a diagnostic. These are stated assumptions, not a recovered HERMES systematic
covariance. The published systematic column already includes normalization;
no second normalization term is added or subtracted.

Errors are conditional on the central unpolarized fit and the cell-constant
assumption. Weight-model variations are reported as sensitivity shifts, not
as additional one-sigma errors. The public rounded GD11 parameter covariance
is audited separately; an indefinite rounded covariance is not repaired
silently or turned into a variance with an absolute value.
The default controls use ALLM97 proton weights (a stated GD11-ratio hybrid for
the deuteron), R1990 in the unpolarized conversion, W2 > 4 weights, and doubled
quadrature order. Their data and Herwig shifts are recorded separately; the
quadrature comparison is a numerical check.

The Herwig statistical errors use the available cell errors with an explicitly
stated independent-cell approximation. The production summary does not supply
a full cross-cell Monte Carlo covariance. No covariance chi-square or precision
claim follows from this reconstructed comparison.

## Running and outputs

The helper reads an existing production summary and writes a separate output
directory. It launches no events and does not modify the campaign manifest,
generation signature, raw YODA, original reference snapshot or existing plots.

```
python3 scripts/hermes_apar_integrated.py \
  --summary campaigns/experimental/HERMES_2007_I726689/hermes_pd_born_3m_20260930_v4/postprocess/summary.json \
  --output campaigns/experimental/HERMES_2007_I726689/hermes_pd_born_3m_20260930_v4/derived-integrated-data
```

`integrated.json` stores the weights, projected values, statistical covariance
and assumptions. `integrated.csv` provides values and separate uncertainty
components. `index.html` links the proton/deuteron plots in both integration
directions and their Q2 > 4 controls.
