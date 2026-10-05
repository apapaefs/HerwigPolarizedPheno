# HERMES Born longitudinal asymmetry and integrated projections

The active `HERMES_2007_I726689` workflow now compares the direct longitudinal
asymmetry with the unfolded Born proton and deuteron values in Tables XI and
XII of [the HERMES paper](https://arxiv.org/abs/hep-ex/0609039), exposed as
[HEPData Table 7](https://doi.org/10.17182/hepdata.11211.v1/t7) and
[Table 8](https://doi.org/10.17182/hepdata.11211.v1/t8). Each target has 45
published cells in 19 x slices. The measured-asymmetry column is retained in
the raw source for verification; the comparison uses only the Born column.

## Estimator and corrections

For every ordinary normalized bin, first subtract the stored nonnegative
NEGNLO magnitude from POSNLO. Then form

```text
sigma_UU = (PP + PM + MP + MM) / 4
sigma_LL = (PP + MM - PM - MP) / 4
A_parallel = sigma_LL / sigma_UU
```

Proton outputs use the proton component. Deuteron outputs use
`sigma_UU^d=(sigma_UU^p+sigma_UU^n)/2` and
`sigma_LL^d=0.925*(sigma_LL^p+sigma_LL^n)/2`, before division. The D-state
factor therefore appears once in the simulated longitudinal numerator. It
does not modify the published data. Independent order/component variances
and the numerator/denominator covariance enter the ratio uncertainty with
their corresponding coefficient products.

The direct estimator uses no `1/D`, R1990 conversion, or subtraction of
`eta*A2`. The ordinary histograms are filled independently of the validity
of the empirical depolarization model. The existing A1 estimator continues
to use its inverse-D moments. The common card selects QCD shower interactions;
the reference is already unfolded to electromagnetic Born level. No experimental polarization, dilution,
QED/detector unfolding, tensor, or normalization correction is repeated.
The paper's 5.2% proton and 5% deuteron normalization uncertainties are
already included in the published systematic errors.

Direct comparison avoids the analysis-level A_parallel-to-A1 conversion.
It retains the generator's leading-power physics: this change adds no
explicit g2, target-mass, genuine twist-3 or other power-correction model.
The deuteron remains a free-p/n impulse approximation with the fixed .925
spin factor and no binding, Fermi-motion, off-shell, shadowing or tensor
structure-function treatment.

## Physical cells and support

The exact numerical internal cuts come from Appendix C, printed page 97,
of the collaboration-listed
[Ehrenfried analysis thesis](https://www.desy.de/~w3hermes/05-LIB/m.ehrenfried.thesis.ps.gz).
They agree with the final paper's
[original Figure 4 EPS](https://www.desy.de/~w3hermes/TRANS/lara.g1_long-4.eps.gz)
within its coordinate precision. The final paper drops the first thesis x
cell and sets the outer Q2 window to 0.18–20 GeV2. The snapshot records
both these physical boundaries and the rounded x labels used by HEPData.
The legacy midpoint-based Q2 estimates are not used by this comparison.

The existing generator selection, 27.6 GeV e+, `1<=Q2<=20`,
`0.1<y<=0.91`, `W2>3.24`, `0.04<=theta<=0.22`, is retained. All 37
supported cells per target lie wholly above Q2=1. The other eight cells
per target lie below that boundary and are explicitly masked in theory
YODA, summary JSON/CSV and controls, while their published points remain
visible. There is no partially filled cell presented as a full-cell
prediction. The Q2=1 PDF boundary remains exploratory.

Panels `/HERMES_2007_I726689/d07-xNN-y02` and `d08-xNN-y02` show the
proton/deuteron cells against Q2, at each published mean Q2. The reference
markers and both total and optional statistical bars use those means;
the theory steps retain their physical bin edges.

## Paper panels of all Born cells

The campaign `plot` and `full` stages automatically create separate proton
and deuteron figures with five columns and four rows. Nineteen panels show
the fixed-x slices in increasing x; the twentieth tile contains the legend.
Each figure contains all 45 published Born A_parallel cells, including the
eight cells below Q2=1 GeV2 with data only. The horizontal axis is logarithmic
Q2, with the same 0.18–20 GeV2 range in every panel. Vertical limits are
shared within each row and match between the two targets; they include zero
and the full experimental and Monte Carlo error bars.

Black data markers are placed at the published mean Q2. Their inner bars
show the statistical error, and their outer bars show the quadrature sum of
the statistical and published combined systematic errors. The published
systematic error already includes the normalization uncertainty. Orange
steps show the nominal Herwig LL/UU ratio in the physical cell boundaries,
with its statistical uncertainty as a band. Empty Monte Carlo cells retain
their data markers without a theory prediction. No depolarization factor,
additional target correction or unpolarized-fit weighting is applied by
this renderer.

The vector PDF and PNG outputs are named
`AParallel_P_BornCells_5x4` and `AParallel_D_BornCells_5x4`. They are embedded
in both campaign and analysis HTML indexes and saved under
`plots/HERMES_2007_I726689/born-cell-panels/<input-checksum>/`, alongside
`born-cells.json`, `born-cells.csv` and a portable `index.html`. The numerical
snapshot records the original cell values, uncertainties, geometry, theory
support and plotting settings. A complete cache is reused only after all
recorded output checksums pass; changed inputs or incomplete outputs produce
a separate revision. Existing campaign outputs are retained.

These figures use the postprocessed cell results, so a completed v4 campaign
needs only its `plot` stage repeated to add them. No new event generation or
postprocessing is required. The implementation is in
[`hermes_born_cell_plots.py`](../scripts/hermes_born_cell_plots.py).

## Both integrated directions

Each projection forms the ratio after integrating the ordinary helicity
cross sections. It does not average the asymmetries of the source cells.

| Integration | Horizontal axis | Proton path | Deuteron path |
| --- | --- | --- | --- |
| Q2 from 1 to 20 GeV2 | x, existing 15 bins | `AParallel_Q2GT1` | `AParallelD_Q2GT1` |
| Q2 from 4 to 20 GeV2 | x, same bins | `AParallel_Q2GT4` | `AParallelD_Q2GT4` |
| x from .0212 to .9 | Q2: 1, 1.5, 2, 3, 4, 6, 8, 12, 20 GeV2 edges | `AParallelQ2_Q2GT1` | `AParallelDQ2_Q2GT1` |
| x from .0212 to .9 | Q2: 4, 6, 8, 12, 20 GeV2 edges | `AParallelQ2_Q2GT4` | `AParallelDQ2_Q2GT4` |

These projection panels are predictions. The published asymmetries and
statistical covariance do not provide the unpolarized yield/cross-section
weights needed to construct an experimental ratio of integrals. A
statistical average of the published ratios would be a different observable.
Data overlays are therefore supplied for the original cells only.

An additional, explicitly model-assisted reconstruction of both integrated
experimental directions is now available through
[`hermes_apar_integrated.py`](../scripts/hermes_apar_integrated.py). It uses
independent GD11 unpolarized cross-section fits, propagates the full published
statistical covariance, and applies the same weights to the Herwig cell
results. The normal campaign `plot` and `full` stages now add its separate
derived plots to the Rivet gallery automatically, preserving this campaign's
event-level projections. See
[`hermes-apar-integrated-data.md`](hermes-apar-integrated-data.md) for the
within-cell assumptions, systematic bound and fit-domain sensitivity.

## References and reproducibility

`data/experimental/HERMES_2007_I726689/born-apar-reference.json` records all
90 Born values, means, errors, boundaries and support masks. Pinned raw
HEPData Tables 18 and 19 supply each target's full 45x45 statistical
correlation matrix, converted to covariance using the published Born
statistical errors. Systematic and joint p/d covariance are not supplied;
no chi-square is computed by this extension. The combined runtime reference YODA
contains 40 panels and 120 points, including the previous 30 A1 points.
Preparation and plotting regenerate it from the pinned snapshots in the
ignored `campaign/reference/` directory, which Rivet searches first. Generated
reference YODAs are not committed; the historical tracked file is unchanged.

Validate the pinned source reconstruction offline with

```bash
python3 scripts/hermes_born_reference_data.py
```

`--build` regenerates the normalized Born snapshot from the pinned raw
sources. The descriptor pins the supplemental snapshot checksum; campaign
signatures also cover all its raw source bytes. Changing reference data,
cell geometry, analysis or outputs requires a fresh campaign tag and new
`.run` files. Existing v3 YODAs cannot provide the new cell/Q2 histograms.

## Campaign preparation

The new histograms share the same generated events and do not enlarge the
nominal 16-run P/N × four-helicity × two-order matrix. A dry-run of the
matched Odysseus budget succeeds with 1,600 shards and 26.4M events:

```bash
python3 scripts/run_experimental_campaign.py prepare \
  --measurement HERMES_2007_I726689 \
  --tag hermes_pd_born_3m_20260930_v4 \
  --posnlo-events 3000000 --negnlo-events 300000 \
  --shards 100 --jobs 100 --seed-base 660726689 --dry-run
```

The new nominal tag is **prepared** in the separate Odysseus source area
`/home/apapaefs/Projects/Herwig/validation/hermes-born-apar-20260930/pheno`:
all 16 logical `.run` files are nonempty, all 1,600 shards are planned,
and no production YODAs exist. The active `herwig/pol` executable and all
seven recorded runtime artifacts, including installed HwMEDIS, match both
the completed proton production and v3. Only the private Rivet plugin was
built. Earlier v3 manifests, events and source areas remain unchanged.

Fresh `hermes_pd_born_smoke_20260930_v4` completed all 60 jobs, requesting
100 events each across five prediction families. Saved outputs were
postprocessed and rerendered with the final runner. All 287 tests pass
on Odysseus as well as locally. The exact nominal preparation plan also
passes a dry-run through the thermalguard-compatible launcher. No nominal
production event command has been run.

To launch the prepared nominal production, use the matched environment and
options:

```bash
cd /home/apapaefs/Projects/Herwig/validation/hermes-born-apar-20260930/pheno
source /etc/profile.d/modules.sh
module purge
module load herwig/pol
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1

python3 campaigns/launchers/hermes-born-apar/run_validation_campaign.py full \
  --measurement HERMES_2007_I726689 \
  --tag hermes_pd_born_3m_20260930_v4 --profile central \
  --posnlo-events 3000000 --negnlo-events 300000 \
  --shards 100 --jobs 380 --seed-base 660726689 \
  --progress-interval 5 --max-listed 32
```

The nominal measurement signature is
`949bc6907f789aee2e33207f85d83896fa16d9f53216010da366ce6829e66cc1`.
The preparation recorded 100 workers. `--jobs` controls runtime concurrency
and may be changed, including to 380, when launching or resuming the same
tag. The stored preparation configuration remains unchanged; the execution
history records the actual worker count. Event counts, shard count, seeds,
physics inputs and prediction families remain locked. The earlier runner
incorrectly treated a worker-count change as an incompatible manifest.
The earlier [preparation guide](hermes-deuteron-preparation-20260930.md)
retains the historical v3 setup and its normalization audit.

## Validation

An independent comparison with paper Tables XI/XII checked 990 numerical
fields across 90 rows: x limits, x/y/Q2 means, measured/Born asymmetries,
and both errors. All agree exactly. Six HEPData-omitted mean-x entries are
explicitly recovered from the paper. Both statistical correlation matrices
are symmetric, unit-diagonal and positive definite. Numerical thesis cuts
are checked against the published Figure 4 geometry.

Focused fixtures check independent ratio gradients, order signs, target
coefficients, ordinary-only loading, both axes, unsupported-cell masks,
direct-00 controls, LO and uncertainty variations, and published-mean
plot coordinates. The final suite passes 287 tests with native YODA available.
A private Rivet 4.1.2/YODA 2.1.2 plugin passes 916 independent raw-moment,
normalization and error checks on 45 on-shell proton and 45 neutron events:
all 37 supported cells fill once, and the eight unsupported cells remain
empty. Actual-YODA postprocessing passes 4,208 independent checks across
60 controlled streams, all five families and all 46 selectors; the maximum
scaled ratio/error difference is 4.20e-15. Six representative Born and
projection plots were rendered and inspected, including published means,
statistical bars and both horizontal axes. Fixture arithmetic is not a
claim of physical agreement with data.
