# HERMES polarized-DIS experimental campaign

This workflow compares HerwigPol with the HERMES 15-bin proton and deuteron
virtual-photon asymmetry measurements in Table 14 of
[HEPData record 11211](https://doi.org/10.17182/hepdata.11211.v1/t14).
It also compares direct longitudinal asymmetries with the 45 unfolded Born
cells per target in Tables 7 and 8, and provides both Q2-integrated x and
x-integrated Q2 predictions. The cell boundaries, support masks, projection
definitions and source validation are documented in
[the Born A_parallel report](hermes-born-apar-20260930.md).
It is independent of the general DIS validation campaign. Independent proton
and neutron samples supply both targets, with four separately generated
physical helicities for each component. The proton prediction remains
`d14-x01-y01`; the deuteron prediction is `d14-x01-y02`.

The current v5 configuration fixes missing physical generation support and
uses a rectangular primary aperture with a shared-event ring control. See
[the 6 October preparation guide](hermes-v5-preparation-20261006.md) for
the exact matched production budget and validation status. V4 remains a
historical comparison, not a corrected prediction.

## Physics definition

The analysis reconstructs the fixed-target laboratory kinematics from the
incoming positron momentum `k`, target-nucleon momentum `P`, and outgoing
positron momentum `k'`:

```text
q = k - k'
Q2 = -q2
x = Q2 / (2 P.q)
y = P.q / P.k
W2 = (P + q)2
theta = angle(k, k')
```

The event selection is

- a 27.6 GeV positron incident on a proton or neutron at rest;
- `1 <= Q2 <= 20 GeV2`;
- `0.1 < y <= 0.91`;
- `W2 > 3.24 GeV2`;
- `0.04 <= theta <= 0.22 rad`;
- primary aperture `abs(theta_x)<0.17` and `0.04<abs(theta_y)<0.14 rad`;
- `0.0212 <= x <= 0.9` for historical one-dimensional projections,
  with the thesis-precision x edges for native Born cells.

The ring control omits only the projected-angle requirements and is filled
on the same events. Rectangle and ring results retain separate acceptance
identities in summaries, compact plots and GD11 reconstruction weights;
their MC errors are correlated. The generator is intentionally looser:
`MinW2=2 GeV2` and `Miny=.095`. Its massless light-cone W2/y variables are
not identical to the physical quantities above. Equal numerical cuts in v4
removed accepted physical events and require new generation to correct.

The same quantities are filled for a nested `Q2 > 4 GeV2` selection. The
published `Q2 > 1 GeV2` view reaches the `QMin = 1 GeV` boundary of both
NNPDF40 NLO and NNPDFpol2.0 NLO, so it is exploratory. The `Q2 > 4 GeV2`
view is the conservative HerwigPol validation region.

For the A1 estimator, the full three-fit R1990 average supplies `R(x,Q2)`
in the longitudinal depolarization factor
`D` defined by `A_parallel = D (A1 + eta A2)`. The prediction neglects
`eta*A2` and does not model `g2`. HERMES's final published `A1` extraction
uses fitted nonzero `g2` through Eqs. (22), (40), and (41); this prediction
therefore retains an approximation to that extraction. The paper's 0.54% and
1.9% average `g2` terms in Eq. (20) are not bounds on the omitted `eta*A2`
term in `A1`, particularly at high `x`. Setting `g2=0` would not set `A2=0`.
Raw ordinary, `1/D`-weighted, and
`1/sqrt(D)` covariance-proxy histograms retain the within-sample covariance
needed for the ratio uncertainty. Direct A_parallel outputs use ordinary
histograms only and do not apply D or an eta*A2 subtraction. This observable
change does not add explicit g2, target-mass or twist-3 physics to Herwig.

The active DIS generator stores both POSNLO and NEGNLO as nonnegative
magnitudes. For each target component and helicity, HERMES applies the
explicit `combination.order_coefficients` map `POSNLO=+1`, `NEGNLO=-1`
to normalized bins before the physical combinations. Independent variances
add with squared coefficients; the ordinary/weighted covariance also retains
both orders because each numerator/denominator sign product is positive.
No sign is inferred from individual bin values. The same convention applies
to direct `00`, all NLO comparison families, and PDF/scale variations; LO
retains coefficient `+1`. The physical combinations are

```text
sigma_UU = (PP + PM + MP + MM) / 4
sigma_LL = (PP + MM - PM - MP) / 4
A_parallel = sigma_LL / sigma_UU
A1 = sigma_LL^(1/D) / sigma_UU
```

The proton view uses the proton component alone. The deuteron view combines
the normalized proton and neutron cross sections before forming either ratio:

```text
sigma_UU^d = (sigma_UU^p + sigma_UU^n) / 2
sigma_LL^d = 0.925 * (sigma_LL^p + sigma_LL^n) / 2
sigma_LL^(1/D),d = 0.925 * (sigma_LL^(1/D),p + sigma_LL^(1/D),n) / 2
```

Thus the deuteron target coefficients are `P=N=0.5` for the denominator and
`P=N=0.4625` for each longitudinal numerator. Covariances use products of
the corresponding numerator and denominator coefficients; independent
proton/neutron variances use squared coefficients. Averaging the two target
asymmetries would give a different estimator and is not used. The proton and
deuteron outputs share proton events, so their Monte Carlo errors are not
independent, although each output retains its own ratio covariance.

The primary acceptance has identity `rectangle_intersect_polar_ring`.
`RingControl_` selectors retain the polar-ring result on the same events.
For each control, the summary includes the ring-minus-primary difference,
its paired statistical error and the primary/control covariance for
A_parallel and, where applicable, A1. The compact acceptance gallery shows
these cell differences with their paired errors; the two standalone error
bars must not be treated as independent.

The spin reduction `0.925 = 1 - 1.5*omegaD`, with `omegaD=0.05`, follows
[Eq. (23) of the inclusive HERMES paper](https://arxiv.org/pdf/hep-ex/0609039).
This per-nucleon impulse approximation omits binding, Fermi motion,
off-shell, shadowing, and tensor structure-function effects. The published
deuteron asymmetry already includes the experimental tensor-asymmetry
correction described in Sec. IV; it is not applied again. The central D-state
factor is fixed; the plotted Monte Carlo band does not include a nuclear-model
uncertainty.

The postprocessor also writes `(PP-MM)/(PP+MM)` and
`(PM-MP)/(PM+MP)` photon-parity residuals for each target view, using the
unpolarized target mixture. A zero denominator masks the bin;
it is never replaced by zero or infinity.

## Reference data

The inspectable snapshot is
`data/experimental/HERMES_2007_I726689/reference.json`. It records the exact
15 bins, their edges, published mean `x` and `Q2`, central values, statistical errors,
and the experimental, parameterization, and evolution systematic components
from Table XXII of the publication for both targets. Central values and
combined systematics are cross-checked against the official HERMES `A1p_a15`
and `A1d_a15` archive members. Both have the same published mean `x` and `Q2`.
The source archive and both member SHA-256 checksums are pinned in the
snapshot and measurement descriptor. The published experimental systematic
column already includes the proton/deuteron normalization uncertainty.

The paper's 15-point A1 values are Q2 averages carrying evolution
uncertainties. Its explicit evolution and covariance-weighted averaging
derivations concern g1; the complete A1 sequence is not separately
reconstructed here. Current predictions integrate selected event cross
sections and do not reproduce the experimental averaging prescription. The
`Q2 > 4 GeV2` curves are internal controls; the published `Q2 > 1 GeV2`
reference is not reused as a measurement with that new cut. The snapshot also
does not provide a full experimental inter-bin covariance matrix.

The equation, correction-factor, and normalized-order audit is documented in
[the specific HERMES validation report](hermes-paper-validation-20260930.md).
The supplementary `born-apar-reference.json` pins all 90 Born points and
two 45x45 statistical covariance matrices. The merged reference YODA has
40 panels and 120 points. Experimental projections are not inferred from
these ratios without the unpolarized yield weights.

Campaign execution is offline. To explicitly download and revalidate the
official source, refresh the cached archive, and regenerate the deterministic
Rivet reference file, run

```bash
python3 scripts/run_experimental_campaign.py fetch-data \
  --measurement HERMES_2007_I726689
```

## Generator configuration

The shared card and the sixteen nominal thin cards are under
`cards/experimental/HERMES_2007_I726689/`. They use

- fixed-target `e+ p` and `e+ n` luminosities at 27.6 GeV;
- `NNPDF40_nlo_pch_as_01180` and `NNPDFpol20_nlo_as_01180`;
- pure photon exchange and native DIS-window POWHEG generation;
- independent `PP`, `PM`, `MP`, and `MM` physical helicities;
- POSNLO and NEGNLO as separate runs;
- the POWHEG real-spin vertex and shower spin correlations;
- the full QCD shower, remnants, hadronization, and decays;
- no MPI and no QED shower radiation.

The ordinary and polarized PDF objects are bound to both `p+` and `n0`.
ThePEG's LHAPDF interface obtains neutron distributions by isospin exchange
`u <-> d` and `ubar <-> dbar` in the proton PDF sets. A separate neutron
LHAPDF set is not required. The incoming `n0` particle must still be present
in the neutron card and accepted by the Rivet analysis.

An opt-in comparison profile adds 44 thin cards below the nested
`comparisons/` directory. `--comparisons` runs five prediction families. The
bold text below is the exact plot and manifest label:

- **`NLO+PS (polarized; full spin)`:** `PP`, `PM`, `MP`, and `MM`, each split into POSNLO and
  NEGNLO, with the exact POWHEG real-emission hard vertex and shower spin
  correlations enabled.
- **`NLO+PS (unpolarized beams)`:** a direct `00` sample split into POSNLO and NEGNLO.
  It uses the same polarized Herwig machinery and PDFs but sets both physical
  beam polarizations to zero. Postprocessing compares this direct
  `sigma_00` with the nominal `(PP+PM+MP+MM)/4`; it does not construct a
  meaningless spin asymmetry from the single `00` sample.
- **`LO+PS (polarized)`:** four physical helicities generated with
  `MEDISNCPol`, the ordinary QCD shower, hadronization, and decays, but no
  POWHEG hard-emission handler. It deliberately keeps the nominal NLO PDF
  inputs so the comparison isolates the hard-process order rather than also
  changing the PDF fit. `--lo-events` controls its per-helicity count and
  defaults to the POSNLO count.
- **`NLO+PS (polarized; Born spin only)`:** four physical helicities,
  each split into POSNLO and NEGNLO, with
  `UsePOWHEGRealSpinVertex No` and shower `SpinCorrelations Yes`. This removes
  the exact `2 -> 3` hard vertex used for real-emission event-record spin
  propagation while retaining the Born spin density, polarized hard-process
  weights, POWHEG emission generation, and shower spin correlations. The
  switch does not alter the real matrix element or POWHEG Sudakov.
- **`NLO+PS (polarized; shower spin off)`:** four physical helicities, each
  split into POSNLO and NEGNLO, with both
  `UsePOWHEGRealSpinVertex No` and shower `SpinCorrelations No`. The polarized
  PDFs, beam helicities, NLO hard-process weights, POWHEG hardest-emission
  generation, QCD shower kinematics, hadronization, and decays remain active;
  only shower-level spin-density propagation is removed.

Each family runs both target components. The expanded profile therefore
contains 60 logical jobs: 16 nominal, 4 direct unpolarized, 8 polarized LO,
16 with only the real-emission spin vertex disabled, and 16 with all shower
spin correlations disabled. Both target views are derived from these same
samples; no additional deuteron event-generation job is needed.

Load the polarized Herwig/Rivet environment before preparing a campaign. The
runner refuses to prepare unless the active Herwig prefix contains `HwMEDIS`
and `FixedTargetLuminosity`, the two PDF sets are installed, and the Rivet
toolchain is usable. The resolved executables, libraries, environment paths,
plugin compiler, configuration, hashes, seeds, and job states are recorded in
the campaign manifest.

## Commands

Run from the source workspace and load the campaign environment first:

```bash
cd /Users/apapaefs/Projects/HerwigPolarizedPheno
source /opt/homebrew/opt/modules/init/bash
module purge
module load herwig/pol
```

The primary entry point is `scripts/run_phenomenology_campaign.py`; it
delegates this fixed-target measurement to the experimental runner. The
stage commands below are valid through either entry point. See
[the v5 preparation note](hermes-v5-preparation-20261006.md)
for the current nominal production tag, preparation command and launch plan.
The September notes retain historical configurations.

Discover registered measurements:

```bash
python3 scripts/run_experimental_campaign.py list
```

Inspect the resolved sixteen-job plan without writing anything:

```bash
python3 scripts/run_experimental_campaign.py prepare \
  --measurement HERMES_2007_I726689 --tag check --dry-run
```

Inspect the expanded 60-job comparison plan:

```bash
python3 scripts/run_experimental_campaign.py prepare \
  --measurement HERMES_2007_I726689 \
  --tag comparison-check \
  --comparisons \
  --dry-run
```

Run a complete 100-event-per-logical-job smoke campaign:

```bash
python3 scripts/run_experimental_campaign.py full \
  --measurement HERMES_2007_I726689 \
  --tag hermes_pd_smoke_20261006_v5_example --smoke --jobs 4
```

Add `--comparisons` to smoke-test all five prediction families. Smoke mode
uses 100 events for each of the 60 logical jobs, including LO:

```bash
python3 scripts/run_experimental_campaign.py full \
  --measurement HERMES_2007_I726689 \
  --tag hermes_pd_comparison_smoke_20261006_v5_example \
  --comparisons \
  --smoke \
  --jobs 4
```

The production defaults are 100,000 POSNLO and 10,000 NEGNLO events per
target component and physical helicity, four concurrent jobs, and one shard
per logical job. The nominal 16-job campaign requests 880,000 events. These
settings can be changed with `--posnlo-events`, `--negnlo-events`, `--jobs`, `--shards`,
and `--seed-base`. Expanded campaigns also accept `--comparisons` and
`--lo-events`. The individual stages are also public commands:
`prepare`, `campaign`, `postprocess`, and `plot`.

Plots show only **`NLO+PS (polarized; full spin)`** by default, including
when the campaign contains the four additional comparison families. The
Monte Carlo statistical uncertainty is drawn as a stepped band spanning each
theory bin, centred on the nominal prediction. Experimental points have one
total uncertainty bar by default. Add `--plot-data-components` only when the
nested statistical-only and total experimental bars are wanted. To overlay
every postprocessed comparison family explicitly, add `--plot-comparisons`
to `plot` (or to `full` together with `--comparisons`):

```bash
python3 scripts/run_experimental_campaign.py plot \
  --measurement HERMES_2007_I726689 \
  --tag comparison-production \
  --plot-comparisons
```

For a label/style-only correction to Rivet `.plot` metadata on an already
complete campaign, `plot --allow-plot-metadata-refresh` permits re-rendering
only if the runner can reconstruct the stored generation signature from the
historical `.plot` file at the campaign's pinned Git commit.  Any concurrent
change to the descriptor, analysis, cards, reference data, or support files is
still refused, and the refresh hashes are written to the manifest.

The generation configuration is immutable. Runtime concurrency (`--jobs`)
may change when launching or resuming an existing tag; its actual value is
recorded in execution history while the preparation configuration is retained.
Event counts, shards, seeds, inputs and prediction families stay locked.
The proton-only historical tags,
including `hermes_prod_300k_20260721`, retain their original scope. Do not
relabel those manifests, append neutron jobs, or attempt a plot-only refresh
to obtain the deuteron prediction. Changed analysis, cards, reference data,
and output target combinations require a fresh tag and regenerated `.run`
files. A nominal-only tag also cannot be expanded by adding `--comparisons`.
The 300,000/30,000 comparison example below is an optional
configuration; the prepared Odysseus run matches the existing HERMES
3,000,000/300,000 nominal production budget in the preparation note.
A 300,000 POSNLO / 30,000
NEGNLO / 300,000 LO comparison production run with ten shards per logical job
contains 10,980,000 events across 600 shard jobs:

```bash
source /opt/homebrew/opt/modules/init/bash
module purge
module load herwig/pol

python3 scripts/run_experimental_campaign.py full \
  --measurement HERMES_2007_I726689 \
  --tag hermes_pd_spin_suite_300k_20261006_v5_example \
  --comparisons \
  --posnlo-events 300000 \
  --negnlo-events 30000 \
  --lo-events 300000 \
  --shards 10 \
  --jobs 4 \
  --seed-base 3726689 \
  --progress-interval 5 \
  --max-listed 32
```

During the Herwig stage, the runner uses the same live tracker as the main
validation campaign. An interactive terminal is refreshed in place with the
completed/running/pending/failed shard totals, logical runs that still have
work, and active-shard event progress, seeds, and runtimes. The same state is
written atomically for external monitoring at

```text
campaigns/experimental/<measurement>/<tag>/monitor/status.txt
campaigns/experimental/<measurement>/<tag>/monitor/status.json
```

The default refresh interval is five seconds and at most twelve logical/active
rows are displayed. Override these without changing the immutable campaign
configuration, for example with `--progress-interval 10 --max-listed 8`. A
negative progress interval suppresses intermediate terminal refreshes while
still writing the initial and final tracker states. Preparation,
postprocessing, plotting, completion, and campaign failure are also reflected
in the monitor phase.

Completed non-empty shard YODA files are recovered from the filesystem and
are not rerun. A failed shard is not silently reused: after inspection, resume
with `--recover-failed` (also accepted as
`--rerun-failed-random-seed`) to assign a fresh unused seed. Configuration
changes require a new tag.

All generated inputs, `.run` files, logs, work products, YODA files,
summaries, and HTML plots are below
`campaigns/experimental/<measurement>/<tag>/`, which is ignored by
Git. Adding another measurement requires a JSON descriptor, Rivet analysis,
reference snapshot, and card family; the runner itself has no HERMES-specific
job-matrix or postprocessing branches.

For expanded campaigns, `postprocess/` contains one prediction YODA per
family, the legacy nominal `summary.json` and `summary.csv`, a
`comparison-summary.json`, and per-family CSV tables. Rivet overlays the four
polarized predictions on `A1` and `A_parallel`; the direct unpolarized sample
appears in `SigmaUU`, accepted-kinematics diagnostics, and the dedicated
`UnpolarizedClosure` observable.

The nominal and comparison summary rows identify four output selections:

| Selection | Target | A1 output |
| --- | --- | --- |
| `Q2GT1` | proton | `/HERMES_2007_I726689/d14-x01-y01` |
| `Q2GT4` | proton | `/HERMES_2007_I726689/A1_Q2GT4` |
| `D_Q2GT1` | deuteron | `/HERMES_2007_I726689/d14-x01-y02` |
| `D_Q2GT4` | deuteron | `/HERMES_2007_I726689/A1d_Q2GT4` |

The descriptor's per-output target maps select the appropriate mixture from
the same raw `SigmaX`, `SigmaOverD_X`, and `CovarianceProxy_X` histograms.
Summary JSON records the resolved maps in `output_target_combinations`.
Accepted-event diagnostics keep their existing proton target convention;
their common names are not deuteron mixtures.
