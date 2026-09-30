# HERMES polarized-DIS experimental campaign

This workflow compares HerwigPol with the HERMES 15-bin proton and deuteron
virtual-photon asymmetry measurements in Table 14 of
[HEPData record 11211](https://doi.org/10.17182/hepdata.11211.v1/t14).
It is independent of the general DIS validation campaign. Independent proton
and neutron samples supply both targets, with four separately generated
physical helicities for each component. The proton prediction remains
`d14-x01-y01`; the deuteron prediction is `d14-x01-y02`.

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
- the published `0.0212 <= x <= 0.9` range.

The same quantities are filled for a nested `Q2 > 4 GeV2` selection. The
published `Q2 > 1 GeV2` view reaches the `QMin = 1 GeV` boundary of both
NNPDF40 NLO and NNPDFpol2.0 NLO, so it is exploratory. The `Q2 > 4 GeV2`
view is the conservative HerwigPol validation region.

The R1990 fit supplies `R(x,Q2)` in the longitudinal depolarization factor
`D` defined by `A_parallel = D (A1 + eta A2)`. The prediction sets `A2=0`
and does not model `g2`; [the HERMES paper](https://arxiv.org/pdf/hep-ex/0609039)
reports averaged `g2` contributions of about 0.54% for the proton and 1.9% for
the deuteron in Eq. (20). Raw ordinary, `1/D`-weighted, and
`1/sqrt(D)` covariance-proxy histograms retain the within-sample covariance
needed for the ratio uncertainty.

For each target component and helicity, normalized POSNLO and NEGNLO bins
are added before the physical combinations

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
cd /Users/apapaefs/Projects/HerwigPol/.worktrees/hermes-deuteron-20260930
source /opt/homebrew/opt/modules/init/bash
module purge
module load herwig/pol
```

The primary entry point is `scripts/run_phenomenology_campaign.py`; it
delegates this fixed-target measurement to the experimental runner. The
stage commands below are valid through either entry point. See
[the 30 September preparation note](hermes-deuteron-preparation-20260930.md)
for the prepared nominal production tag and launch command.

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
  --tag hermes_pd_smoke_20260930_v1 --smoke --jobs 4
```

Add `--comparisons` to smoke-test all five prediction families. Smoke mode
uses 100 events for each of the 60 logical jobs, including LO:

```bash
python3 scripts/run_experimental_campaign.py full \
  --measurement HERMES_2007_I726689 \
  --tag hermes_pd_comparison_smoke_20260930_v1 \
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

Campaign configuration is immutable. The proton-only historical tags,
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
  --tag hermes_pd_spin_suite_300k_20260930_v1 \
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
