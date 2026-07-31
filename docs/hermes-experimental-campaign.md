# HERMES polarized-DIS experimental campaign

This workflow compares HerwigPol with the HERMES 15-bin proton
virtual-photon asymmetry measurement in Table 14 of
[HEPData record 11211](https://doi.org/10.17182/hepdata.11211.v1/t14).
It is independent of the general DIS validation campaign and uses four
separately generated physical-helicity samples.

## Physics definition

The analysis reconstructs the fixed-target laboratory kinematics from the
incoming positron momentum `k`, target-proton momentum `P`, and outgoing
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

- a 27.6 GeV positron incident on a proton at rest;
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
and does not model `g2`; the HERMES paper reports that the `g2` term averaged
over its proton bins is about 0.54%. Raw ordinary, `1/D`-weighted, and
`1/sqrt(D)` covariance-proxy histograms retain the within-sample covariance
needed for the ratio uncertainty.

For each helicity, normalized POSNLO and NEGNLO bins are added before the
physical combinations

```text
sigma_UU = (PP + PM + MP + MM) / 4
sigma_LL = (PP + MM - PM - MP) / 4
A_parallel = sigma_LL / sigma_UU
A1 = sigma_LL^(1/D) / sigma_UU
```

The postprocessor also writes `(PP-MM)/(PP+MM)` and
`(PM-MP)/(PM+MP)` photon-parity residuals. A zero denominator masks the bin;
it is never replaced by zero or infinity.

## Reference data

The inspectable snapshot is
`data/experimental/HERMES_2007_I726689/reference.json`. It records the exact
15 bin edges, published mean `x` and `Q2`, central values, statistical errors,
and the experimental, parameterization, and evolution systematic components
from Table XXII of the publication. Its central values and combined
systematics are cross-checked against the official HERMES `A1p_a15` archive
member. The source archive and member SHA-256 checksums are in both the
snapshot and measurement descriptor.

Campaign execution is offline. To explicitly download and revalidate the
official source, refresh the cached archive, and regenerate the deterministic
Rivet reference file, run

```bash
python3 scripts/run_experimental_campaign.py fetch-data \
  --measurement HERMES_2007_I726689
```

## Generator configuration

The shared card and the eight nominal thin cards are under
`cards/experimental/HERMES_2007_I726689/`. They use

- fixed-target `e+ p` luminosity at 27.6 GeV;
- `NNPDF40_nlo_pch_as_01180` and `NNPDFpol20_nlo_as_01180`;
- pure photon exchange and native DIS-window POWHEG generation;
- independent `PP`, `PM`, `MP`, and `MM` physical helicities;
- POSNLO and NEGNLO as separate runs;
- the POWHEG real-spin vertex and shower spin correlations;
- the full QCD shower, remnants, hadronization, and decays;
- no MPI and no QED shower radiation.

An opt-in comparison profile adds 22 thin cards below the nested
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

The expanded profile therefore contains 30 logical jobs: 8 nominal, 2 direct
unpolarized, 4 polarized LO, 8 with only the real-emission spin vertex
disabled, and 8 with all shower spin correlations disabled.

Load the polarized Herwig/Rivet environment before preparing a campaign. The
runner refuses to prepare unless the active Herwig prefix contains `HwMEDIS`
and `FixedTargetLuminosity`, the two PDF sets are installed, and the Rivet
toolchain is usable. The resolved executables, libraries, environment paths,
plugin compiler, configuration, hashes, seeds, and job states are recorded in
the campaign manifest.

## Commands

Discover registered measurements:

```bash
python3 scripts/run_experimental_campaign.py list
```

Inspect the resolved eight-job plan without writing anything:

```bash
python3 scripts/run_experimental_campaign.py prepare \
  --measurement HERMES_2007_I726689 --tag check --dry-run
```

Inspect the expanded 30-job comparison plan:

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
  --measurement HERMES_2007_I726689 --tag smoke01 --smoke --jobs 4
```

Add `--comparisons` to smoke-test all five prediction families. Smoke mode
uses 100 events for each of the 30 logical jobs, including LO:

```bash
python3 scripts/run_experimental_campaign.py full \
  --measurement HERMES_2007_I726689 \
  --tag comparison-smoke01 \
  --comparisons \
  --smoke \
  --jobs 4
```

The production defaults are 100,000 POSNLO and 10,000 NEGNLO events per
physical helicity, four concurrent jobs, and one shard per logical job. They
can be changed with `--posnlo-events`, `--negnlo-events`, `--jobs`, `--shards`,
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

Campaign configuration is immutable. In particular, do not add
`--comparisons` while resuming the existing nominal-only
`hermes_prod_300k_20260721` tag; use a fresh tag. A 300,000 POSNLO / 30,000
NEGNLO / 300,000 LO comparison production run with ten shards per logical job
contains 5,490,000 events across 300 shard jobs:

```bash
source /opt/homebrew/opt/modules/init/bash
module purge
module load herwig/pol

python3 scripts/run_experimental_campaign.py full \
  --measurement HERMES_2007_I726689 \
  --tag hermes_prod_spin_suite_300k_20260721 \
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
