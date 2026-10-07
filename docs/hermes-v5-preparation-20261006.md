# HERMES v5: physical support, aperture and R1990

This revision addresses the 6 October audit of `HERMES_2007_I726689` against
[HERMES hep-ex/0609039](https://arxiv.org/pdf/hep-ex/0609039). The changes
are prepared on canonical `HerwigPolarizedPheno/main`. They require a fresh
campaign tag and regenerated `.run` files. Existing v4 events, manifests and
plots remain historical products and cannot supply the omitted phase space.

## Physics changes

The physical Rivet cuts remain Q2 in [1,20] GeV2, y in (.10,.91],
W2>3.24 GeV2 and theta in [.04,.22] radians. The primary selection intersects
this polar-angle envelope with the spectrometer aperture
`abs(theta_x)<.17`, `.04<abs(theta_y)<.14` radians. The former polar-angle
ring remains a control filled on the same events, with distinct acceptance
identities in summaries and compact plots. Automatic GD11 reconstruction
uses weights matching the primary rectangle; the integrator separately
exposes a ring option for manual acceptance sensitivity studies. This
is a geometrical aperture model, not a detector simulation. The descriptor
identity is `rectangle_intersect_polar_ring`; control selectors and YODA
objects use the `RingControl_` prefix. The primary paper galleries retain
one set of experimental points. A separate compact acceptance gallery shows
proton/deuteron cell comparisons and ring-minus-primary differences, plus
the x-binned direct asymmetry and A1 proxy for each aperture.

Generation uses `MinW2=2 GeV2` and `Miny=.095`; the upper y and Q2 limits
remain .91 and 20 GeV2. For a=M/(2E), the active fixed-target generator uses
`xi=x/(1+a*x)`, `W2_gen=W2_physical-M^2+a*Q2` and
`y_gen=y_physical*(1+a*x)/(1+a)`. Its old `MinW2=3.24`, `Miny=.1`
therefore cut away physical events accepted by Rivet. The looser generation
window restores room below those physical boundaries; Rivet applies the
measurement cuts. An event at x=.7,Q2=6 GeV2 illustrates the old defect:
physical W2=3.4518 GeV2 passes the paper cut while generator W2=2.6734
GeV2 failed the old generator cut.

The A1 proxy now uses the complete `R1990=(Ra+Rb+Rc)/3` prescription in
[Whitlow's thesis](https://www.slac.stanford.edu/pubs/slacreports/reports11/slac-r-357.pdf),
Sec. 5.3.4, Eqs. (5.31)–(5.35), also described in
[the E143 R1998 paper](https://arxiv.org/pdf/hep-ex/9808028). The old fit B
is itself a published parameterization, but does not equal this average.
This alters inverse-D moments and the alternate-R reconstruction control;
it does not alter the direct ordinary LL/UU estimator. Central GD11 weights
still use R1998 for the unpolarized cross section.

## Factors and interpretation retained

Normalized POSNLO minus the stored NEGNLO magnitude is combined before
helicity and target ratios. UU is `(PP+PM+MP+MM)/4`, LL is
`(PP+MM-PM-MP)/4`. Direct Born A_parallel is LL/UU, with no D or eta*A2
operation. The A1 proxy is inverse-D LL divided by ordinary UU, with D
applied once. Deuteron UU uses .5(P+N), and LL uses .4625(P+N): the .925
D-state factor occurs once in theory. No extra experimental polarization,
dilution, tensor, unfolding or normalization correction acts on data.

The 90 Born reference points, published means and both 45x45 statistical
correlation matrices are retained unchanged. Thirty-seven cells per target
pass the Q2 support requirement; eight lower-Q2 cells remain data-only.
Passing that mask alone does not establish complete generator coverage,
as v4 demonstrated. Q2>4 and the chosen Q2 projection edges are internal
analysis controls, with no corresponding new published measurement.

The generator remains leading power without an explicit finite-Q2 g2,
target-mass, higher-twist or nuclear-dynamics model. The A1 proxy omits
eta*A2 and does not reproduce HERMES's Q2 averaging prescription. The paper
explicitly derives g1 evolution; its A1 averages carry evolution errors
but the full A1 averaging sequence has not been independently reconstructed.
Direct Born data also retain the experiment's unfolding assumptions.

NNPDFpol2.0 includes this HERMES dataset in its fit, so comparison is a
generator consistency check rather than an independent PDF validation.
NNPDF4.0's DIS fit cuts are Q2>=3.5 and W2>=12.5 GeV2; evaluation below
them is beyond fitted DIS kinematics, even where the PDF grid is defined.
The central profile has no PDF or scale band. The fixed omegaD=.05 does not
propagate its quoted .01 uncertainty (a .015 uncertainty on the .925 factor).

GD11 projections remain model-assisted weighted cell combinations with full
experimental statistical covariance. Split cells require a constant
asymmetry assumption. Unknown systematic correlations use a stated bound;
fit-model and acceptance controls are sensitivity shifts, not one-sigma
uncertainties. Shared-event rectangle/ring errors are correlated. The
postprocessor records `ring_minus_primary_a_parallel`, its `_stat` error
and `primary_ring_a_parallel_covariance` in control summary rows, with
corresponding A1 fields where applicable. Difference errors use the paired
covariance; subtracting independent error bars would be incorrect.

## Preparation and production command

Use canonical `main` on Odysseus. The old v4 validation workspace is
preserved unchanged:

```bash
cd /home/apapaefs/Projects/HerwigPolarizedPheno
source /etc/profile.d/modules.sh
module purge
module load herwig/pol
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1
```

The nominal generation matrix remains 16 logical runs: P/N, four helicities
and POSNLO/NEGNLO. Rectangle and ring share those events. Matching v4's
3M/300k per component/helicity and 100 shards gives 1,600 shard jobs and
26.4M requested events. V5 uses seed base `760726689` and tag
`hermes_pd_born_3m_20260930_v5`.

Preparation performs preflight, builds the private Rivet plugin and runs
Herwig read; it does not run production events:

```bash
python3 scripts/run_phenomenology_campaign.py prepare \
  --measurement HERMES_2007_I726689 \
  --tag hermes_pd_born_3m_20260930_v5 --profile central \
  --posnlo-events 3000000 --negnlo-events 300000 \
  --shards 100 --jobs 380 --seed-base 760726689 \
  --progress-interval 5 --max-listed 32
```

Adding `--dry-run` prints the plan but returns before runtime preflight and
does not establish that the executable, libraries, PDFs or plugin work.
The production command below is supplied for a later user launch. Preparing
v5 and running the bounded validation smoke do not authorize this full run:

```bash
python3 scripts/run_phenomenology_campaign.py full \
  --measurement HERMES_2007_I726689 \
  --tag hermes_pd_born_3m_20260930_v5 --profile central \
  --posnlo-events 3000000 --negnlo-events 300000 \
  --shards 100 --jobs 380 --seed-base 760726689 \
  --progress-interval 5 --max-listed 32
```

Retain the original tag's immutable physics, event budget, shards and seeds
on a resume. Only runtime concurrency may change. An incompatible manifest
must not be overwritten or relabelled; use a fresh tag for changed physics.

## Validation required before production

Validation has two separate purposes. Source/reference and arithmetic tests
check all acceptance identities, signed normalization, target factors,
inverse-D moments, support masks, cached-source pins and paired selection
consistency. A real-event smoke must additionally demonstrate that the
installed runtime generates the newly available physical phase space.

In particular, inspect the dedicated low-W2 diagnostic between the physical
3.24 GeV2 threshold and the old effective threshold near 4 GeV2. The
`Accepted_W2Fine_Q2GT1` and `RingControl_Accepted_W2Fine_Q2GT1`
histograms have 0.1 GeV2 bins from 3.24 to 6.04 GeV2. Require actual
accepted entries below the old threshold, including high-x cell 43, and
verify that events failing physical cuts do not leak into primary bins.
Compare cell 43 accepted Q2 distributions and means with v4, retaining MC errors
and the rectangle/ring distinction. Independently verify normalized
POSNLO-NEGNLO sums and the .925 deuteron combination from resulting YODAs.
Render the direct 5x4 figures and both integrated directions to check
acceptance labels, masks and matched GD11 weighting.

The previous GD11 geometry exercise gave cell-43 mean Q2 about 6.31 for
the full ring, 7.10 with the old generator hole and 7.07 in v4 events. These
numbers diagnose the old loss; **6.31 is not a required v5 MC result**.
The PDF, NLO cross section, primary rectangular aperture and sample
fluctuations can change the corrected mean. Do not tune a cut or force a
mean to reproduce that auxiliary unpolarized model. Model-estimated missing
fractions and constant-A1 asymmetry shifts are also not measured Herwig
corrections.

## Bounded validation results and targeted support probe

The native-runtime suite passed all 370 tests. The broad real-event smoke
completed all 160 shards: 100k POSNLO/10k NEGNLO events per target component
and helicity, ten shards and at most 80 concurrent jobs, for 880k requested
events. All seven inspected low-W2 bins from 3.24 to 3.94 GeV2 are
populated for both components and both acceptances. This establishes that the revised generator can produce the previously
excluded physical region.

Its cell-43 support diagnostics are:

| Free-nucleon component | Acceptance | Mean Q2 [GeV2] | Statistical error [GeV2] | Raw entries | Effective entries |
| --- | --- | ---: | ---: | ---: | ---: |
| Proton | Rectangle | 6.31387 | 0.12175 | 45 | 44.77 |
| Proton | Ring | 6.29387 | 0.10276 | 58 | 57.70 |
| Neutron | Rectangle | 6.20900 | 0.20996 | 16 | 15.99 |
| Neutron | Ring | 6.28112 | 0.16113 | 22 | 21.98 |

The moments combine all four helicities into normalized UU, with POSNLO
minus the stored NEGNLO magnitude and event-weighted shard normalization.
Raw entries count events across those streams; effective entries are the
squared signed normalized yield divided by its independent-stream variance.
These are component support checks, before a proton/neutron deuteron mixture.
The broad sample does not meet the predeclared requirement of at least 100
effective entries in cell 43 for each component and each acceptance. Its
means are therefore sparse support diagnostics, not completion of the
precision check. The checker also requires each ring mean to be below
6.8 GeV2, an analyst-chosen diagnostic threshold for moving away from the
old cutoff-biased result, not a HERMES measurement or exact model prediction.

The complete postprocessing and gallery smoke passed, producing 154 PNGs.
The proton 5x4 Born panel and compact acceptance comparison were inspected
visually; labels, masks, logarithmic point placement and layout were intact.
Both reconstructed projection directions and paired acceptance differences
were generated automatically.

A separate targeted probe has been prepared with tag
`hermes_v5_highx_20k_20261006`, seed base `780726689`, 20k POSNLO/2k
NEGNLO events per component and helicity, four shards and at most 32 jobs:
176k requested events in 64 shards. The probe copies the 17 run/common cards
under `campaigns/validation-source/hermes-v5-highx-20261006/cards` and
sets only its generator `X2Min=.5`. The saved cloned descriptor
`measurement.json` and `support_probe` metadata record this restricted
validation setup. Canonical production cards remain unchanged.

Cell 43 has minimum generator light-cone fraction about .657, safely above
the probe's .5 cut. The probe thus increases the event density at high x
without removing this cell's physical acceptance. Its purpose is to repeat
the low-W2 and cell-43 moment checks at adequate statistics. It is not a
full-range HERMES prediction: do not use its missing low-x cells in the paper
gallery, merge it with the broad sample, or substitute it for production.
Its mean is still assessed with MC errors, without tuning to the auxiliary
GD11 number.

Nominal v5 preparation completed on Odysseus: the private plugin and all
16 non-empty `.run` files exist; the manifest has 1,600 planned jobs and
26.4M requested events, with zero production YODA outputs. The active
executable is `/home/apapaefs/Projects/Herwig/Herwig-Pol/bin/Herwig`,
with `HwMEDIS.so.7.1.0` from that prefix. Physics signature:
`0737c5a6c3e2df41554bfb4eb6a08e827f153ca98249d49e4dc37ec775739fd3`.
The production plugin SHA256 is
`c3c7a1a9362028d8002dff478515f1dbf44104f924b537711a4be935210d43d3`.

The targeted probe is still running, with slow NEGNLO generation. Its
precision gate is **pending**, not passed. The generated validation pipeline
will write `support-validation.json` inside its campaign directory when
all shards finish. The report includes raw-file hashes, normalized means,
MC errors, effective counts and explicit failures. The standalone check is:

```bash
python3 scripts/check_hermes_generation_support.py \
  campaigns/experimental/HERMES_2007_I726689/hermes_v5_highx_20k_20261006 \
  --output campaigns/experimental/HERMES_2007_I726689/hermes_v5_highx_20k_20261006/support-validation.json
```

Check that report before launching the nominal production. A failed or
incomplete check requires investigation; do not lower its thresholds to
obtain a pass. No full v5 production has been launched.
