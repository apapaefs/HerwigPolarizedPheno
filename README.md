# Herwig Polarized Phenomenology

This repository contains the experimental-data phenomenology layer built on
Herwig Polarized: experimental and diagnostic Rivet analyses,
a reference-only open-charm measurement, the
registry-driven campaign runner,
checksum-pinned reference inputs, postprocessing and plotting code, focused
tests, and reproducible figure-generation tools. The manuscript
[*Phenomenological Investigations of Polarized Collisions in Herwig 7*](https://github.com/apapaefs/Phenomenological-Investigations-of-Polarized-Collisions-in-Herwig-7)
has its own paper-only repository.

The generator implementation remains an external dependency. This repository
does not vendor Herwig, ThePEG, generated campaigns, installed libraries, or
ordinary YODA/plot products.

## Repository layout

- `analyses/rivet/{dis,pp}/`: the twenty-one data-linked Rivet analyses and their
  vendored or deterministically generated reference YODA inputs, plus the
  internal `MC_POLDIJETS` and particle-level `MC_POLJETSHAPES` shower-spin
  measurements;
- `cards/`, `config/`, and `data/`: Herwig cards, measurement registries, raw
  provenance inputs, normalized snapshots, and checksums;
- `scripts/run_phenomenology_campaign.py`: primary campaign entry point;
- `scripts/run_mc_poljetshapes_campaign.py`: pinned combined jet-shape runner;
- `scripts/run_experimental_campaign.py`: backward-compatible fixed-target
  engine used by the primary runner;
- `campaigns/control/compatibility-corrected-20260819/`: tracked, runtime-
  locked controller for the compatibility-corrected production handoff;
- `campaigns/control/sidis-tranches-20260828/`: tracked pilot/central SIDIS
  controller with immutable statistical gating;
- `scripts/tests/`: self-contained campaign, reference, postprocessing, and
  Rivet-source tests;
- `docs/`: workflow documentation, the experimental-compatibility audit, and
  the future EIC polarized-DIS study concept;
- The paper source, bibliography and selected figure assets live in the
  [paper repository](https://github.com/apapaefs/Phenomenological-Investigations-of-Polarized-Collisions-in-Herwig-7).

Generated work is written below `campaigns/` and ignored by Git. The tracked
compatibility controller is the sole exception; its generated runtime reports
remain ignored.

## Quick start

Use a shell with the compatible Herwig Polarized, Rivet, YODA, LHAPDF, and
PDF sets described in [`docs/runtime-compatibility.md`](docs/runtime-compatibility.md).

List the registered measurements:

```bash
python3 scripts/run_phenomenology_campaign.py list
```

Run the focused tests and build the Rivet plugin:

```bash
make test
make rivet
```

Run a small end-to-end campaign, for example:

```bash
python3 scripts/run_phenomenology_campaign.py full \
  --measurement HERMES_2007_I726689 \
  --tag local-smoke \
  --smoke \
  --jobs 4
```

The primary runner also supports COMPASS inclusive DIS, polarized COMPASS and
HERMES identified-hadron SIDIS, unpolarized COMPASS pion/hadron and kaon
multiplicities, STAR weak-boson and jet measurements, and the diagnostic
PHENIX prompt-photon comparison. See
[`docs/phenomenology-campaign.md`](docs/phenomenology-campaign.md)
for profiles, physics families, PDF replicas, scale points, sharding, resume,
and recovery.

The loose 510 GeV `MC_POLDIJETS` spin-on/spin-off control, its cross-section
and asymmetry definitions, and ready-to-run commands are documented in
[`docs/mc-poldijets-workflow.md`](docs/mc-poldijets-workflow.md).
The stable-particle splitting-plane, energy-correlator, quadrupole, and
multijet-rate extension is documented in
[`docs/mc-poljetshapes-workflow.md`](docs/mc-poljetshapes-workflow.md).

The three COMPASS SIDIS reference-data, estimator, and campaign contracts are
documented in
[`docs/compass-sidis-workflow.md`](docs/compass-sidis-workflow.md).
The six-analysis SIDIS tranche/charge-ratio suite, schema-6 target contract,
source hierarchy and production gates are documented in
[`docs/sidis-tranches-workflow.md`](docs/sidis-tranches-workflow.md).
The lower-dimensional COMPASS transverse-momentum result and the COMPASS/HERMES
unpolarized azimuthal diagnostics, including ready-to-run commands, are
documented in
[`docs/sidis-diagnostics-workflow.md`](docs/sidis-diagnostics-workflow.md).

The COMPASS open-charm asymmetries from arXiv:1211.6849v2 are registered as
published reference data: 45 points and nine data-only plots, with simulation
disabled pending a validated massive-charm calculation and analysis response.
See [`docs/compass-open-charm.md`](docs/compass-open-charm.md).

To automatically scan the latest compatible results across all registered
analyses and browse every plot in one searchable offline front-end, run
`python3 scripts/build_results_browser.py --tar` on the campaign host. See
[`docs/results-browser.md`](docs/results-browser.md) for selection rules and options.

To browse completed SIDIS campaigns offline and export a compact tarball of
plots, cuts, data provenance and numerical CSVs without campaign logs or event
files, see [`docs/sidis-results-browser.md`](docs/sidis-results-browser.md).

The synchronized 2026-08-19 corrected-production workflow, including the
STAR 510 GeV generator-cut gate and the explicit no-production preparation
boundary, is documented in
[`campaigns/control/compatibility-corrected-20260819/README.md`](campaigns/control/compatibility-corrected-20260819/README.md).

Build the manuscript with `make paper` in the separate paper repository.
Copy only selected final figures and their provenance there; keep the
analysis, generation and plotting programs here.

## Scientific status

These implementations provide exact, checksum-checked reference-data
transcription and conditional truth/generator-level comparisons. They are not
full reproductions of detector, trigger, background, response, or complete
experimental covariance chains. In particular:

- all fixed-target production made before the compatibility audit must be
  regenerated with the corrected prompt-lepton and bin definitions;
- the HERMES low-`Q^2` projection still needs authoritative per-`x` cell
  boundaries;
- STAR weak-boson and jet results are explicitly generator-level proxies;
- the STAR 510 GeV result needs a newly signed covariance-corrected campaign;
- PHENIX prompt photons remain diagnostic-only because the available hard
  process is incomplete for a publication-level comparison;
- measurements already entering NNPDFpol2.0 are closure tests, not independent
  validations of that fit.

The authoritative measurement-by-measurement assessment is
[`docs/experimental-analysis-compatibility-audit.md`](docs/experimental-analysis-compatibility-audit.md).

## Future EIC study

The proposed next phenomenology project is a study of longitudinal spin
asymmetries in EIC jet production at NLO+parton-shower accuracy. Its aim is to
identify polarized-PDF-sensitive observables that remain stable under QCD
radiation, hadronization, and DIS reconstruction, without performing a PDF
fit, reweighting, or projected PDF-set analysis.

The scoped physics programme, measurement hierarchy, uncertainty budget,
validation requirements, and implementation roadmap are documented in
[`docs/future-eic-study/README.md`](docs/future-eic-study/README.md).

## Repository ownership and provenance

This is the canonical repository for Rivet analyses, generator cards,
reference data, campaign controllers, postprocessing, plotting and tests.
Herwig/ThePEG source remains an external dependency provided by
[herwigdispol](https://github.com/apapaefs/herwigdispol).

The phenomenology code/data snapshot was imported from `apapaefs/HerwigPol`
commit `822911203d431daa2a026fc697aabf5eeeb060b2`. Its inherited paper history
is preserved, but code and paper now advance independently. Do not push code
branches to the paper repository or mirror their `main` refs.

See [the separation audit](docs/repository-separation-audit.md) for the
version choices and [the repository workflow](docs/overleaf-mirroring.md)
for publication and Overleaf instructions.
