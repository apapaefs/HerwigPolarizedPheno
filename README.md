# Herwig Polarized Phenomenology

This repository contains the experimental-data phenomenology layer built on
Herwig Polarized: nineteen Rivet analyses, the registry-driven campaign runner,
checksum-pinned reference inputs, postprocessing and plotting code, focused
tests, and the JHEP manuscript *Phenomenological Investigations of Polarized
Collisions in Herwig 7*.

The generator implementation remains an external dependency. This repository
does not vendor Herwig, ThePEG, generated campaigns, installed libraries, or
ordinary YODA/plot products.

## Repository layout

- `analyses/rivet/{dis,pp}/`: eighteen data-linked Rivet analyses plus the
  internal `MC_POLDIJETS` shower-spin measurement and their metadata;
- `cards/`, `config/`, and `data/`: Herwig cards, measurement registries, raw
  provenance inputs, normalized snapshots, and checksums;
- `scripts/run_phenomenology_campaign.py`: primary campaign entry point;
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
- `main.tex`, `references.bib`, and `figures/`: the Overleaf-compatible paper
  bundle, intentionally kept at repository root.

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

The three COMPASS SIDIS reference-data, estimator, and campaign contracts are
documented in
[`docs/compass-sidis-workflow.md`](docs/compass-sidis-workflow.md).
The five first/second-tranche analyses, schema-6 target contract, source
hierarchy and production gates are documented in
[`docs/sidis-tranches-workflow.md`](docs/sidis-tranches-workflow.md).

The synchronized 2026-08-19 corrected-production workflow, including the
STAR 510 GeV generator-cut gate and the explicit no-production preparation
boundary, is documented in
[`campaigns/control/compatibility-corrected-20260819/README.md`](campaigns/control/compatibility-corrected-20260819/README.md).

Build the paper with:

```bash
make paper
```

This requires `latexmk` and a suitable TeX installation.

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

## Provenance and paper mirroring

The phenomenology code/data snapshot was imported from `apapaefs/HerwigPol`
commit `822911203d431daa2a026fc697aabf5eeeb060b2`. The paper history is retained
through audited manuscript commit `f6af80a`, descending linearly from the
original Overleaf import. The repository `main` branch is mirrored to the
legacy paper GitHub repository used by Overleaf; see
[`docs/overleaf-mirroring.md`](docs/overleaf-mirroring.md).
