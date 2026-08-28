# SIDIS implementation assessment and roadmap

**Updated:** 2026-08-28
**Repository:** `HerwigPolarizedPheno`
**Authority:** versioned experimental releases and frozen paper-source audits

## Current status

The repository now contains nineteen data-linked Rivet analyses. Nine are
dedicated COMPASS/HERMES SIDIS implementations: the previously available
COMPASS 2009 and two 2017 analyses, plus the five first/second-tranche analyses
and the high-z charge-ratio addition listed below. HERMES 2019 longitudinal
SIDIS remains available separately.

| Tranche | Analysis | Observable | Source authority | Central jobs | Status |
|---|---|---|---|---:|---|
| 1 | `COMPASS_2026_I3096394` | corrected isoscalar charged-hadron, pion and kaon `dM/dz` | HEPData v1 | 4 | implemented; primary modern isoscalar result |
| 1 | `COMPASS_2025_I2840545` | proton charged-hadron, pion and kaon `dM/dz` | HEPData v1 | 2 | implemented |
| 1 | `COMPASS_2010_I862410` | proton identified-hadron `A1(x)` | paper TeX/PDF | 8 | implemented with frozen extraction audit |
| 2 | `HERMES_2013_I1208547` | H/D pion and kaon multiplicities in five 3D binnings | full HERMES archive; HEPData projections as checks | 4 | implemented, 6,592 retained target/species cells |
| 2 | `COMPASS_2018_I1624692` | isoscalar `d2M/(dz dPhT2)` | HEPData v1 | 4 | implemented, 4,664 released cells |
| next | `COMPASS_2020_I1788430` | isoscalar high-z antiproton/proton and K-/K+ multiplicity ratios | paper TeX/PDF | 4 | implemented, 67 published cells |

The complete COMPASS 2018 HEPData v1 submission contains 4,664 numerical
rows, 2,332 per charge, rather than the 4,918 quoted in the earlier planning
assessment. This is a source-level discrepancy, not a removable-cell choice.
All released rows are retained and the missing 254 claimed cells are not
invented.

Each implementation vendors the complete source release, normalized JSON,
checksums, and an audit or inventory. Reference YODA is generated
deterministically from the normalized snapshot by `fetch-data` or `prepare`;
under the current workspace guardrail the new 2020 YODA is not committed.
Multiplicity analyses
store event-aggregated hadron numerators, inclusive-DIS denominators and
same-event covariance proxies. Shards are combined within POSNLO and NEGNLO,
the signed contributions are added at normalized-bin level, target components
are combined, and only then is the ratio divided by its published density
width. Missing cells and non-positive denominators are masked.

For the 2020 charge ratios, the event-aggregated numerator and denominator are
the negative- and positive-hadron yields rather than a hadron and inclusive-DIS
yield. Their same-event covariance is retained, the P/N isoscalar sum is formed
before division, and no density width is applied. Non-positive positive-hadron
denominators are masked.

The descriptor interface is schema version 6. A SIDIS descriptor may request
a nonempty subset of proton/neutron target components and maps each published
target to explicit component weights in `postprocess_config.target_outputs`.
Per-dataset density widths and the longitudinal target scale are data-driven.
Schema-5 descriptors keep their previous behaviour. Unpolarized SIDIS varies
only the unpolarized PDF axis; non-central polarized-PDF selectors are rejected.

## Source hierarchy and uncertainty treatment

- Modern COMPASS multiplicities use the version-pinned HEPData submissions.
  The 2026 corrected release supersedes the 2017 pion/kaon records for the
  primary isoscalar comparison; the 2017 implementations remain reproducible.
- COMPASS 2010 uses the paper TeX/PDF as numerical authority. Its asymmetries,
  per-x statistical correlation blocks, systematic decomposition and published
  deuteron corrections are frozen in a row-by-row audit. The proton prediction
  has longitudinal scale one and a shared 6% multiplicative nuisance.
- HERMES 2013 uses the checksum-pinned 25,937,678-byte VM-subtracted archive
  (`e54ae9e96fb417f43c7845e11319977c21b4c0f7349f00ca987658e00b181e1b`)
  for the five 3D binnings. Dense released statistical covariances are retained
  within each target/species/binning after a frozen normalization of source
  transpose differences and nine one-sided decimal/exponent corruptions. The
  audit records raw/normalized matrix hashes and every replacement; no
  diagonal is changed. Point-to-point systematics are not correlated and
  overlapping binnings are never combined into one fit. The 104 official
  projections are formed from separately integrated generator numerator and
  denominator yields.
- COMPASS 2018 compares primarily to the vector-meson-subtracted values while
  retaining corrections and unsubtracted values as provenance. Its released
  statistical and systematic uncertainties are diagonal.
- COMPASS 2020 has no official HEPData submission as of 2026-08-28. The
  checksum-pinned arXiv:2003.11791 TeX/PDF tables are the numerical authority;
  the arXiv:1802.00584 TeX source is also pinned because the new paper inherits
  its kaon selection. All 67 rows are frozen in an extraction audit. The
  quoted systematic-correlation range 0.7--0.8 is represented by coefficient
  0.75 within each published table: a `sqrt(0.75)` correlated amplitude and a
  `0.5` diagonal amplitude. The overlapping proton x-z and z-momentum tables
  are never combined into one fit.

These are truth/particle-level generator comparisons, not detector-level
reproductions. The isoscalar/deuteron predictions use explicit free-proton and
free-neutron components; no nuclear binding, transport, detector, radiative,
or experimental correction is applied to Herwig.

## Campaign readiness and event estimates

The tracked controller is
`campaigns/control/sidis-tranches-20260828`. It supplies runtime verification,
compact dry runs, immutable preparation, launch/recovery, postprocessing,
plotting, packaging and a statistics assessor. Production remains an explicit
user action.

| Analysis | Pilot POS/NEG per target-helicity | Central POS/NEG floor | Shards per logical job | Central total |
|---|---:|---:|---:|---:|
| `COMPASS_2026_I3096394` | 3M / 0.3M | 30M / 3M | 300 | 66M |
| `COMPASS_2025_I2840545` | 3M / 0.3M | 30M / 3M | 300 | 33M |
| `COMPASS_2010_I862410` | 1M / 0.1M | 3M / 0.3M | 100 | 13.2M |
| `HERMES_2013_I1208547` | 3M / 0.3M | 30M / 3M | 300 | 66M |
| `COMPASS_2018_I1624692` | 10M / 1M | 100M / 10M | 1,000 | 220M |
| `COMPASS_2020_I1788430` | 10M / 1M | 100M / 10M | 1,000 | 220M |

The pilots total 64.9M events. The starting central suite is 618.2M events in
11,800 shards. Extrapolation from completed HERMES production and sparse-bin
smokes suggests about 7 hours of pure generation at 100 Odysseus workers;
12--30 hours is the safer envelope once high-dimensional and rare-high-z Rivet
work and I/O are included. The assessor requires no masked primary bins and
recommends samples
for which at least 90% of bins have `sigma_MC <= 0.5 sigma_exp,total` and every
finite bin has `sigma_MC <= sigma_exp,total`. It preserves the 10:1 POS/NEG
ratio and never recommends below the floors. A bin with nonzero experimental
data but a zero prediction and zero Monte Carlo variance is treated as
unresolved rather than infinitely precise, so it triggers the doubled pilot
instead of an invalid `1/sqrt(N)` extrapolation.

The unpolarized paper profiles contain 103 variations; COMPASS 2010 contains
203. At the central floors the combined paper programme would require about
65.0 billion events, so it must be staged one analysis at a time after the
central gate.

## Next analysis candidates

The strongest recommended unpolarized addition, `COMPASS_2020_I1788430`, is
now implemented. Ratios remove the inclusive-DIS normalization and directly
probe charge and species dependence in fragmentation at the edge of the
current acceptance. `COMPASS_2018_I1652831` remains the logical next
source-consistency cross-check: it supplies the earlier kaon-only result that
the 2020 release extends, but should not supersede the later measurement.

After those ratios, the practical hierarchy is:

1. lower-dimensional COMPASS transverse-momentum multiplicities as shower,
   recoil, intrinsic-transverse-momentum and hadronization diagnostics;
2. COMPASS/HERMES unpolarized azimuthal harmonics as diagnostic-only Fourier
   moments, with no claim of a controlled Cahn/Boer--Mulders/twist-3 model;
3. older longitudinal asymmetries only where complete numerical releases and
   correction conventions can be frozen reproducibly.

Transverse-target observables, Collins/Sivers/pretzelosity measurements,
dihadron transversity, polarized fragmentation or hyperon spin transfer, and
nuclear-medium attenuation/broadening remain outside the present physics
model. They require new distributions, fragmentation dynamics, spin-aware
hadronization, and/or nuclear transport rather than another Rivet histogram.

## Collaborator material

The three collaborator-supplied files in `SIDIS_ByFrank/` are tracked
byte-for-byte unchanged. `provenance-manifest.json` freezes their sizes and
SHA-256 checksums. `.DS_Store` and generated figures are deliberately excluded.

## Validation outcome and remaining production gates

The implementation gate passed locally. All version/checksum and normalized-
source reconstructions passed, including a full eigenvalue validation of the
40 HERMES covariance matrices. The integrated 155-test suite passed with two
unrelated environment-dependent skips. All twenty Rivet sources (the nineteen
data-linked analyses plus the internal `MC_POLDIJETS` control) built in the
aggregate plugin and all twenty metadata lookups succeeded. The phenomenology
paper built with no undefined citation or reference.

Six central `full --smoke` campaigns completed their exact 4/2/8/4/4/4 job
matrices with 100 events per logical job. Every raw job succeeded, all
postprocessors wrote their flattened outputs, and every plot stage completed.
The HERMES smoke alone produced all 1,408 expected reference/prediction panels
as both PDF and PNG. The COMPASS 2020 smoke contained all nine raw estimators
per target/sign job and produced 19 reference/prediction objects and 19 plots
in both formats. Empty low-statistics regions are not interpreted as physics.
The smoke runtime was Herwig devel with Rivet 4.1.2 and installed
`HwMEDIS.7.so` SHA-256
`07e78a9d313f3a66d98befeca6bfd0760499a5ba00bd7ba9bb1306c78ce24e8c`.

Pilot, central and paper production are intentionally not part of this
implementation. At synchronization time the new 800M-event `MC_POLDIJETS`
campaign was actively using the Odysseus primary checkout, so that checkout
was deliberately left at its pinned production commit `bbd4561`. The detached
`HerwigPolarizedPheno-sidis-tranches-20260828` worktree was fast-forwarded to
published `main`; it passed the runtime/source/checksum/registry gate, the full
HERMES covariance validation, and all three campaign dry runs. The active
runtime was not rebuilt and no SIDIS production was prepared or launched.
Local `main`, both GitHub remotes, and the detached Odysseus worktree match;
the primary checkout must be fast-forwarded only after its active Herwig and
campaign-runner processes finish. Its untracked STAR archive remains unchanged
with SHA-256
`9917945085736df4f086e2901fde7dac15d9e65f348741450f0f2d4e04fa67fc`.
