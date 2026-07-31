# Figure provenance

The figures under `figures/paper-central-3m-20260723/` were copied without
numerical modification from the `paper_central_3m_20260723` campaign in the
former HerwigPol `DISPOL/campaigns/` run area. Generated campaigns are not
vendored in this repository; the retained figure files and this manifest are
the manuscript provenance record.

The displayed prediction uses the nominal polarized full-spin family, the
central NNPDF4.0 and NNPDFpol2.0 members, and the central hard scale.  The DIS
samples contain 3,000,000 POSNLO and 300,000 NEGNLO events per physical
helicity and target component.  The proton--proton samples contain 3,000,000
LO events per channel and physical helicity.

| Paper directory | Measurement | Primary overlays |
|---|---|---:|
| `dis/` | Fixed-target inclusive DIS, retained on disk but withdrawn pending prompt-lepton regeneration | 0 |
| `hermes-legacy/` | HERMES 2007 Born-level low-\(Q^2\) projection, retained on disk but withdrawn because the \(Q^2\)-cell boundaries are unverified | 0 |
| `hermes-sidis/` | HERMES 2019 identified-hadron SIDIS, retained on disk but withdrawn pending corrected processing and a fresh three-dimensional run | 0 |
| `star-weak/` | STAR 2019 \(W^\pm\) and \(Z/\gamma^*\) | 5 |
| `star-jets-200/` | STAR 2021 inclusive jets and dijets at 200 GeV | 4 |
| `star-jets-510/` | STAR 2022 inclusive jets and dijets at 510 GeV, retained on disk but withdrawn pending a newly signed covariance-corrected campaign or validated migration | 0 |

The pre-audit COMPASS 2010 overlay is also retained under `dis/` but is not
included in the manuscript because its event histogram used incorrect
\(x\)-bin edges.  The manuscript therefore contains 9 primary data-overlay
predictions.
`star-weak/` also carries the five corresponding pull panels retained in the
draft.  Diagnostic closure plots, HERMES azimuthal moments, alternative
comparison families, and plots without experimental reference points are not
included.
