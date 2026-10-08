# COMPASS open-charm asymmetries

`COMPASS_2013_I1204782` registers the measured photon–nucleon asymmetries
from [arXiv:1211.6849v2](https://arxiv.org/abs/1211.6849v2), published in
[Phys. Rev. D 87 (2013) 052018](https://doi.org/10.1103/PhysRevD.87.052018).
The registration is **reference data only**. It supplies 45 published points,
nine PNG/PDF panels and the measurement contract. It does not yet supply a
Rivet event implementation, runnable Herwig cards or a simulated comparison.

## Browse and verify

Run from the repository root:

```bash
python3 scripts/build_results_browser.py
```

Open the new index printed by the command, search for `1204782`, and select
**View 9 plots**. The card and every figure say **Reference data only**.
**Analysis info** links to the paper, the complete CSV and JSON, and describes
the cuts, estimator and remaining simulation requirements. Building the data
panels needs matplotlib; viewing them needs only a browser.

To create a small browser containing only this measurement:

```bash
python3 scripts/build_results_browser.py --measurement COMPASS_2013_I1204782
```

Validate the pinned data, including an independent extraction of all 45 rows
from the PDF (the latter needs Poppler's `pdftotext`):

```bash
python3 scripts/compass_open_charm_reference.py check --verify-pdf-text
```

`rebuild` deterministically regenerates `reference.json` and the source
manifest from the pinned CSV and definition. It refuses altered raw inputs.
The main runner's `fetch-data --measurement COMPASS_2013_I1204782` validates
the vendored inputs offline. Campaign preparation, generation, postprocessing
and prediction plotting fail before dispatch for this reference-only entry.

## Numerical contract

The source is the exact **v2 PDF**, checksum-pinned alongside a manually
transcribed CSV. The CSV was checked visually against the tables and then
compared row-by-row with independently extracted PDF text. There are three
tables, each containing five transverse-momentum bins and three laboratory
energy bins:

| Table | Published sample group | Points |
|---|---|---:|
| 5 | Untagged D0→Kπ, D* tagged Kπ and tagged Ksubπ, combined | 15 |
| 6 | D* tagged Kππ0, with the π0 unobserved | 15 |
| 7 | D* tagged Kπππ | 15 |

Charge conjugates are included. The five pT intervals are 0–0.3, 0.3–0.7,
0.7–1, 1–1.5 and >1.5 GeV/c; the energy intervals are 0–30, 30–50 and
>50 GeV. `null` upper edges in JSON and empty upper-edge fields in CSV mean
open-ended intervals. The plots use a categorical pT axis, avoiding an
invented finite upper edge. All tabulated means of y, Q², pT, energy and D
are retained; they are weighted by the experiment's squared signal weight.

Statistical and systematic uncertainties remain separate. The pale outer
bars show their quadrature sum for display; this does not construct a
covariance matrix. Estimates outside [-1, 1], including Table 7's
7.03 ± 4.74 ± 0.71, are retained without clipping. No fit or chi-squared is
computed. The paper's LO/NLO fitted gluon-polarisation values and calculated
analysing powers are not treated as additional independent measurements.

## Boundary of the simulation comparison

The data use 160 GeV/c positive muons on polarised proton and deuteron targets
over 2002–2007, with results combined in the published tables. The dominant
mechanism is photon–gluon fusion into massive charm. The paper reports mean
Q² approximately 0.6 GeV²; the quoted observed Q² and y ranges must not be
turned into extra fiducial cuts. Its hard scale is based on the charm
transverse mass, rather than Q² alone.

There is no validated backend for this measurement configured in this
repository. Reusing the ordinary massless DIS cards, or selecting charm
hadrons appearing in their showers, would not establish the required
massive-charm muoproduction prediction. No installed generator was changed.

Before enabling simulations, implement and validate:

- Massive charm production with the appropriate lepton/photon polarisation
  transfer and helicity-dependent cross sections throughout the low-Q²
  region; validate each claimed perturbative contribution separately.
- A justified proton/deuteron comparison. The published tables do not give
  target-specific mixture fractions; the SIDIS deuteron factor cannot simply
  be inherited.
- Reconstructed D0 candidates and D* tags, the unobserved π0 treatment,
  Table 4's channel-specific cuts, particle identification and candidate
  priorities. A true D0 four-vector is not automatically equivalent to the
  partially reconstructed Kπ system.
- A documented treatment of the acceptance, signal purity and weights.
  Equation (2) gives a finite-muon-mass depolarisation factor. The measured
  asymmetry is A_gammaN = A_muN / D, extracted with
  w_S = P_mu f D s/(s+b), and is not the existing inclusive/SIDIS A1 estimator.
  The table means do not reconstruct the event-level estimator or detector
  response.
- A Rivet implementation, independent normalization/helicity checks and a
  low-statistics end-to-end validation with the confirmed active runtime.

The separate `config/reference/` registry makes this boundary explicit.
Promoting this entry to a runnable analysis requires replacing that contract
with the validated implementation, rather than toggling an enable flag.
