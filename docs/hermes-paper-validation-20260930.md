# HERMES inclusive analysis audit, 30 September 2026

This audit checks `HERMES_2007_I726689` proton/deuteron Table 14 predictions
against [HERMES hep-ex/0609039v2](https://arxiv.org/pdf/hep-ex/0609039),
the official HERMES archive and the active Odysseus generator. The starting
canonical GitHub commit is `c4e9de1f311937db0f4573accfbba8f28a30e3bf`.
The audit distinguishes correction-factor closure from a complete
reproduction of the experiment's finite-Q2 extraction and averaging.

## 6 October follow-up: generation support and acceptance

The factor and reference checks below remain useful, but they did not validate
physical phase-space coverage. The later audit found a generation-window
defect in v3/v4: equal numerical generator and Rivet cuts do not select equal
kinematics for a massive fixed target. Some of the 37 cells above Q2=1 were
only partly generated. The previous conclusion of full support is withdrawn.

For a target at rest, put a=M/(2E). The active generator's light-cone fraction
and massless DIS window obey

```text
xi = x / (1 + a*x)
W2_generator = W2_physical - M^2 + a*Q2
y_generator = y_physical * (1 + a*x) / (1 + a)
```

Thus `MinW2=3.24` removed physical events accepted by Rivet near its W2 cut,
and `Miny=.1` removed a narrow part of the physical lower-y boundary. A v4
example, x=.7 and Q2=6 GeV2, has W2_physical=3.4518 GeV2 but
W2_generator=2.6734 GeV2. Preserving the lepton momentum through the shower
does not cure this difference in variable definitions. V5 uses generation
`MinW2=2 GeV2`, `Miny=.095`, retaining the exact physical cuts in Rivet.
Missing events cannot be restored by replotting or reweighting old YODAs.

V5 also makes the rectangular spectrometer aperture, intersected with the
published polar-angle envelope, the primary selection. The polar-angle ring
is retained as a control on the same events. The former fit-B-only R choice
is replaced by the full three-fit R1990 average for the A1 proxy. See the
[v5 preparation and validation guide](hermes-v5-preparation-20261006.md).
Historical runtime counts and hashes below describe the September audit;
they are not validation evidence for the revised selection.

## Verdict of the September factor audit

The lepton depolarization factor and deuteron spin factor are each applied
once. Published reference values are not corrected a second time. The
reference data and the numerical scalar Rivet cuts match the selected paper
entries. This does not establish generator coverage or the detector aperture.

A separate NLO normalization error was found: the active generator emits
NEGNLO as a positive magnitude, while the old experimental loader added it
to POSNLO. HERMES now explicitly subtracts normalized NEGNLO central values
and adds independent variances. This correction applies to every HERMES
physical helicity, P/N component, NLO comparison family and variation.
Fresh immutable source signatures/tags are required; the v1/v2 preparations
and their smoke curves are historical and must not be interpreted as corrected
NLO predictions. The earlier smoke validation established execution, not the
missing order sign.

The corrected prediction remains an approximation to the published A1:
eta*A2 and a consistent finite-Q2 g2 treatment are absent, and its event
averaging does not reproduce the published Q2-averaged A1 prescription.
This audit does not establish quantitative agreement with the measured curves.

## Factors from event generation to the plotted ratio

| Factor or operation | Location and result |
| --- | --- |
| Lepton/target polarization | Cards generate unit physical helicities. The PDF polarization enters as Pz*Deltaq/q with ordinary PDFs in the flux. No experimental beam/target polarization or target-gas dilution factor is applied to the theory. |
| Lepton angular response | The photon Born factor is 1+ell^2+2*Pl*Pq*ell, ell=2/y-1. Its raw longitudinal asymmetry contains D_LO=2*ell/(1+ell^2), rather than being a photon A1 already. |
| R1990 and finite-gamma depolarization | Rivet computes D once from reconstructed x,y,Q2 and fills ordinary w and weighted w/D. The postprocessor divides the weighted LL by ordinary UU and performs no further D division. |
| Covariance | Rivet fills w/sqrt(D); its squared-weight moment is sum(w^2/D). Ratio covariance uses products of numerator/denominator coefficients. This proxy is statistical and adds no central-value correction. |
| NLO order sign | The installed DISBase clips POS=max(raw,0), NEG=max(-raw,0). The explicit HERMES map is POSNLO=+1, NEGNLO=-1 for stored normalized magnitudes. Signs square in independent variances and ordinary/weighted covariance. |
| Deuteron wavefunction | Per-nucleon UU coefficients are P=N=.5. LL and weighted LL coefficients are P=N=.4625=.925/2, so Eq. (23)'s 1-1.5*omegaD is applied once. Metadata's d_state_factor=.925 is descriptive and is not another multiplication. |
| Experimental corrections | References already include beam/target normalization, background subtraction, QED/detector unfolding and the deuteron tensor correction. The reader/writer copies published A1 directly. No data-side D, .925 or extra tensor correction is applied. |
| Experimental normalization uncertainty | Table XXII already includes 5.2% proton/5% deuteron normalization in the experimental systematic component. It is not added again. |

The signed-order map multiplies stored bins without inferring their sign or
taking absolute values. A generator that stores an already-signed negative
stream would require its own explicit convention. The confirmed current
Odysseus installation stores positive magnitudes.

## Equation and selection checks

Printed pp. 5-6, Eqs. (19)-(22) define sigma_LL as antiparallel-minus-parallel over 2
and Aparallel=sigma_LL/sigma_UU. Equal incoming helicity labels in the collision
frame correspond to opposite physical spin directions, so the four-helicity
combination (PP+MM-PM-MP)/4 has the paper's sign. No additional positron sign
is required for photon exchange.

With gamma^2=4*M^2*x^2/Q2, define

```text
epsilon = (1-y-gamma^2*y^2/4)/(1-y+y^2/2+gamma^2*y^2/4)
D = [1-(1-y)*epsilon]/[1+epsilon*R]
eta = epsilon*gamma*y/[1-(1-y)*epsilon]
A1 = (g1-gamma^2*g2)/F1
A2 = gamma*(g1+g2)/F1
Aparallel = D*(A1+eta*A2)
```

This epsilon form reproduces the implemented D algebra independently. Combining
the paper's Eqs. (20), (37), (40)-(41) closes the last identity to floating-point precision.
HERMES explicitly uses R1990 for its A1 extraction (Sec. VI B, printed p. 16).
The R parameterization is an empirical analysis conversion; it is not a
recalculation of generator-native FL/F1. The September implementation used
Whitlow fit B alone. Fit B is a published form, but the thesis definition of
R1990 is `(Ra+Rb+Rc)/3`; v5 uses that full average. This change affects the
inverse-D A1 estimator, not direct Born Aparallel. Mean-point diagnostics
support the revised convention but do not establish the exact experimental
implementation at every point.

The 27.6 GeV positron beam, Q2 in [1,20] GeV2, y in (.1,.91], W2>3.24GeV2,
theta in [.04,.22] rad and x in [.0212,.9] match the selected Q2>1 publication
region and Fig. 4's scalar envelope. The body text's angular unit typo does not
replace the radian values labelled in Fig. 4. The implementation selects the
highest-energy prompt outgoing positron, consistent with the paper's leading
lepton selection for this QED-radiation-free generator setup. The September
implementation modeled the scalar acceptance envelope. V5 additionally
intersects it with the rectangular aperture, while still omitting detector
response. Cards use
QCD shower interactions and apply no detector/QED radiative correction to the
already Born-unfolded data. The Q2>4 outputs are internal controls with no
separate published reference.

## Independent numerical checks

All 30 P/D reference rows exactly match the pinned official A1p_a15/A1d_a15
members for mean x,Q2,central value,statistical error and combined systematic.
All 240 independently checked Table XXII numerical fields, including its three
systematic components, match the snapshot. All 16 selected x-bin edges match
the published ranges. Separate-systematic quadrature versus the rounded
combined archive values differs by at most 4.58e-5.

An independent rational event fixture varies D within each bin. For P, its
(UU,A1,D) strata are (8,.3,.25) and (2,.6,.75); for N they are (3,-.2,.4) and
(17,.1,.8). Physical helicity event weights are UU*(1+h*D*A1). The actual
extracted C++ fill body and Python postprocessor produce:

| Target | UU | LL | Aparallel | Weighted-LL/UU |
| --- | ---: | ---: | ---: | ---: |
| Proton |10 |1.5 |.15 |.36 |
| Deuteron |15 |1.21175 |.0807833333 |.1449166667 |

Independent Poisson event-count derivatives also reproduce both ratio
variances and ordinary/weighted covariance. The deuteron wrong alternatives
are .3646041667 for a second D division, .1324316940 for division by an
UU-weighted mean D, and .1340479167 for a second .925 multiplication. All
are distinguished by the fixture.

Permanent focused tests are in `scripts/tests/test_hermes_depolarization_audit.py`.
They include compiled extraction of the actual R1990/D/fill functions when a
compiler is available, independent equations and event-count errors, all 30
actual reference values and optional real-YODA P/D verification. Separate
order-magnitude tests check central-value subtraction and unchanged variance
addition, including direct00 and comparison outputs.

## Limits of paper-level validation

HERMES uses a fitted nonzero g2 to obtain final A1 (Sec. V A, printed p. 14;
Sec. VI B, Eqs. (40)-(41), printed p. 16). Neglect of A2 in the Appendix A radiative
kernels is a different step. This implementation estimates an event-weighted
A1+eta*A2 and interprets it as A1 by neglecting eta*A2. Setting g2=0 does
not set A2=0. At the highest published mean (x=.7248, Q2=12.21 GeV2), the
kinematic eta*gamma is .11844; in the illustrative g2=0 case, this would
shift Aparallel/D relative to A1 by 11.84%. This is not an estimate from the
paper's fitted g2, and the paper's .54%/1.9% g2-term averages in Eq. (20) do not
bound this A1 conversion error.

Table XXII's 15 A1 points are Q2-averaged results carrying evolution
uncertainties. The paper explicitly describes g1 evolution in Eq. (38) and
covariance-weighted averaging in Appendix B; it does not separately specify
the complete A1 averaging/evolution sequence. Our accepted-event
cross-section-weighted proxy does not reproduce that experimental
prescription. Correlated statistical averaging is not simply weighting by
the inverse diagonal variances. A faithful precision reproduction needs an
explicit, consistent finite-Q2 g1/g2/unpolarized baseline and matching bin
prescription. No paper-model numerical g2 correction has been silently
added, and no target-mass/higher-twist treatment is inferred from a massive
nucleon in the event record.

The deuteron model is free (p+n)/2 with fixed omegaD=.05, without binding,
Fermi motion, off-shell effects, shadowing or tensor structure functions. The
Monte Carlo error band does not include these nuclear uncertainties or a
cross-target covariance matrix. The original A1 snapshot does not provide
its full experimental inter-bin covariance; the later Born snapshot does
retain both published 45x45 statistical matrices. Q2>1 reaches the selected
PDF boundary; Q2>4 is an analyst-selected internal control, not a published
measurement or a guarantee that power corrections are negligible.

Other descriptors using the historical inclusive order sum are unchanged by
this scoped HERMES repair and require a separate signed-stream audit before
using their old curves as quantitative NLO predictions.

## Historical September runtime validation and replacement preparation

The active Odysseus `herwig/pol` environment passed all **247 regression tests**,
including the independent factor, reference and order-magnitude tests. Fresh
`hermes_pd_smoke_20260930_v3` completed **60/60 jobs** (100 requested events
per job), then postprocessed and plotted all five prediction families. All
four NLO families record `POSNLO=+1, NEGNLO=-1`; LO records `LO=+1`. Summary
numbers are finite or explicitly masked. Both primary P/D plots were inspected.
An independent reconstruction from the 16 normalized nominal input YODAs
reproduces all 60 P/D summary rows: 416 central-value and statistical-error
checks agree within a maximum scaled difference of 4.59e-15. This calculation
imports no production postprocessor. These small samples establish execution
and arithmetic, not data agreement.

Replacement nominal production `hermes_pd_3m_20260930_v3` is **prepared** with
16 nonempty logical `.run` files and 1,600 planned shards. It matches existing
HERMES budgets: 3M POSNLO/300k NEGNLO per target and helicity, 100 shards each,
26.4M total requested events and 100 workers. Seed base is `560726689`.
The thermalguard-compatible launcher passes its exact preparation dry-run.
No production events have been launched. Older v1/v2 products are retained.
See [the preparation guide](hermes-deuteron-preparation-20260930.md) for the
isolated Odysseus path and complete launch command.

Both fresh tags record measurement signature
`0556fc25305904a1391d535be6b6d0b740dc1d968962e93a9005067f5f3e92dd`.
The actual installed `HwMEDIS.so.7.1.0` SHA256 remains
`33d696b1e16b4d3d3f14d7840040de7d009d4c2c5d8a863daf390e14114e3eda`,
matching the previous production runtime. The private Rivet plugin is rebuilt;
no installed generator or historical manifest is modified by the sign repair.
