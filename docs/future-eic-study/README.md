# Future EIC polarized-DIS study

## Status and purpose

This document is a concept note for a possible extension of the Herwig
Polarized phenomenology programme to the Electron-Ion Collider (EIC). It is
not yet an implemented campaign, a detector study, or a polarized-PDF
projection.

The proposed first study is

> **Longitudinal spin asymmetries in EIC jet production at NLO+parton-shower
> accuracy.**

Its central question is deliberately narrower than a global-fit forecast:

> Which polarized-DIS observables at the EIC are simultaneously sensitive to
> proton helicity structure, perturbatively stable, robust under showering and
> hadronization, and realistically reconstructible?

The result would be an observable and uncertainty study. Existing polarized
PDF sets and replicas would be compared as alternative theory inputs, but no
set would be fitted, reweighted, profiled, or projected.

## Why this is a useful next step

Inclusive polarized DIS is the clean baseline for quark-helicity information,
while jet and dijet production adds direct sensitivity to QCD radiation and,
through boson--gluon fusion, to the polarized gluon distribution. The EIC
offers a broad lever arm in Bjorken `x` and `Q^2`, variable beam energies, and
high polarization. It therefore provides an opportunity to ask not only where
an asymmetry is large, but where it survives the full chain from a hard NLO
calculation to experimentally defined objects.

This question is particularly appropriate for Herwig Polarized. The study can
test the interplay of polarized hard processes, POWHEG matching, the
spin-aware initial-state shower, hadronization, and DIS reconstruction. It
also extends the existing phenomenology without pretending that a generator-
level calculation replaces an experimental response model or a global PDF
analysis.

## Recommended scope

### Collision process and accuracy

The baseline process is longitudinally polarized neutral-current DIS,

```text
e- p -> e- + X,
```

including the available photon and `Z` contributions. The primary predictions
should use NLO QCD matched to the parton shower with POWHEG, followed by
hadronization. LO+PS and fixed-order results should be retained as controlled
comparisons rather than mixed into the central prediction.

The initial analysis should concentrate on a perturbative DIS region with
explicit cuts on `Q^2`, `y`, and the hadronic final state. Low-`Q^2`
photoproduction, resolved photons, and multiparton interactions should remain
outside the baseline scope unless they are introduced and validated as a
separate extension.

Two beam-energy settings would give a useful first lever arm, for example
`5 x 41 GeV` and `18 x 275 GeV`, with an intermediate setting such as
`10 x 100 GeV` considered if it adds distinct coverage. These are starting
points for optimisation, not fixed machine assumptions. Every result should
state the beam energies, luminosity scenario, beam polarizations, and analysis
cuts explicitly.

### Physical helicity samples

The production campaign should generate four independently sampled unit-
helicity configurations,

```text
PP, PM, MP, MM,
```

where the first and second labels denote the lepton and proton helicities in a
documented convention. Each sample must be normalized before helicity
combinations are formed. For unit polarizations, define

```text
sigma_UU = (sigma_PP + sigma_PM + sigma_MP + sigma_MM) / 4,
sigma_LL = (sigma_PP - sigma_PM - sigma_MP + sigma_MM) / 4,
A_LL     = sigma_LL / sigma_UU.
```

The generator and PDF sign conventions must be checked against an analytic or
fixed-order benchmark before physics conclusions are drawn. Signed NLO
contributions must be accumulated and combined before ratios are formed.

Correlated optional event weights can be useful diagnostics for isolating
changes in matrix elements, Born spin density, or shower evolution. They must
not be described as substitutes for independently generated physical
helicity-conditioned samples.

## Measurement programme

The observables below form a hierarchy. The inclusive measurements establish
normalization and helicity conventions; jets and dijets provide the principal
polarized-PDF sensitivity; exclusive radiation observables test whether the
theory description is sufficiently controlled to support those measurements.

### 1. Inclusive polarized-DIS baseline

Measure or predict `sigma_UU`, `sigma_LL`, and `A_LL` in bins of

- `(x, Q^2)`;
- `(x, y)` at fixed or grouped `Q^2`;
- `Q^2` in several `x` regions.

The inclusive baseline tests the quark-helicity combinations, electroweak
implementation, normalization, and DIS reconstruction. It should include
event yields and statistical reach for clearly labelled benchmark luminosity
and polarization scenarios. Showing the polarized cross section alongside
the asymmetry is essential: cancellations in `A_LL` can otherwise conceal
large or mismodelled corrections in its numerator and denominator.

### 2. Breit-frame jets and dijets

The core of the study should be inclusive-jet and dijet spin asymmetries,
preferably beginning in the Breit frame. Candidate measurements include

- `A_LL^jet(pT_B, eta_B)` in broad `(x, Q^2)` regions;
- `A_LL^jet(x, Q^2)` with a resolved jet requirement;
- `A_LL^dijet(xi, M_jj, Delta eta, Delta phi)`;
- cross sections and asymmetries in boson--gluon-fusion- and
  QCD-Compton-enriched categories;
- matched lab-frame and Breit-frame definitions;
- a controlled jet-radius scan.

A useful partonic momentum estimator is

```text
xi = x_B * (1 + M_jj^2 / Q^2).
```

Its correlation with the incoming parton's momentum fraction should be
quantified at hard-process, showered, and hadron levels rather than assumed.
The study should determine whether `xi`, `M_jj`, or a two-dimensional
observable provides the best compromise between polarized-gluon sensitivity,
event rate, migration, and perturbative stability.

These measurements are the most directly relevant to polarized-PDF
determination. They need not be inserted into a fit to establish that
relevance: one can identify bins in which existing polarized-PDF predictions
separate more strongly than the combined non-PDF theory and anticipated
measurement uncertainties.

### 3. Exclusive radiation and the spin-aware shower

Exclusive observables are valuable both as measurements in their own right
and as theory-control observables for the jet analysis. Candidate quantities
include

- the helicity-resolved three-to-two jet fraction and `A_LL^3j`;
- the third-jet transverse momentum or hardness relative to the leading
  dijet system;
- the azimuth of the third jet relative to the hard-dijet plane;
- cumulative `cos(2 phi)` moments;
- energy--energy correlations in selected DIS frames;
- infrared-safe groomed angularities or radial jet profiles, provided that
  detector robustness can be demonstrated.

The principal generator comparisons should be

```text
LO+PS        versus NLO+PS,
parton level versus shower level versus hadron level,
full spin    versus spin-disabled evolution,
full spin    versus Born-density-only diagnostics.
```

The existing Herwig Polarized study found evidence that the Born spin-density
initialization can affect exclusive radiation, while it did not resolve a
separate effect from the accepted real-emission shower increment. The EIC
extension should therefore make the full-spin versus spin-disabled comparison
primary. Any interpretation of smaller differences between internal spin
modes should follow the observed numerical resolution rather than be assumed
in advance.

### 4. Reconstruction and measurement readiness

At least three standard DIS reconstructions should be compared where
applicable:

- electron reconstruction;
- Jacquet--Blondel reconstruction;
- double-angle reconstruction.

For each leading observable, compare truth-level and reconstructed
definitions, quantify migrations, and repeat the comparison with QED radiation
enabled and disabled. A lightweight response parametrization may be used to
study robustness, but it must be labelled as a sensitivity model rather than
an EIC detector simulation.

The useful outcome is not simply the smallest truth-level uncertainty. It is a
ranking of measurements that retain polarized-PDF discrimination after QCD
radiation, hadronization, QED radiation, and realistic kinematic
reconstruction are applied.

## A sensitivity study without a PDF projection

The analysis can answer PDF-relevant questions without producing a projected
PDF set. For each observable bin, compare predictions made with existing
polarized PDF sets or replicas and retain the sources of uncertainty
separately. An indicative pairwise separation measure is

```text
S_i(a,b) = |A_i(a) - A_i(b)|
           / sqrt(delta_theory,i^2 + delta_exp,i^2).
```

Here `delta_theory` excludes the polarized-PDF difference being tested and
`delta_exp` is an explicitly stated expected measurement uncertainty. This
quantity is only a bin-ranking diagnostic. It is not a likelihood, does not
account automatically for inter-bin correlations, and must not be presented
as a fit constraint.

The main synthesis could instead be an **observable-quality map**. Each bin or
measurement is classified as one or more of

- polarized-PDF dominated;
- perturbatively limited;
- shower sensitive;
- nonperturbative;
- reconstruction limited;
- statistically limited;
- measurement ready.

This provides a concrete deliverable for both phenomenology and experimental
planning: it identifies which measurements merit a more complete detector
study or eventual inclusion in a global analysis, and which require better
theory control first.

## Uncertainty budget

The following components should be evaluated and displayed separately before
any combined band is formed:

- polarized-PDF uncertainty or replica spread;
- unpolarized-PDF uncertainty;
- renormalization- and factorization-scale variation;
- POWHEG radiation-scale variation;
- supported shower starting-scale and recoil variations;
- hadronization and tune dependence;
- jet-radius dependence;
- QED-radiation dependence;
- reconstruction or response-model dependence;
- Monte Carlo statistical uncertainty.

Beam polarimetry and relative-luminosity uncertainties belong in the expected
measurement precision, not in the generator theory band. Correlations and
cancellations between `sigma_UU` and `sigma_LL` should be preserved wherever
the variation procedure allows them.

## Validation strategy

Validation should proceed from inclusive quantities to exclusive ones:

1. Verify helicity conventions and the inclusive polarized DIS cross section
   against an analytic or independent fixed-order calculation.
2. Validate inclusive jet and dijet distributions at fixed order, including a
   comparison with available polarized DIS jet calculations.
3. Establish the fixed-order, LO+PS, and NLO+PS hierarchy before studying
   nonperturbative corrections.
4. Compare hard-process, showered, and hadron-level results for every proposed
   headline observable.
5. Check the four independent helicity samples against correlated diagnostic
   weights only after their separate statistical and normalization properties
   are understood.
6. Require closure tests for helicity algebra, signed-weight handling, empty
   denominators, bin migration, and reproducibility metadata.

The historical `MC_DIS_BREIT` material identified in the repository audit is
a useful source of observable definitions and generator checks, but it should
be treated as a seed rather than inherited as an experimental measurement.
Its cuts, frame definitions, normalization, and weight semantics must be
audited before reuse.

## Minimum figure set

A compact first paper could be built around the following figures:

1. EIC `(x, Q^2)` coverage and event yields for the selected energy settings.
2. Inclusive `sigma_UU`, `sigma_LL`, and `A_LL` in representative bins.
3. Jet and dijet `A_LL` versus `pT_B` and `xi` or `M_jj`.
4. Fixed-order, LO+PS, and NLO+PS comparisons.
5. Shower and hadronization correction factors for the headline observables.
6. Full-spin, Born-density-only, and spin-disabled comparisons for an
   exclusive radiation observable.
7. Reconstruction and QED-radiation migrations.
8. The observable-quality map, including the spread from existing polarized
   PDFs and the non-PDF uncertainty budget.

The paper should quote expected statistical precision for one or more
transparent benchmark scenarios, but should stop before pseudodata fitting,
reweighting, or a new PDF error band is produced.

## Implementation roadmap

### Phase 0: definition audit

- Audit the historical `MC_DIS_BREIT` analysis and any associated cards.
- Fix the helicity, frame, jet, and reconstruction conventions in writing.
- Select the baseline cuts and two beam-energy settings.
- Identify an independent fixed-order benchmark.

### Phase 1: minimal generator study

- Add a dedicated EIC Rivet analysis and Herwig cards.
- Generate low-statistics `PP`, `PM`, `MP`, and `MM` smoke samples.
- Validate normalization and helicity combinations.
- Produce inclusive and leading-jet baseline distributions.

### Phase 2: fixed-order and shower validation

- Add dijet and exclusive-radiation observables.
- Compare fixed order, LO+PS, and NLO+PS.
- Test spin-aware, Born-density-only, and spin-disabled configurations.
- Establish parton-to-shower and shower-to-hadron corrections.

### Phase 3: uncertainty campaign

- Add scale, shower, hadronization, jet-radius, and PDF ensembles.
- Compare energy settings and optimize the final binning.
- Produce the first observable-quality map.

### Phase 4: reconstruction and reach

- Implement electron, Jacquet--Blondel, and double-angle variants.
- Quantify QED-radiation and response-model migrations.
- Add benchmark luminosity and polarization scenarios.
- Finalize the measurement ranking and paper-level figures.

A possible repository layout, to be created only when implementation begins,
is

```text
analyses/rivet/dis/        EIC Rivet analysis and plot metadata
cards/future-eic/          beam-energy and helicity-specific Herwig cards
config/future-eic/         observable, variation, and campaign registries
campaigns/                 ignored generated products
```

Tests should cover the helicity algebra, per-sample normalization, ordering of
signed NLO combinations, zero or statistically empty denominators, and the
configuration signatures needed to reproduce each result.

## Focused follow-up: charged-current DIS

A strong but distinct second study would be

> **Polarized charged-current DIS at the EIC: NLO+PS predictions and hadronic
> reconstruction.**

For `e- p -> nu + X`, the absence of a measured final-state lepton makes the
Jacquet--Blondel variables, missing transverse momentum, and QCD/QED migrations
central to the physics analysis. Single-spin charged-current asymmetries can
provide complementary flavour separation to neutral-current DIS. This topic
is sufficiently different in reconstruction and electroweak structure that it
should be a separate focused paper or a later phase, rather than an extra
section added to the initial jet study.

If supported by the available generator and nuclear inputs, proton and
effective-neutron targets could then be compared. Such an extension would
require its own assessment of nuclear and spectator-tagging assumptions.

## Explicit non-goals

The initial study should not claim to provide

- a polarized-PDF fit, reweighting exercise, profiling, or projected PDF set;
- a full EIC detector simulation or detector-response-equivalent prediction;
- a TMD, transversity, or fragmentation-function extraction;
- an experimentally validated description of unresolved low-`Q^2`
  photoproduction or multiparton interactions;
- physical helicity samples constructed solely from correlated optional
  weights;
- a resolved accepted-real-emission spin effect unless it is numerically
  observed;
- validation of an existing polarized PDF merely because it is compared with
  its own fitted data.

These boundaries keep the central claim precise: the study identifies robust,
PDF-relevant polarized EIC measurements and establishes the theory control
needed before a true global-fit or detector-level projection is attempted.

## Starting references

- EIC Yellow Report: [arXiv:2103.05419](https://arxiv.org/abs/2103.05419).
- Herwig Polarized phenomenology study:
  [arXiv:2606.23845](https://arxiv.org/abs/2606.23845).
- Polarized inclusive DIS jet calculations at NNLO:
  [arXiv:2005.10705](https://arxiv.org/abs/2005.10705).
- Charged-current polarized EIC phenomenology:
  [arXiv:1309.5327](https://arxiv.org/abs/1309.5327).
- Recent global polarized-PDF context:
  [arXiv:2503.11814](https://arxiv.org/abs/2503.11814).

These references are entry points. The implementation phase should add the
specific fixed-order, reconstruction, detector-performance, and EIC beam-
scenario references used to define the final campaign.
