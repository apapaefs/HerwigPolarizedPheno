# Audit of the `HwPolNotesNew` data-comparison studies

> **Historical-disposition note.**  This document records why studies were
> retained or replaced.  It is not a certification of paper-level
> compatibility.  The later
> [experimental-analysis compatibility audit](experimental-analysis-compatibility-audit.md)
> supersedes its compatibility claims and documents result-changing repairs,
> invalidated products, and irreducible detector/generator differences.

## Scope and decision rule

This audit compares the historical analyses in `HwPolNotesNew` with the
version-pinned experimental workflow used for the phenomenology paper.  A
historical study is retained only when it corresponds to a published
measurement that is not already represented by the same observable, binning,
and fiducial definition.  Theory-only and generator-validation distributions
remain useful, but they are not presented as experimental comparisons.

The nominal predictions always use independently generated unit-helicity
samples, labelled (PP,PM,MP,MM).  Cross sections are normalized and shards
are combined before helicity sums are formed.  For DIS at NLO, POSNLO and
NEGNLO are also added at normalized-bin level before any ratio is taken.

## Summary

| Historical study | Experimental content | Audit result | Final disposition |
|---|---|---|---|
| HERMES 2007 | Table 7 Born (A_\parallel(x,Q^2)), 45 points in 19 (x) slices | Same publication as the primary Table 14 analysis, but a different projection and a lower (Q^2) reach | Retained as `HERMES_2007_I726689_LEGACY`; low-(Q^2) points are marked and excluded from quantitative tests |
| COMPASS 2010 proton | Table 1 v2 (A_1^p(x)), 15 points at nominal 160 GeV | Distinct from the later 200 GeV proton and 160 GeV deuteron data already present | Added as `COMPASS_2010_I843494` |
| STAR 2019 (W^\pm/Z) | (W^\pm A_L), publication-only (W^\pm A_{LL}), and integrated (Z/\gamma^* A_L) | Distinct polarized-(pp) measurements; historical implementation was incomplete | Rebuilt as `STAR_2019_I1708793` |
| PHENIX 2023 photons | Inclusive and isolated direct-photon cross sections and isolated-photon (A_{LL}) | Distinct polarized-(pp) measurement; historical implementation covered only a leading isolated photon | Rebuilt as `PHENIX_2023_I2033856` |
| `MC_WJETPOL` | Fiducial (W+)jet theory distributions | No matched experimental table or detector-level result was identified | Excluded from the data-comparison paper; retained as future theory work |
| `MC_DIS_BREIT` | Breit-frame DIS dijet and validation distributions | Generator/POLDIS validation rather than a matched experimental result | Excluded from the data-comparison paper; retained as future theory work |

## HERMES 2007

### Relationship to the primary HERMES comparison

Both analyses refer to the HERMES inclusive longitudinal proton data in
Phys. Rev. D 75 (2007) 012007, but they are not duplicates at observable
level.

| Feature | Primary comparison | Historical projection retained here |
|---|---|---|
| HEPData table | Table 14 | Table 7 |
| Observable | (A_1^p(x)) | Born (A_\parallel^p(x,Q^2)) |
| Number of points | 15 | 45 |
| Organization | one-dimensional (x) bins | 19 (x) slices with one to three (Q^2) bins |
| (Q^2) range | (1<Q^2<20\,\mathrm{GeV}^2) | (0.18<Q^2<20\,\mathrm{GeV}^2) |
| Depolarization treatment | inverse-(D) weighting, R1990, with (\eta A_2) neglected | direct (A_\parallel=\sigma_{LL}/\sigma_{UU}) |
| Quantitative status | primary exploratory comparison at the PDF boundary | measured data; only the generator curve below (Q^2=1\,\mathrm{GeV}^2) is a display-only extrapolation |

The historical Rivet code booked Table 7 correctly and reproduced the main
lepton cuts, but it accumulated cross sections rather than an asymmetry and
assumed that an external two-sample combination would supply the spin
observable.  The rebuilt version reconstructs the fixed-target laboratory
kinematics locally and uses all four physical helicities.  It applies
(0.18<Q^2<20\,\mathrm{GeV}^2), (0.1<y<0.91),
(W^2>3.24\,\mathrm{GeV}^2), and (0.04<\theta<0.22).  The implemented
per-(x) (Q^2)-sub-bin boundaries were not found in the paper, HEPData, or
official archive during the later compatibility audit and therefore remain
unverified.

The PDF, factorization, renormalization, and POWHEG-emission scales are frozen
at 1 GeV whenever the generated scale would be lower.  This prescription is a
controlled extrapolation, not a claim of perturbative reliability.  The eight
points whose published mean (Q^2) is below (1\,\mathrm{GeV}^2) are
identified in the reference data, visually distinguished, and omitted from
pulls, goodness-of-fit summaries, and numerical conclusions.  The published
statistical and systematic errors, including the stated normalization content
of the systematic uncertainty, remain separate in the vendored record.

## COMPASS proton measurements

The 2010 proton result is distinct from both current COMPASS analyses:

| Dataset | Beam and target | Observable and bins | Principal selection |
|---|---|---|---|
| COMPASS 2010 | nominal 160 GeV (mu^+p) | 15-bin (A_1^p), (0.004<x<0.7) | (Q^2>1\,\mathrm{GeV}^2), (0.1<y<0.9) |
| COMPASS 2016 | 200 GeV (mu^+p) | 17-bin (A_1^p) extending to (x=0.0025) | published high-energy proton selection |
| COMPASS 2017 | 160 GeV (mu^+d) | 15-bin (A_1^d) | deuteron impulse-approximation workflow |

The historical 2010 code omitted the published (Q^2>1\,\mathrm{GeV}^2)
cut and contained a dimensional typo in its (Q^2) conversion.  It also
stored a cross section in a histogram named (A_1).  The replacement uses
HEPData Table 1 version 2, finite-muon laboratory kinematics, R1998, and
neglected $(\eta A_2)$ without a separate $g_2$ model.  It imposes no $(W)$
cut because none is reported for this table.  A fixed 160 GeV beam is an
explicit approximation to the published
140--180 GeV beam-energy selection.

All three COMPASS datasets use NNPDF40 NLO for unpolarized densities and
NNPDFpol2.0 NLO for helicity differences.  The proton predictions use four
independent helicity samples.  The deuteron prediction keeps proton and
neutron target components separate through generation and combines them only
in the impulse approximation, including the configured effective
polarization.  Statistical, point-to-point systematic, Monte Carlo, PDF, and
hard-scale effects are not collapsed into a single error source.

## STAR (W^\pm) and (Z/\gamma^*)

The historical STAR analysis captured the basic lepton-(\eta) idea, but it
did not define prompt dressed leptons explicitly, contained a missing energy
unit on the upper transverse-momentum cut, and provided neither the
publication-only (A_{LL}) values nor the (Z/\gamma^*) result.  It also did
not implement the two-beam symmetrization required for a single-spin
asymmetry at a symmetric collider.

The replacement uses prompt dressed electrons and positrons with dressing
radius (R=0.1), (25<E_T<50\,\mathrm{GeV}), and six bins spanning
(-1.5<\eta<1.5).  HEPData Tables 17 and 18 supply the six (W^+) and six
(W^-) (A_L) points at their published coordinates.  The two beam
estimators are formed separately and combined as

\[
A_L(\eta)=\frac12\left[A_L^A(\eta)+A_L^B(-\eta)\right].
\]

(A_{LL}) uses all four helicity cross sections, with positive and negative
pseudorapidity bins folded at normalized-yield level.  The six
publication-only (W^\pm A_{LL}) points and the integrated
(Z/\gamma^*\) value (A_L=-0.04\pm0.07) are pinned to the checksum of the
publication PDF and record their page and table provenance.  The dilepton
selection requires isolated opposite-charge dressed electrons with
(E_T>14\,\mathrm{GeV}),
(|\eta|<1.1), and (70<m_{ee}<110\,\mathrm{GeV}).

This is a simplified truth proxy, not the STAR fiducial analysis: the later
audit identified different barrel/endcap regions and missing isolation,
signed-balance, opposite-side-energy, trigger, tracking, charge-sign, purity,
and background logic.

These are LO hard-process predictions with NNPDF40 NLO and NNPDFpol2.0 NLO
inputs.  The NLO PDFs do not make the (W/Z) matrix elements NLO accurate.
QCD and QED showers, spin-correlated decays, hadronization, remnants, and MPI
are included nominally.  MPI-off and unpolarized-beam samples are modeling and
closure variations.  The appropriate diagnostic is beam exchange,
(\sigma_{h_1h_2}(\eta)=\sigma_{h_2h_1}(-\eta)), rather than a parity
identity that weak interactions need not obey.

## PHENIX prompt photons

The historical PHENIX code selected only the leading prompt photon and only
the isolated spectrum.  It did not produce the inclusive spectrum, did not
apply the published charged and neutral detector thresholds, and did not
retain the global experimental uncertainties as named components.

The replacement fills every prompt photon with
(|\eta|<0.25) in the exact 18 Table 1 bins.  It provides inclusive and
isolated invariant cross sections using

\[
E\frac{\mathrm d^3\sigma}{\mathrm d p^3}
=\frac{1}{2\pi p_T\,\Delta\eta}
  \frac{\mathrm d\sigma}{\mathrm d p_T},\qquad \Delta\eta=0.5.
\]

Isolation uses (R=0.5), a cone energy below 10% of the photon energy, a
0.2 GeV/c charged-particle transverse-momentum threshold, and a 0.3 GeV
neutral-energy threshold.  The isolated (A_{LL}) is reported in all seven
Table 2 bins.  The 10% luminosity uncertainty on the cross sections, the
(3.9\times10^{-4}) relative-luminosity shift on (A_{LL}), and the 6.6%
polarization scale uncertainty are separate global components.

As for STAR, the prompt-photon hard process is LO despite the NLO PDF inputs.
The nominal full event simulation includes QCD and QED showers,
hadronization, remnants, and MPI; an MPI-off prediction is maintained as a
modeling variation.  Because this is a parity-conserving channel, the output
also contains (PP-MM), (PM-MP), and both single-spin closure diagnostics.
The generated direct-(gamma+jet) sample omits the separately normalized QCD
dijet fragmentation-photon event class, so these predictions are
diagnostic-only despite the close truth-observable implementation.

## Uncertainty and provenance contract

The HEPData payloads, normalized JSON, and Rivet reference YODA are vendored
with source-version identifiers and SHA-256 checksums.  Refreshing a source
requires the byte checksum, table identity, column headers, row count, and
normalized snapshot to remain unchanged; a revised upstream table therefore
fails rather than silently changing the analysis.

For the paper profile, a common hard-scale envelope uses
(\mu/\mu_0=0.5,1,2).  In electroweak Born production this is a
factorization-scale variation; in DIS and prompt-photon production it is the
supported correlated renormalization/factorization variation.  Polarized and
unpolarized NNPDF replicas are varied independently, using a common member in
all helicities and signed NLO components, and their replica variances are
added in quadrature.  Hard-scale, PDF, MPI, experimental-global, and Monte
Carlo uncertainties remain identifiable.  Empty denominators are masked, and
near-zero asymmetries are assessed with pulls rather than data/theory ratios.

## Excluded theory studies

`MC_WJETPOL` defines useful charged-lepton and jet distributions for a
polarized (W+)jet study, but the historical directory does not identify a
matching published dataset, fiducial correction, or covariance model.  It is
therefore not relabelled as a data comparison.

`MC_DIS_BREIT` is a detailed Breit-frame dijet and event-generator validation
analysis, including optional correlated event-weight diagnostics.  It is
scientifically distinct from the inclusive HERMES and COMPASS measurements
and has no matched experimental reference in this project.  Both studies are
kept as candidates for a later theory/validation paper.
