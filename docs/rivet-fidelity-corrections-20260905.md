# Rivet fidelity corrections — 5 September 2026

The five actionable findings from the review of 21 data-linked registry entries
have been implemented. The numerical experimental points, errors and bin edges
remain unchanged. This record distinguishes verified code corrections from
the remaining limits of generator comparisons.

## Corrected selections

`COMPASS_2017_I1444985`, `COMPASS_2017_I1483098`,
`COMPASS_2025_I2840545` and `COMPASS_2026_I3096394` now apply the same
z-bin-dependent restriction to their numerator and DIS denominator:

\[
\frac{\sqrt{12^2+m_h^2}}{z_{\min}}<\nu=
\frac{P\cdot q}{M}<\frac{\sqrt{40^2+m_h^2}}{z_{\max}}.
\]

Energies and masses here are in GeV. Pions and unidentified charged hadrons
use the pion mass; kaons use the kaon mass. This follows
[COMPASS 2025 §3.1](https://arxiv.org/abs/2410.12005), which states that the
restriction was used for the earlier measurements, and the corresponding
sections in the [2017 pion](https://arxiv.org/abs/1604.02695) and
[kaon](https://arxiv.org/abs/1608.06760) papers.

The unsupported 10–120 mrad cut on generated vertex angles has been removed.
The papers apply RICH entrance geometry to reconstructed acceptance samples
and kinematic cuts to generated samples. Interpreting that prescription as
correcting out RICH geometry is an inference from the published construction;
no magnetic transport or detector replay is claimed. The generated angle
remains available as a diagnostic.

With the nominal fixed 160 GeV beam, five old pion/hadron cells and four old
kaon cells per species have no support. They are explicitly masked in YODA,
CSV, JSON and goodness-of-fit calculations. Their status is
`outside_nominal_beam_support`, rather than a finite prediction. The 2025
and 2026 releases have no such unsupported cells. No beam profile is invented.

The [STAR 200 GeV](https://arxiv.org/abs/2103.05571) and
[510 GeV](https://arxiv.org/abs/2110.11020) inclusive selections now keep at
most two highest-pT eligible jets. Eligibility precedes the cap; the separate
dijet selection and primary 510 GeV 13.1 GeV threshold are preserved.

## Corrected estimators

COMPASS 2014 now stores 16 phi yields and the full-phi epsilon sum per
kinematic cell/harmonic. Postprocessing fits bin averages of
`1, cos(phi), cos(2phi), sin(phi)` on bins 1–14 (zero-based), excluding the
two bins adjacent to zero, then divides each cosine amplitude by mean epsilon.
The sine is a nuisance parameter. This implements the published function and
exclusions in [§§5–6](https://arxiv.org/abs/1401.6284).

COMPASS 2013 now fits the bin-integrated exponential directly to signed
yields in all 36 pT² bins over 0.1–0.85 GeV in pT. It never takes logs or
deletes negative bins. Fits require positive variance in all bins, nonsingular
covariance and two-sigma positive yield in each of four contiguous regions.
These support checks are declared Monte Carlo quality policy. Chi-squared,
degrees of freedom and a poor-shape flag (chi-squared/ndof > 5) are retained;
bad model agreement alone does not erase a fit.

Both estimators retain within-event covariance. A triangular histogram filled
with `sqrt(input_i * input_j)` stores `sum(w_event² input_i input_j)` in
its **variance**, including negative-NLO events. Target coefficients enter
covariances squared. Shards combine within a contribution; normalized POSNLO
and NEGNLO contributions are added, followed by target sums and finally the
fit. Phi errors also propagate the correlation with mean epsilon. Slope
errors use local linear propagation of the two fit parameters.

The [2015 COMPASS slope erratum](https://doi.org/10.1140/epjc/s10052-014-3255-y)
retains the slopes and conclusions. It identifies a trigger-acceptance bias
in integrated multiplicities, clarifies 5% point-to-point and up to 40%
normalization uncertainty, and says radiative corrections were not applied
but had negligible effects on the fitted shapes. No unprovided slope
systematic is synthesized.

The papers do not provide enough minimizer settings or fit-input covariance
to assert an exact replay of their numerical extraction. Pure-function and
signed-spectrum closure tests validate this implementation. Closure by
refitting released experimental spectra remains outstanding; attempted
HEPData retrieval for record 61432 was unavailable during this work.

## Reproduction and migration

Load the matching Herwig/Rivet environment, then run:

```sh
make test
make rivet check-rivet
make check-rivet-fidelity
```

The dedicated fixture checks run actual plugins on momentum-conserving events:
three eligible jets give two inclusive entries; a leading jet outside eta does
not consume a slot; the COMPASS nu counterexample is rejected in both yields;
a valid 5 mrad hadron survives; two hadrons in one event reproduce every
stored within-cell covariance entry. The C++ fixture also checks both nu
boundaries and invariance under a common boost. Numerical tests cover exact
harmonics, excluded bins, finite-difference error propagation, target sums,
integrated exponentials, negative bins, sparse tails and singular covariance.

The clean test suite generates the intentionally untracked COMPASS 2020
reference YODA in a temporary directory, so it needs no prior fetch-data run.

Validation completed with `herwig/pol`, Rivet 4.1.2 and YODA 2.1.2:

- All 24 plugin registrations compile and load; 208 Python tests pass.
- All nine event fixtures pass, including covariance serialization to YODA.
- Fresh 100-event-per-job campaigns pass generation and postprocessing for
  COMPASS 2013 (four jobs, 48 prediction objects), 2014 (four jobs, 112 objects)
  and 2025 (two jobs, 246 objects). The 2025 HTML/plot workflow also completes.
- Every one of the 6,348 experimental points in the six affected reference
  snapshots has identical coordinates, values and uncertainties.

The slope smoke masks all 368 cells; the azimuthal smoke masks 478 of 480.
This is expected with such small samples. The smoke verifies the execution
and missing-information behavior, not statistical convergence or data agreement.
Installed executable/library paths and hashes are recorded in the isolated
campaign manifests; no installed Herwig library was rebuilt.

Use **new immutable campaign tags for all eight changed plugins**. Their
analysis/reference/support-file signatures change, and the runner refuses
to postprocess old manifests against the new definitions. Old flat YODA lacks
the omitted selection information or new fit inputs. Reanalyse retained
events or generate fresh campaigns; changing a signature by hand cannot fix
old predictions. Historical campaigns are retained as provenance.

## Remaining experimental qualifications

The STAR outgoing-hard-parton traversal still does not establish equivalence
to the experiment's ISR-plus-FSR jet definition. HERMES 2007 legacy Q² cells
remain inferred; PHENIX lacks a separate hard-QCD fragmentation-photon
sample. Fixed beam energies, nuclear effects, diffractive/weak-decay content,
detector response and unavailable covariance remain measurement-dependent
limits. These corrections do not promote those comparisons to precision
experimental reproductions, or change the project conventions for correlated
RIVETWEIGHTS showers and generator-native POWHEG scales.
