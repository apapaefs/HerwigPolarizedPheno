# Experimental-analysis compatibility audit

Audit date: 2026-07-26; corrected-production protocol added 2026-08-19

## Scope and meaning of compatibility

This is the authoritative paper-to-code audit for the experimental
comparisons implemented in this repository. It covers the ten registry entries
matched to published measurements:

- HERMES 2007 inclusive proton `A1` and the retained 45-point
  `A_parallel` projection;
- COMPASS 2010 and 2016 proton `A1`, and COMPASS 2017 deuteron `A1`;
- HERMES 2019 identified-hadron SIDIS;
- STAR 2019 weak bosons, STAR 2021 and 2022 jets; and
- PHENIX 2023 direct photons.

`MC_DIS_BREIT`, `MC_DIS_PS`, and `MC_WJETPOL` are generator or theory
analyses, not implementations of an experimental measurement, and are
outside this compatibility claim.

Three different claims must not be collapsed into the phrase “compatible
with the paper”:

1. **Reference-data compatibility** means that the published central values,
   coordinates, errors, covariance, and provenance have been transcribed
   correctly.
2. **Truth-observable compatibility** means that the particle-level
   observable, binning, cuts, normalization, and helicity algebra match a
   quantity to which the experiment unfolded.
3. **Experimental-analysis compatibility** means that the trigger,
   acceptance, reconstruction, background subtraction, response correction,
   and correlated uncertainty model have also been reproduced or forward
   folded.

No implementation in this workspace can presently be certified as 100%
compatible in the third sense.  Several experiments publish unfolded
observables without a reusable detector response, trigger emulator, or full
systematic covariance.  The defensible target is therefore exact reference
data plus a clearly labelled truth-level comparison, with every
generator-dependent approximation exposed.

## Executive compatibility matrix

| Registry entry | Reference data | Truth observable | Experimental analysis | Audit disposition |
|---|---|---|---|---|
| `HERMES_2007_I726689` | Exact selected Table-14 values and 15 `x` bins | Conditional | Not reproduced | Regenerate with the prompt-lepton definition; retain only as a Born-level truth comparison with `eta A2` neglected |
| `HERMES_2007_I726689_LEGACY` | Exact selected Table-7 values | Unverified `Q2` cells | Not reproduced | Auxiliary only; obtain the authoritative 45-cell boundaries before claiming exact bin integration |
| `COMPASS_2010_I843494` | Exact HEPData v2 table | Correct after this audit | Not reproduced | **All pre-audit predictions are invalid and require regeneration** |
| `COMPASS_2016_I1357198` | Exact selected averaged `A1p` table, with a documented source typo correction | Conditional | Not reproduced | Regenerate with the prompt-lepton definition; retain as a fixed-energy truth proxy |
| `COMPASS_2017_I1501480` | Exact selected averaged `A1d` table | Conditional nuclear proxy | Not reproduced | Regenerate with the prompt-lepton definition and retain only with the free-`p+n` qualification |
| `HERMES_2019_I1698889` | Exact checksum-pinned APS tables and published statistical correlations | Conditional after the deuteron and sparse-3D repairs | Not reproduced | Regenerate all projections with the prompt lepton, including 3D at `x=0.4`; exact detector acceptance is absent |
| `STAR_2019_I1708793` | Central values are correct; complete combined covariance is unavailable | Simplified truth proxy | Not reproduced | Qualitative theory comparison only, not a literal STAR fiducial analysis |
| `STAR_2021_I1850855` | Exact HEPData inventory and 36-point primary covariance | Generator-dependent parton proxy | Not reproduced | Retain only as a hard-parton modeling comparison |
| `STAR_2022_I1949588` | Exact HEPData inventory; luminosity covariance repaired by this audit | Generator-dependent parton proxy | Not reproduced | Produce a newly signed campaign after the covariance repair; retain only as a modeling comparison |
| `PHENIX_2023_I2033856` | Exact selected tables and global-error metadata | Fiducial photon observable is close, generated process is incomplete | Not reproduced | Diagnostic only; exclude from publication-grade comparisons |

“Not reproduced” does not mean that detector cuts should be applied again to
an unfolded result.  It means that equivalence to the experiment's corrected
observable cannot be demonstrated without the response and correction chain.

## Result-changing defects found and repaired

### COMPASS 2010 `x` bins

The Rivet source inserted an unpublished edge at `x=0.014` and omitted
`x=0.500`.  Because both the incorrect and correct definitions contain 15
bins, ordinal plotting silently compared different intervals rather than
failing.  The corrected sequence is

```text
0.004, 0.005, 0.006, 0.008, 0.010, 0.020, 0.030, 0.040,
0.060, 0.100, 0.150, 0.200, 0.250, 0.350, 0.500, 0.700
```

Every `COMPASS_2010_I843494` prediction generated before this repair,
including `nominal_scales_300k_20260722` and
`paper_central_3m_20260723`, is incompatible with the published bins and
must not be plotted or quoted.  The event histograms cannot be repaired by
relabeling: the old definition split `0.010--0.020` and merged
`0.350--0.700`, so a fresh event run is required.

### HERMES SIDIS deuteron normalization

The published HERMES `A_parallel` values already divide the measured
deuterium asymmetry by the effective deuteron factor `f_D=0.926`.  A theory
calculation formed from free proton and neutron cross sections must therefore
be compared directly to those nucleon-level values.  The old postprocessor
multiplied the longitudinal free-`p+n` numerator by `0.926` a second time,
suppressing every deuteron prediction by 7.4%.  The extra theory factor has
been removed.  This arithmetic repair alone could reuse event histograms, but
the simultaneous prompt-lepton correction means that a compatibility-certified
result requires fresh events.  All deuteron postprocessing, figures, pulls,
and goodness-of-fit values must be regenerated.

### HERMES SIDIS sparse three-dimensional cells

The official three-dimensional APS tables use dense published row numbers for
species-dependent sparse subsets.  The Rivet histogram instead uses the full
canonical 81-cell grid

```text
theory_bin = ix * 9 + iz * 3 + ipT
```

with `ix=0..8`, `iz=0..2`, and `ipT=0..2`.  The old parser treated the
published row number as this canonical cell and consequently paired data
with the wrong Monte Carlo bin.  The normalized reference now preserves the
raw number as `published_row` and derives `flat_bin` from the row's mean
`x`, `z`, and `pT`.  All 477 three-dimensional points require regenerated
reference YODA, plots, pulls, and goodness-of-fit values.

The main paper gives the final three-dimensional `x` split as `0.4`, while
the APS text-file headers print `0.45`; the last-bin means near `x=0.44`
support the paper value.  The normalized record treats `0.45` as a documented
header typo and the Rivet grid now uses `0.4`.  Consequently pre-audit
three-dimensional event histograms also require a fresh run; they cannot be
repaired by postprocessing alone.

### HERMES SIDIS azimuthal observable

The optional generator diagnostic formed a helicity-summed unpolarized yield
moment but `--include-diagnostics` assigned it the reference path of the
published spin-asymmetry cosine-fit amplitude.  That made Rivet overlay two
different observables.  The published points are now labelled
`A_parallel_cosphi_amplitude`, with their literal APS source-column label
retained separately.  The generator output is labelled
`UnpolarizedCosPhi_*`, written only below the generator diagnostics path, and
has no reference path or overlay.  A data comparison remains disabled until
the ten-bin `A_parallel(phi)` fit is implemented.

### STAR 510 GeV relative luminosity

For the 510 GeV jet result the HEPData point systematic already contains the
absolute `4.7e-4` relative-luminosity component, while the published
correlation tables explicitly exclude relative luminosity and polarization.
The old normalized covariance used the full systematic and then profiled
`4.7e-4` again.  The repaired point-to-point component is

```text
sqrt(max(reported_systematic_total^2 - (4.7e-4)^2, 0))
```

and the luminosity nuisance is profiled exactly once.  The reported total is
retained as metadata.  The old jet histograms are arithmetically sufficient
for this reference-only repair, but the current immutable-signature guard
cannot certify a reference-only migration of a pre-audit manifest.  Until a
narrow verified migration tool exists, the safe supported workflow is a new
tag and event run; correlated goodness-of-fit summaries must be regenerated.

The corrected-production protocol additionally requires three independent
50-million-event-per-helicity generator-cut samples at 3, 4, and 5 GeV. The
3-versus-4 GeV comparison gates the nominal 4 GeV production, while 5 GeV is
reported as a stress test. This test addresses sensitivity to the generator
cut; it does not promote the hard-shower-parton observable to an experimental-
analysis-level jet definition.

### DIS PDF-replica variance

The fixed-target uncertainty path used deviations from PDF member zero
instead of deviations from the replica sample mean.  It now uses the sample
variance about the replica mean, consistently with the polarized-`pp`
campaign path.  Central-only production plots are unchanged, but any DIS PDF
band made with the old function must be regenerated.

### Scattered-lepton origin

The fixed-target Rivet analyses previously selected the highest-energy stable
same-PID lepton.  A fully decayed generator event can contain same-sign
leptons from hadron decays, so this did not encode the experimental
target-vertex requirement.  HERMES and COMPASS now select the highest-energy
**prompt** same-PID final-state lepton.  With QED showering disabled in the
Born-level comparison cards this is the appropriate generator-level signal
lepton convention.  Because old YODA shards do not retain alternative
candidate histories, numerical equivalence cannot be certified after the
fact: all pre-audit fixed-target production predictions must either be rerun
or pass a dedicated event-level old/new equivalence study.

### Immutable campaign compatibility

The postprocessors previously checked campaign signatures during preparation
and generation but not when reading completed shards.  A stale COMPASS 2010
histogram could therefore be relabelled with the corrected bins, and a
HERMES 3D histogram with the `0.45` boundary could be read against the new
`0.4` reference.  Experimental and phenomenology postprocessing and plotting
now reject any manifest whose analysis/reference/card signature differs from
the current definition.  Reusing arithmetic-compatible subsets requires a
new, explicitly designed conversion rather than a silent override.

The final diagnostic smoke test also exposed a Rivet/YODA plotting edge case:
the upstream script generator did not create intermediate directories for
nested paths such as `.../DIAGNOSTICS/UnpolarizedCosPhi_*`, then caught its own
`FileNotFoundError` and returned success after rendering only the primary
plots.  The local safe wrapper now creates the nested parents, resolves the
root plotting style at arbitrary namespace depth, and converts any suppressed
per-object generation exception into a nonzero command status.  A generated
low-statistics diagnostic with non-finite automatic axis limits is now
explicitly skipped rather than assigned invented limits or allowed to abort
the remaining plot set.

## Analysis-by-analysis findings

### HERMES 2007 inclusive proton

The `x` bins, `Q2`, `y`, `W2`, and angular cuts, R1990 implementation, and
inverse-depolarization estimator agree with the selected comparison.  The
data, however, were unfolded for detector smearing, external and internal QED
radiation, Bethe-Heitler background, and bin migration.  HERMES used a
world-data parametrization of `g2` in extracting `A1`, while the Monte Carlo
comparison evaluates `A_parallel/D` and thereby neglects the `eta A2` term.
It is therefore inaccurate to claim that the code models both `A2=0` and
`g2=0`: it does not calculate either quantity, but uses the explicit
approximation `eta A2 -> 0`.

The publication supplies non-diagonal statistical correlations from the
unfolding and a correlated normalization component.  The selected project
snapshot currently uses point errors only.  Table-14 values also involve the
publication's `Q2` combination/evolution convention, whereas the Monte Carlo
forms a cross-section-integrated bin ratio.  These differences preclude an
exact extraction-level claim.

The legacy 45-point `A_parallel` data are genuine low-`Q2` measurements; only
the Monte Carlo prediction below `Q2=1 GeV2` is an extrapolation with a scale
floor.  The hard-coded per-`x` `Q2` boundaries are not supplied in the paper,
HEPData, or official archive inspected in this audit and must be treated as
unverified until an authoritative bin map is located.

Primary source: [HERMES, Phys. Rev. D 75 (2007)
012007](https://arxiv.org/abs/hep-ex/0609039).

### COMPASS inclusive proton and deuteron

After the 2010 edge repair, the selected `x` bins and main cuts agree with the
three chosen averaged-`A1` tables.  The finite-muon depolarization and `eta`
formulae and R1998 parametrization agree with the publications.  The
comparison neglects `eta A2`, as done for the quoted approximation.  In the
2016 HEPData payload, row 5 prints the upper edge as `0.001`; the normalized
record verifies that raw value and applies the documented `0.001 -> 0.010`
correction implied by the `x=0.009` mean and the next bin.

The generator uses a fixed nominal beam energy, while the 160 GeV samples
accept 140--180 GeV and the 200 GeV sample accepts 185--215 GeV.  The
experiments use event weights involving beam polarization, dilution,
depolarization, acceptance, and time-dependent target cells.  Polarization
and dilution factors must **not** be multiplied into a unit-helicity theory
asymmetry, but their response-weighted kinematic sampling is not reproduced
by a fixed-energy truth integral.

The COMPASS deuteron target was polarized `6LiD`.  The implemented
`0.5(p+n)` model with the configured effective longitudinal polarization is a
free-nucleon impulse approximation; binding, Fermi motion, shadowing,
off-shell, tensor, and other nuclear effects are absent.  Published
systematic components are also reduced to per-point combined errors.

Primary sources:
[COMPASS 2010](https://arxiv.org/abs/1001.4654),
[COMPASS 2016 proton](https://arxiv.org/abs/1503.08935), and
[COMPASS 2017 deuteron](https://arxiv.org/abs/1612.00620).

### HERMES 2019 identified-hadron SIDIS

The checksum-pinned APS archive, point inventories, published statistical
correlations, `x/z/pT` grids, main DIS and hadron cuts, momentum/PID ranges,
and multiplicity-weighted yield construction are reproduced.  Signed NLO
components and independent helicities are combined in the correct order.
Overlapping projections remain statistically separate because HERMES did not
publish cross-projection covariance.

The ideal Rivet acceptance is not the HERMES spectrometer acceptance.  HERMES
retained a rectangular track region approximately
`|theta_x|<170 mrad` and `40<|theta_y|<140 mrad`; the implementation uses a
circular lepton envelope `0.04<theta<0.22` and no hadron angular cut.
Year-dependent acceptance and unfolding response are not available as a
reusable model, so exact fiducial equivalence remains unattainable.

The APS column labelled `2<cos(phi)>` is the cosine amplitude of the
spin-dependent asymmetry.  HERMES first formed `A_parallel(phi)` in ten
azimuthal bins and then fitted constant, `cos(phi)`, and `cos(2phi)` terms.
The implementation instead sums the helicity samples and forms the
unpolarized yield moment `2<cos(phi)>`.  These are different observables:
the current Fourier moment is diagnostic-only and must not be overlaid on,
or labelled as, the HERMES fitted moment.  The combined point systematic is
treated diagonally because its full covariance is unpublished; the quoted
6.6% proton and 5.7% deuteron polarization components are correlated, so
current projection goodness-of-fit values are approximate rather than exact.

Additional model limits are the free-`p+n` deuteron, absence of a dedicated
exclusive/diffractive vector-meson component, ideal truth PID, and unresolved
weak-decay track-origin convention.

Primary source: [HERMES, Phys. Rev. D 99 (2019)
112001](https://arxiv.org/abs/1810.07054).

### STAR 2019 weak bosons

The helicity algebra, beam exchange and `eta` reflection, prompt dressed
leptons, displayed `25<ET<50 GeV` condition, and the truth-level
`Z/gamma*` selection are reasonable theory definitions.  They do not
reproduce the STAR `W` analysis:

- STAR used four barrel regions plus one-sided endcap acceptance, not a
  contiguous six-bin `|eta|<1.5` truth histogram;
- local and cone isolation, signed transverse-momentum balance,
  opposite-side energy, trigger, tracking, charge-sign, purity, and
  background logic are absent; and
- the complete covariance for the combined data is not public.  The plotted
  HEPData point systematics omit the 3.3% polarization scale and other
  correlated content present in the publication's totals.

The resulting curves are simplified LO+PS truth predictions evaluated near
the published points.  Pulls using only statistical and plotted point
systematics are diagonal diagnostics, not a complete STAR goodness of fit.

Primary source: [STAR, Phys. Rev. D 99 (2019)
051102](https://arxiv.org/abs/1812.04817).

### STAR 200 and 510 GeV jets

The jet radii, main kinematic cuts, topology definitions, selected data
inventories, display coordinates, and published covariance blocks are
implemented carefully.  The current hard-parton object is nevertheless not
the STAR parton correction:

- STAR clustered hard-scattered partons including initial- and final-state
  radiation, excluding underlying event and beam remnants.  The helper seeds
  only outgoing signal-vertex lines, so it includes their final-state shower
  descendants but omits initial-state-radiation emissions.
- STAR selected triggered detector jets, applied underlying-event,
  reconstruction, and trigger-bias corrections, and quoted results at
  PYTHIA+GEANT-derived representative parton coordinates.  The code applies
  the detector event-bin intervals directly to Herwig parton jets.  Those
  intervals are not a published generator-independent parton-bin selection.
- The STAR parton correction is itself tied to PYTHIA 6 Perugia12,
  embedding, matching, and detector simulation.  “Parton level” is not a
  universal cross-generator object.

Exact equivalence would require forward folding through the STAR response or
a published interpolation prescription.  Until a remnant-safe ISR recovery
is validated, these analyses must be described as Herwig hard-shower-parton
proxies.  The predictions are LO+PS with NLO PDFs, not the NLO theory used in
global analyses.  The 4 GeV generator cut also requires stability checks for
the first analysis bins.

Primary sources:
[STAR 200 GeV jets](https://arxiv.org/abs/2103.05571) and
[STAR 510 GeV jets](https://arxiv.org/abs/2110.11020).

### PHENIX 2023 direct photons

The Rivet observable itself reproduces the selected truth definition well:
all prompt photons, `|eta|<0.25`, the exact 18 cross-section and seven
`A_LL` bins, the `1/pT` invariant-cross-section factor, and the published
isolation cone and particle thresholds.

The generated sample is not the published direct-photon process.  The cards
contain the direct `gamma+jet` matrix element plus a QED shower but omit the
QCD dijet event class in which a hard quark fragments or showers to a
non-decay photon.  Isolation reduces but does not eliminate that component.
The publication's NLO calculation also contains a photon fragmentation
function and a fragmentation scale.  Inclusive cross sections are therefore
definitely incomplete, and isolated cross sections and `A_LL` are not exact.

The recorded 10% cross-section luminosity, `3.9e-4` `A_LL` shift, and 6.6%
polarization scale are global uncertainties.  Generic diagonal pulls that
omit them must not be interpreted as a complete goodness of fit.  PHENIX
should remain excluded from publication-grade conclusions until a separately
normalized fragmentation-photon sample with overlap removal, or a validated
NLO prompt-photon implementation, is available.

Primary source: [PHENIX, Phys. Rev. Lett. 130 (2023)
251901](https://arxiv.org/abs/2202.08158).

## Monte Carlo comparison contract

The following rules apply to every publication comparison:

1. Compare Born-corrected DIS data to a QED-off Born-level lepton definition,
   and do not reapply experimental beam-polarization or dilution factors to
   unit-helicity theory.
2. State whether particle definitions are prompt, dressed, stable, or
   ancestry restricted.  A stable-particle PID match alone is not a signal
   lepton definition.
3. Keep a detector-level selection, an unfolded particle/parton target, and a
   generator proxy as separate concepts.  Reusing the same numeric bin edges
   does not make them equivalent.
4. Treat beam spectra, nuclear targets, exclusive channels, weak decays, QED
   radiation, MPI, and fragmentation as physics-model choices, not hidden
   implementation details.
5. Add normalized POSNLO and NEGNLO bins before forming asymmetries.  Form
   independent-helicity observables from independently normalized samples.
   Correlated-weight shower results remain shared-history estimators and must
   not be relabelled as independent helicity-conditioned showers.
6. Preserve published covariance and named global nuisances.  Do not add a
   global uncertainty to a point covariance and profile it again.
7. Do not combine overlapping projections without cross-projection
   covariance, and do not call a diagonal pull plot a full goodness of fit.
8. Label comparisons to data already included in NNPDFpol2.0 as closure
   tests, not independent PDF validation.
9. Keep generator-native Herwig scale variations distinct from the
   fixed-`Q2` central-scale convention used by POLDIS.

## Publication and regeneration gate

Before any audited comparison is used in a paper:

- rebuild all Rivet plugins and regenerate `.run` files after interface or
  persistence changes;
- rerun every fixed-target production comparison with the prompt scattered
  lepton definition, including COMPASS 2010 with corrected binning and HERMES
  SIDIS 3D with the corrected `x=0.4` boundary;
- produce a newly signed STAR 510 campaign with the repaired covariance (or
  first implement and validate a reference-only manifest migration);
- regenerate any DIS PDF bands made with the old replica variance;
- perform low-statistics end-to-end smoke tests with finite accepted bins;
- run generator-cut stability scans for the first STAR jet bins;
- remove “publication-level” or “exact experimental analysis” language from
  STAR jets and weak bosons;
- keep PHENIX diagnostic-only until process completeness is demonstrated;
  and
- report which covariance and global nuisances enter every quoted pull or
  goodness-of-fit value.

The reference-data and workflow unit suite is necessary but not sufficient:
before this audit it passed while COMPASS 2010 retained the same number of
incorrect bins.  Source-to-reference edge parity, sparse-cell mapping, and
global-nuisance de-duplication are now explicit regression-test targets.

## Validation record

The final audited tree passed 175 Python tests in the active Herwig/Rivet
environment.  All ten experimental Rivet analyses compiled together into one
plugin, and an independent review found no remaining priority-0, priority-1,
or priority-2 defect.  Source/JSON/YODA edge parity, all 801 HERMES normalized
points, sparse-cell uniqueness, and positive-semidefinite STAR covariance were
checked independently.

Fresh final-signature smoke campaigns completed all eight COMPASS 2010 shards
and all sixteen HERMES SIDIS shards, normalized their signed contributions,
and produced the primary HTML reports.  The HERMES diagnostics pass generated
170 scripts and rendered 165 nonempty PNG/PDF pairs.  Five partially masked
parity diagnostics were explicitly skipped for non-finite automatic limits;
all six finite `UnpolarizedCosPhi_*` diagnostics rendered below the nested
diagnostic index, while no published cosine-fit object or data overlay entered
those scripts.  The corrected 20-page manuscript also builds without undefined
citations/references, LaTeX errors, or overfull boxes.
