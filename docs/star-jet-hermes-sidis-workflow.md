# STAR jet and HERMES SIDIS workflow

## Scope

Three measurements are available through the phenomenology campaign:

- `STAR_2021_I1850855`: inclusive-jet and dijet $A_{LL}$ at 200 GeV;
- `STAR_2022_I1949588`: inclusive-jet and topology-resolved dijet $A_{LL}$
  at 510 GeV;
- `HERMES_2019_I1698889`: identified charged-hadron SIDIS
  $A_\parallel$ on proton and deuteron targets.

The default outputs are truth- or hard-parton-level data comparisons with the
generator qualifications documented below; they are not detector-response-
equivalent publication analyses.  Alternate STAR projections, closure
observables, the current unpolarized HERMES azimuthal moment, and MPI-on
stable-particle jets are diagnostics and require explicit options.

## Reference data

The STAR inputs vendor all 21 tables from HEPData record
`10.17182/hepdata.104836.v1` and all 20 tables from
`10.17182/hepdata.114778.v1`.  Each record, table, row count, table name,
coordinate, uncertainty label, and checksum is validated before normalized
JSON and Rivet reference YODA are generated.  The 200 GeV primary covariance
contains the 36 independent points from Tables 4--7; the 11-point combined
inclusive projection in Table 8 is display-only.  The 510 GeV covariance is
the complete $63\times63$ matrix for 14 inclusive and 49 dijet points.

At 510 GeV, HEPData reports corrected parton coordinates with horizontal
jet-energy-scale uncertainties, not the event intervals used to select the
jets.  The exact operational intervals from the publication are held in the
separately checksummed `analysis-binning.json`.  Rivet fills those intervals;
postprocessing verifies them and maps the estimates to contiguous display
cells around the published parton coordinates.  The original horizontal
uncertainties remain attached to the reference points.

The HERMES publication has no HEPData record.  Its official APS supplemental
ZIP is vendored with a checksum and a complete 57-file inventory.  The parser
rejects HTML or publisher challenge pages, missing files, reordered grids,
row-count changes, and covariance-dimension changes.  If APS blocks the
download, use:

```bash
python3 scripts/run_phenomenology_campaign.py fetch-data \
  --measurement HERMES_2019_I1698889 \
  --source-file /path/to/DB12134_supplemental.zip
```

## STAR physics definition

Both STAR predictions use 200 or 510 GeV longitudinally polarized proton
beams, four independent physical helicity samples, NNPDF4.0 NLO unpolarized
PDFs, and NNPDFpol2.0 NLO polarized PDFs.  The hard scattering is nevertheless
LO and is labelled `LO+PS`.

The full QCD $2\to2$ matrix element contracts its complex helicity amplitudes
with both incoming spin-density matrices, retaining color-flow interference.
The same contracted terms select diagrams and color flows, and the selected
amplitude tensor is attached to the hard vertex for the spin-correlated
shower.  The unpolarized limit uses $\rho=I/2$.  The hard-scale control maps
$\mu/\mu_0=(0.5,1,2)$ to `ScalePreFactor` $(0.25,1,4)$.

The nominal jets are clustered from terminal descendants of outgoing lines of
the explicit signal-process vertex before hadronization.  Beam-remnant
descendants are excluded and the event is rejected if this provenance cannot
be established.  QCD shower spin correlations, hadronization, and remnants
remain active, while MPI and QED showering are disabled.  The optional
`hadron_mpi_on` family instead clusters stable particles with MPI enabled and
is excluded from default goodness-of-fit results.

At 200 GeV the analysis uses anti-$k_T$ with $R=0.6$, the published central
and forward inclusive regions, and the two leading jets for dijets with
$p_{T,1}>8$ GeV, $p_{T,2}>6$ GeV, $|\eta|<0.8$, and
$\Delta\phi>120^\circ$.  Same-sign and opposite-sign pseudorapidity
topologies are filled separately.

At 510 GeV it uses anti-$k_T$ with $R=0.5$, inclusive jets with
$|\eta|<0.9$, and dijets with $p_{T,1}>7$ GeV, $p_{T,2}>5$ GeV,
$\Delta\phi>120^\circ$, and $|\Delta\eta|<1.6$.  The A--D categories are
forward--forward same sign, forward--central, central--central, and
forward--backward.

For each bin,

\[
\sigma_{UU}=\frac{PP+PM+MP+MM}{4},\qquad
\sigma_{LL}=\frac{PP+MM-PM-MP}{4},\qquad
A_{LL}=\frac{\sigma_{LL}}{\sigma_{UU}}.
\]

The full point-to-point covariance is used for goodness-of-fit.  Relative
luminosity is an additive global nuisance and beam polarization is a
multiplicative global nuisance.  Empty denominators are masked.  Parity,
single-spin, and beam-exchange closures are diagnostic outputs.

For the 510 GeV data, the published per-point systematic already contains the
absolute relative-luminosity uncertainty of (4.7\times10^{-4}), whereas the
published correlation matrices exclude it. The normalized point-to-point
systematic therefore subtracts that component in quadrature before building
the covariance. The positive-semidefinite point-to-point covariance is then
augmented only by Monte Carlo statistical variance, and the luminosity and
polarization nuisances are profiled once. The primary fit restricts the
published 63-point covariance to 59 points: inclusive analysis bins 1--4 are
retained as diagnostic-only outputs, while inclusive bin 5 is the first
primary bin and is validated to start at \\(p_T=13.1\\) GeV. All dijet points
remain primary.

The nominal 4 GeV generator cut is gated by independent 50-million-event per
helicity samples at 3 and 4 GeV. Inclusive bins 5 and 6 and the first two
finite bins of each dijet topology must agree for both `SigmaUU` and `A_LL`
within the larger of the fixed tolerance and three combined Monte Carlo
standard errors. The 5 GeV sample is a non-gating stress test. Distinct scan
cuts must use disjoint initial seeds and the same immutable campaign source
commit. The report separately records the checker, measurement descriptor,
and comparison-policy hashes. Existing scan summaries are consumed read-only;
they are not re-postprocessed under the changed measurement signature.

## HERMES SIDIS physics definition

The SIDIS cards generate a 27.6 GeV positron on a fixed proton or neutron
through pure photon exchange.  They use the NNPDF4.0 NLO and NNPDFpol2.0 NLO
sets, polarized POWHEG NLO+PS, QCD shower spin correlations, hadronization,
and remnants.  QED showering and MPI are disabled because the measurement is
Born-level unfolded.

Kinematics are reconstructed from the fixed-target event:
$Q^2=-q^2$, $x=Q^2/(2P\cdot q)$, $y=P\cdot q/(P\cdot k)$,
$W^2=(P+q)^2$, $z=P\cdot p_h/(P\cdot q)$, as well as $P_{hT}$,
$x_F$, and the laboratory lepton angle.  The selection requires
$Q^2>1$ GeV$^2$, $W^2>10$ GeV$^2$, $y<0.85$, the circular truth-level
lepton envelope $0.04<\theta<0.22$, $x_F>0.1$, and normally
$0.2<z<0.8$.  This is not the rectangular, year-dependent HERMES
spectrometer acceptance, and no hadron angular acceptance is applied.  The
published low-$z$ extension is used only for the applicable $x$--$z$
projection.

Every accepted identified hadron is filled.  Proton targets use
$\pi^\pm$ with $4<p_h<13.8$ GeV.  Deuteron targets use
$\pi^\pm,K^\pm$ with $2<p_h<15$ GeV.  The one-dimensional, coarse
two-dimensional, and nine-by-three-by-three grids are retained.  The main
publication gives the last three-dimensional split as \(x=0.4\).  The APS
text-file headers instead print 0.45, but the reported last-bin means cluster
near \(x=0.44\); the workflow records that discrepancy and uses the
publication value 0.4.

POSNLO and NEGNLO are combined within each helicity and target component
before asymmetries are formed.  Proton and neutron components give

\[
\sigma_{UU}^{d}=\frac{\sigma_{UU}^{p}+\sigma_{UU}^{n}}{2},\qquad
\sigma_{LL}^{d}=
  \frac{\sigma_{LL}^{p}+\sigma_{LL}^{n}}{2}.
\]

The HERMES deuteron asymmetries already include the published
\(1/f_D\) correction with \(f_D=0.926\).  The free-proton-plus-free-neutron
theory comparison therefore does not apply that factor a second time.

Charge differences are formed from normalized yields,

\[
A_\parallel^{h^+-h^-}=
\frac{\sigma_{LL}^{h^+}-\sigma_{LL}^{h^-}}
     {\sigma_{UU}^{h^+}-\sigma_{UU}^{h^-}},
\]

not by subtracting asymmetries.  The primary reference is the published
$A_\parallel$.  Charge-difference $A_1$ is auxiliary.  The current
helicity-summed $2\langle\cos\phi\rangle$ output is an unpolarized yield
moment and does not implement the HERMES spin-asymmetry cosine amplitude
fitted from ten $A_\parallel(\phi)$ bins; it must not be overlaid on that
published moment.
The one-, two-, and three-dimensional projections overlap and are therefore
reported separately; they are never combined into one goodness-of-fit value
without a cross-projection covariance.  The combined tabulated point
systematics already include polarization effects, so the quoted 6.6% proton
and 5.7% deuteron polarization uncertainties are metadata rather than
additional profiled nuisances.

## Campaign commands

Load the polarized Herwig environment, then validate the vendored inputs:

```bash
source /opt/homebrew/opt/modules/init/bash
module purge
module load herwig/pol

python3 scripts/run_phenomenology_campaign.py fetch-data \
  --measurement STAR_2021_I1850855
python3 scripts/run_phenomenology_campaign.py fetch-data \
  --measurement STAR_2022_I1949588
python3 scripts/run_phenomenology_campaign.py fetch-data \
  --measurement HERMES_2019_I1698889
```

Run central smoke checks:

```bash
for measurement in \
  STAR_2021_I1850855 \
  STAR_2022_I1949588 \
  HERMES_2019_I1698889
do
  python3 scripts/run_phenomenology_campaign.py full \
    --measurement "$measurement" \
    --tag smoke_20260722 \
    --smoke --jobs 4
done
```

Run the nominal paper grid without optional modeling families:

```bash
python3 scripts/run_phenomenology_campaign.py full \
  --measurement STAR_2022_I1949588 \
  --tag paper_star510 \
  --profile paper \
  --polarized-pdf-members all \
  --unpolarized-pdf-members all \
  --scales all \
  --lo-events 300000 \
  --shards 10 --jobs 4
```

Use `--families hadron_mpi_on` to generate only that optional STAR modeling
family, or `--families all` for all registered families.
`--include-diagnostics` adds closures and alternate reference projections.
`--plot-comparisons` overlays explicitly generated non-nominal families.
The terminal tracker, resume behavior, immutable manifest checks, and
`--recover-failed` fresh-seed recovery are shared with the other
phenomenology measurements.

## Validation status

The active library has passed central card reads and smoke campaigns for all
three measurements, including postprocessing and data-overlay plots. Focused
tests cover the reference inventories and checksums, covariance construction,
all campaign axes, SIDIS arithmetic, random physical density matrices, and
the exact pre-change unpolarized color/helicity sums for all eight QCD
subprocess topologies. A separate external or analytic fixed-phase-space
helicity oracle for the legacy `MEQCD2to2` class is still a desirable
validation hardening step; no discrepancy is currently known. Production
statistics are intentionally not part of the implementation validation.
