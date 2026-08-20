# Polarized data-phenomenology campaign

## Measurements

The primary entry point discovers both the established fixed-target DIS
registry and the polarized-$pp$ registry.  The implemented experimental
comparisons are:

- HERMES Table 14 $A_1^p(x)$ and the separate Table 7 low-$Q^2$
  $A_\parallel(x,Q^2)$ projection;
- COMPASS 2010 and 2016 proton $A_1^p$, and COMPASS 2017 deuteron
  $A_1^d$;
- HERMES identified-hadron SIDIS $A_\parallel$ for proton and deuteron
  targets;
- STAR $W^\pm A_L$, $W^\pm A_{LL}$, and integrated
  $Z/\gamma^* A_L$;
- STAR inclusive-jet and dijet $A_{LL}$ at 200 and 510 GeV;
- PHENIX inclusive and isolated prompt-photon cross sections and isolated
  $A_{LL}$.

The older experimental entry point remains valid for existing fixed-target
campaigns.  The phenomenology entry point is preferred for new work because
it scans both registries and records process, channel, family, helicity,
perturbative contribution, PDF member, hard scale, MPI state, and shard in a
job identity.

## Profiles and sample labels

`central` is the default profile.  It runs only the nominal physics family:
polarized POWHEG NLO+PS for DIS and SIDIS, and polarized LO+PS for RHIC.
The default STAR jet prediction uses hard shower-parton jets with MPI off;
the established weak-boson and prompt-photon definitions retain their
measurement-specific nominal MPI settings.

The fixed-target comparison profile retains these scientific labels:

- `NLO+PS (polarized; full spin)`;
- `NLO+PS (unpolarized beams)`;
- `LO+PS (polarized)`;
- `NLO+PS (polarized; Born spin only)`;
- `NLO+PS (polarized; shower spin off)`.

The `paper` profile controls the PDF-member and hard-scale grids; it does not
enable unrelated modeling families.  The same PDF member and scale are used
in every physical helicity, target component, and signed NLO contribution.
Use `--families all` or a comma-separated family list to generate optional
modeling samples.  `--comparisons` remains a backward-compatible alias for
`--profile paper --families all`.  Plots remain nominal-only unless
`--plot-comparisons` is requested.

For STAR jets, the optional `hadron_mpi_on` family is labelled
`LO+PS (polarized; hadron level, MPI on)`, and `unpolarized_closure` provides
the unpolarized hard-parton closure check.  These samples are diagnostics and
are excluded from the default goodness-of-fit calculation.  All RHIC hard
predictions remain explicitly labelled LO even though NLO PDF inputs are used.

## Typical commands

After loading a clean polarized environment:

```bash
source /opt/homebrew/opt/modules/init/bash
module purge
module load herwig/pol
```

list and validate the available measurements:

```bash
python3 scripts/run_phenomenology_campaign.py list
python3 scripts/run_phenomenology_campaign.py fetch-data \
  --measurement STAR_2019_I1708793
```

run a central smoke campaign:

```bash
python3 scripts/run_phenomenology_campaign.py full \
  --measurement PHENIX_2023_I2033856 \
  --tag phenix_smoke \
  --smoke --jobs 4
```

The new jet and SIDIS measurements use the same interface:

```bash
python3 scripts/run_phenomenology_campaign.py full \
  --measurement STAR_2022_I1949588 \
  --tag star510_smoke \
  --smoke --jobs 4

python3 scripts/run_phenomenology_campaign.py full \
  --measurement HERMES_2019_I1698889 \
  --tag hermes_sidis_smoke \
  --smoke --jobs 4
```

inspect a paper-profile plan without launching it:

```bash
python3 scripts/run_phenomenology_campaign.py prepare \
  --measurement STAR_2019_I1708793 \
  --tag star_paper \
  --profile paper --dry-run
```

For a nominal STAR jet paper grid with the two PDF ensembles and hard-scale
points, omit `--families all`:

```bash
python3 scripts/run_phenomenology_campaign.py full \
  --measurement STAR_2021_I1850855 \
  --tag star200_paper \
  --profile paper \
  --polarized-pdf-members all \
  --unpolarized-pdf-members all \
  --scales all \
  --lo-events 300000 --jobs 4
```

STAR jet campaigns also pin the generator-level `JetKtCut` in the immutable
job identity and generated card. The nominal value is 4 GeV; only the audited
3, 4, and 5 GeV scan points are accepted:

```bash
python3 scripts/run_phenomenology_campaign.py prepare \
  --measurement STAR_2022_I1949588 \
  --tag star510-cut3 \
  --profile central --jet-kt-min-gev 3 \
  --lo-events 50000000 --shards 100 --jobs 100
```

The option is rejected for non-jet measurements. STAR postprocessing retains
the normalized unpolarized yield `SigmaUU_<observable>` in the signed summary
alongside each `A_LL`; those yields are used by the independent generator-cut
stability audit and are not experimental cross-section overlays.

Add `--families hadron_mpi_on` only when the stable-particle MPI diagnostic is
wanted.  `--include-diagnostics` exposes alternate reference projections and
closure observables during postprocessing and plotting.
Nested diagnostic namespaces are rendered below their own HTML index.  Any
per-object Rivet/YODA script-generation exception makes the plot command fail,
even if upstream catches it internally.  A low-statistics diagnostic whose
generated axis limits are non-finite is instead skipped explicitly in
`rivet-plot-scripts.log`; its missing bins are never replaced by zeros or an
invented display range.

The HERMES APS supplement can be refreshed from an already downloaded archive
if the publisher blocks automated access:

```bash
python3 scripts/run_phenomenology_campaign.py fetch-data \
  --measurement HERMES_2019_I1698889 \
  --source-file /path/to/DB12134_supplemental.zip
```

The same command and selectors apply to fixed-target measurements.  For
example, a compact HERMES card-validation matrix is

```bash
python3 scripts/run_phenomenology_campaign.py prepare \
  --measurement HERMES_2007_I726689 \
  --tag hermes_variation_readcheck \
  --profile paper \
  --polarized-pdf-members 0,7 \
  --unpolarized-pdf-members 0,9 \
  --scales 0.5,1,2 \
  --posnlo-events 1 --negnlo-events 1 --lo-events 1
```

PDF and scale selectors can restrict a validation run without changing the
profile definition:

```bash
python3 scripts/run_phenomenology_campaign.py full \
  --measurement STAR_2019_I1708793 \
  --tag star_variation_readcheck \
  --profile paper \
  --polarized-pdf-members 0,7 \
  --unpolarized-pdf-members 0,9 \
  --scales 0.5,1,2 \
  --lo-events 100 --jobs 4
```

An immutable manifest prevents options from being changed under an existing
tag.  Successful nonempty shard outputs are reused.  A failed shard is only
resubmitted when failed-shard recovery is requested, and receives a fresh
seed.  The terminal tracker, machine-readable status, generated cards, run
files, logs, YODA, summaries, and plots are written below the ignored campaign
area.

The full measurement, family, helicity, contribution, PDF, scale, and shard
identity remains in the manifest and campaign output path.  The token passed
locally to `Herwig run -t` is deliberately shortened to the shard and attempt
numbers, preventing Rivet and EvtGen output basenames from exceeding the
255-byte per-component filename limit on macOS.

The postprocessed YODA keeps Monte Carlo, PDF, hard-scale, and MPI effects as
named uncertainty components.  Paper plots render the Monte Carlo component
as a bin-aligned stepped band centred on the nominal prediction and the PDF
and hard-scale components as separate shaded bands, rather than replacing
them with a single quadrature band.  Experimental points use one total error
bar by default; `--plot-data-components` adds nested statistical-only bars.
Fixed-target plots show only the nominal
`NLO+PS (polarized; full spin)` prediction by default; passing
`--plot-comparisons` to `plot` or `full` explicitly overlays the additional
postprocessed families.  The low-$Q^2$ HERMES projection additionally marks
display-only points with open squares.

## Combination rules

Four physical helicity samples define every nominal asymmetry.  For DIS,
POSNLO and NEGNLO and then shards are combined at normalized-bin level before

\[
\sigma_{UU}=\frac{PP+PM+MP+MM}{4},\qquad
\sigma_{LL}=\frac{PP+MM-PM-MP}{4}
\]

are evaluated.  Within-sample covariance between the ordinary and
inverse-depolarization fills is retained in $A_1$.  A missing helicity or
signed NLO component stops postprocessing.

STAR single-spin asymmetries combine the two beam estimators after reflecting
the second beam in pseudorapidity.  STAR double-spin bins are folded only
after normalized yields have been formed.  PHENIX invariant cross sections
average the four unit-helicity cross sections; its $A_{LL}$ uses their
longitudinal difference.  Zero denominators are represented by masked bins.

The STAR jet analyses form $A_{LL}$ from normalized $PP$, $PM$, $MP$, and
$MM$ hard-parton yields.  The 510 GeV analysis fills the detector-level event
intervals published in the paper, then displays the estimates at the
HEPData-provided corrected parton coordinates and coordinate uncertainties.
This avoids treating horizontal jet-energy-scale uncertainties as event-bin
boundaries.

For HERMES SIDIS, proton and neutron POSNLO/NEGNLO yields are combined before
forming proton or deuteron asymmetries.  The published deuteron asymmetries
already contain the \(1/f_D\) correction with \(f_D=0.926\), so the free
proton-plus-neutron theory numerator receives no additional factor.
Charge-difference predictions subtract normalized positive- and
negative-hadron yields before taking the ratio; asymmetries are never
subtracted directly.  Overlapping one-, two-, and three-dimensional projections
retain separate covariance-aware goodness-of-fit results.

## Reproducibility limitations

NNPDFpol2.0 includes much of the HERMES, COMPASS, and STAR information used in
the comparisons.  Agreement therefore tests the event-generation and
analysis chain but is not an independent validation of the polarized PDF fit.
The HERMES Table 7 points below $Q^2=1\,\mathrm{GeV}^2$ are measured
Born-level data; only the corresponding generator predictions use a 1 GeV
scale floor and are display-only extrapolations.  The 160 GeV COMPASS
generation is an approximation to the published 140--180 GeV beam selection.  The DIS
$A_1$ predictions neglect the $\eta A_2$ contribution and do not separately
model $g_2$; the deuteron is treated in the configured impulse
approximation.  RHIC hard processes are LO even though NLO PDF sets are used.

The STAR jet prediction is a Herwig hard-shower-parton proxy, not a
publication-level, detector-level, or stable-particle jet prediction.  It
starts from outgoing lines of the explicit signal-process vertex and therefore
contains their final-state shower descendants but not initial-state shower
partons.  It excludes beam-remnant descendants, uses the detector-level event
bin edges as a proxy for the measurement's matched representative parton
coordinates, and fails an event whose hard provenance cannot be established.
The MPI-on stable-particle result is therefore a modeling diagnostic rather
than an uncertainty band on the primary prediction.

The HERMES APS tables provide combined point systematics.  Their quoted 6.6%
proton and 5.7% deuteron polarization uncertainties are retained as metadata,
but are not profiled again as independent nuisances because doing so would
double count them.  The published charge-difference $A_1$ values are
auxiliary diagnostics; the primary comparison is to $A_\parallel$.  The
current helicity-summed $2\langle\cos\phi\rangle$ output is an unpolarized
yield moment, not the HERMES spin-asymmetry cosine amplitude obtained by
fitting $A_\parallel(\phi)$, and must not be overlaid on the published
azimuthal result.
