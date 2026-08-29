# MC_POLDIJETS shower-spin workflow

`MC_POLDIJETS` is an internal Rivet measurement for a statistically efficient
comparison of the polarized Herwig shower with spin correlations enabled and
disabled. It is not a STAR data comparison. It reuses the STAR 510 GeV beam,
PDF, hard-process, MPI, and spin-density setup while replacing the published
STAR selections by a deliberately loose, MC-oriented dijet selection.

## Selection and sample definition

The default analysis clusters anti-$k_T$ jets with $R=0.5$ from terminal
hard-scatter shower partons before hadronization. Beam-remnant descendants are
excluded by the same provenance implementation used for the STAR jet
analyses. The event selection is

- $p_{T,1}>5$ GeV and $p_{T,2}>4$ GeV;
- $|\eta_j|<1.5$;
- no same-side/opposite-side category, $\Delta\phi$ cut, or $\Delta\eta$ cut;
- MPI off and a 3 GeV generator-level `JetKtCut`.

The nominal and `shower_spin_off` families use the same polarized hard
process, NNPDFpol2.0 and NNPDF4.0 inputs, independent `PP`, `PM`, `MP`, and
`MM` samples, analysis selection, and 510 GeV collision energy. Their Monte
Carlo streams have disjoint seeds. The only configured physics switch between
the families is
`/Herwig/Shower/ShowerHandler:SpinCorrelations`; the off family is styled blue.
Calling it a spin-averaged or unpolarized *shower* does not make the hard
process or beams unpolarized.

Theory cross sections use unit beam helicities, as in the STAR 510 campaign.
STAR's measured beam-polarization magnitudes are not multiplied into a
parton-level theory asymmetry.

## Observables

The raw yield histograms include the selected dijet rate, $p_T$ and $\eta$ for
the first, second, third, and fourth jets, $(p_{T,1}+p_{T,2})/2$,
$p_{T,2}/p_{T,1}$, $p_{T,3}/p_{T,1}$, dijet mass, dijet $p_T$,
$p_{T,1}+p_{T,2}$, $\Delta\phi$, $|\Delta\eta|$,
$\Delta R_{12}$, $\Delta R_{13}$, $\Delta R_{24}$, $\Delta R_{14}$,
$(\eta_1+\eta_2)/2$, and $\cos(2\Delta\phi)$. Third-jet observables require
$p_{T,3}>2$ GeV. Fourth-jet observables, $\Delta R_{24}$, and
$\Delta R_{14}$ require $p_{T,4}>2$ GeV. The full $p_T$ histogram ranges are
retained in YODA; the default plots stop at 70, 45, 20, and 15 GeV for jets
one through four, respectively, to avoid a gallery dominated by empty tails.

For every raw observable $X$, normalized independent-helicity histograms are
combined bin by bin as

\[
 \sigma_{UU}=\frac{\sigma^{++}+\sigma^{+-}+\sigma^{-+}+\sigma^{--}}{4},
 \qquad
 \Delta\sigma_{LL}=\frac{\sigma^{++}-\sigma^{+-}-\sigma^{-+}+\sigma^{--}}{4},
\]

and

\[
 A_{LL}=\frac{\sigma^{++}-\sigma^{+-}-\sigma^{-+}+\sigma^{--}}
 {\sigma^{++}+\sigma^{+-}+\sigma^{-+}+\sigma^{--}}.
\]

The primary YODA objects are named `SigmaUU_X`, `DeltaSigmaLL_X`, and
`ALL_X`. For the four jet-$p_T$ spectra, four jet-$\eta$ spectra, and four
requested $\Delta R$ spectra, the postprocessor also writes the directly
normalized helicity cross sections `SigmaPP_X`, `SigmaPM_X`, `SigmaMP_X`, and
`SigmaMM_X`. These are $\sigma^{++}$, $\sigma^{+-}$, $\sigma^{-+}$, and
$\sigma^{--}$, not polarization-diluted detector yields.
`--include-diagnostics` additionally writes the two single-spin combinations
and the `PP`--`MM` and `PM`--`MP` parity residuals. Statistical errors assume
the four helicity samples and the two shower families are independent.

The checked `data/phenomenology/MC_POLDIJETS/analysis-definition.json` file is
an internal observable contract, not experimental data. Consequently this
measurement has no reference YODA, data overlay, pull, or goodness-of-fit.

## Commands

First run the eight-sample smoke test with a new immutable tag:

```bash
module purge
module load herwig/pol

python3 scripts/run_phenomenology_campaign.py full \
  --measurement MC_POLDIJETS \
  --tag mc_poldijets_helicity_jets_smoke_20260829_v1 \
  --families nominal,shower_spin_off \
  --smoke --jobs 8 \
  --plot-comparisons --include-diagnostics
```

A 500,000-event-per-sample pilot was used to project the statistical precision.
In the well-populated ranges $p_{T,1}<45$ GeV, $p_{T,2}<30$ GeV,
$p_{T,3}<15$ GeV, $p_{T,4}<10$ GeV, and away from the last angular endpoint
bins, 100 million events per helicity and shower treatment gives binwise
relative errors below about 3% for the direct helicity cross sections, except
for the fourth-jet $p_T$ spectrum where the worst bin is about 5%. The
corresponding worst absolute $A_{LL}$ errors are about 0.02 for fourth-jet
$p_T$, 0.013 for $\Delta R_{24}$, and below 0.01 for the other requested
spectra.

For a balanced production, use 250 million events per helicity and family.
This is 2 billion generated events in total and projects to about 3.2% in the
limiting fourth-jet helicity bin and an absolute $A_{LL}$ error of about 0.013.
With 500 shards per logical sample, each shard contains 500,000 events:

```bash
python3 scripts/run_phenomenology_campaign.py full \
  --measurement MC_POLDIJETS \
  --tag mc_poldijets_helicity_jets_250m_20260829_v1 \
  --families nominal,shower_spin_off \
  --lo-events 250000000 \
  --shards 500 --jobs 100 \
  --seed-base 7107000 \
  --plot-comparisons --include-diagnostics
```

For a quicker first look, change only the immutable tag, use
`--lo-events 100000000`, and choose a new disjoint seed base; this is 800
million events total. For publication-quality fourth-jet and large-angle
spectra, 500 million events per helicity and family (4 billion total) reduces
the limiting projected errors to about 2.2% on a helicity cross section and
0.009 on $A_{LL}$. The extreme $p_T$ and $\Delta R$ endpoint bins may still
need merging: they were too sparsely populated in the pilot for a reliable
extrapolation.

These estimates describe statistical readability, not guaranteed separation
of the spin-on and spin-off shower models. Resolving a model difference of
size $\delta$ requires the combined statistical uncertainty of the two
independent families to be comfortably smaller than $\delta$. The pilot is
too small to determine that effect size without biasing the estimate.

`--lo-events` is per helicity and per family. Always use a new immutable tag
if any event count, shard count, seed, family, cut, PDF, or analysis definition
changes. To recover failed shards without touching successful nonempty YODA
files, repeat the exact command with `--recover-failed`.

The completed `mc_poldijets_spin_100m_20260828_v1` YODAs predate these newly
booked histograms. They cannot be postprocessed into the jet-rank and direct
helicity spectra, so this extension requires a fresh production tag.
