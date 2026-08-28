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

The raw yield histograms include the selected dijet rate, the first three jet
$p_T$ spectra, $(p_{T,1}+p_{T,2})/2$, $p_{T,2}/p_{T,1}$,
$p_{T,3}/p_{T,1}$, dijet mass, dijet $p_T$, $p_{T,1}+p_{T,2}$,
$\Delta\phi$, $|\Delta\eta|$, $\Delta R$, $(\eta_1+\eta_2)/2$, and
$\cos(2\Delta\phi)$. Third-jet observables require $p_{T,3}>2$ GeV.

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
`ALL_X`. `--include-diagnostics` also writes the two single-spin combinations
and the `PP`--`MM` and `PM`--`MP` parity residuals. Statistical errors assume
the four helicity samples and the two shower families are independent.

The checked `data/phenomenology/MC_POLDIJETS/analysis-definition.json` file is
an internal observable contract, not experimental data. Consequently this
measurement has no reference YODA, data overlay, pull, or goodness-of-fit.

## Commands

First run the eight-sample smoke test:

```bash
module purge
module load herwig/pol

python3 scripts/run_phenomenology_campaign.py full \
  --measurement MC_POLDIJETS \
  --tag mc_poldijets_spin_smoke_20260828_v1 \
  --families nominal,shower_spin_off \
  --smoke --jobs 8 \
  --plot-comparisons --include-diagnostics
```

A useful first production has 100 million events per helicity and shower
treatment, i.e. 800 million generated events in total. With 500 shards per
logical sample, each shard contains 200,000 events:

```bash
python3 scripts/run_phenomenology_campaign.py full \
  --measurement MC_POLDIJETS \
  --tag mc_poldijets_spin_100m_20260828_v1 \
  --families nominal,shower_spin_off \
  --lo-events 100000000 \
  --shards 500 --jobs 100 \
  --seed-base 5107000 \
  --plot-comparisons --include-diagnostics
```

`--lo-events` is per helicity and per family. Always use a new immutable tag
if any event count, shard count, seed, family, cut, PDF, or analysis definition
changes. To recover failed shards without touching successful nonempty YODA
files, repeat the exact command with `--recover-failed`.
