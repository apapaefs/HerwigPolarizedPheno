# Harder particle-jet spin-transfer pilot

`MC_POLJETSHAPES_HARD` is a separate internal measurement. It leaves the frozen
`MC_POLJETSHAPES`, `MC_POLJETSHAPES_LHE`, STAR definitions and completed tags
unchanged. It runs the new particle-level analysis and the existing
`MC_POLDIJETS` companion on each generated event once.

## Physics and sampling

The setup remains LO longitudinally polarized pp at 510 GeV, independent
PP/PM/MP/MM unit-helicity theory samples, hadronization/decays/remnants on,
MPI off, all stable visible particles, anti-kT R=0.5, |eta|<1.5 and the loose
pT1>5, pT2>4 GeV baseline. NLO PDF inputs do not make the hard calculation NLO.

- Intra-jet primary/secondary and hard-plane angles use each jet's own
  half-open pT window: [20,30) or [30,45) GeV.
- Inter-jet primary angles and resolved jet-3/4 planes, spectra and radiation
  rates use the **leading** jet's window. There is no extra balance/topology cut.
- C/A declustering and signed plane transport match the frozen particle-level
  `MC_POLJETSHAPES` definitions. Standard primary/secondary splittings have
  z1,z2>0.1, with both kT>1 or both kT>2 GeV. Hard-sharing variants require
  0.25<z1<0.40, z2>0.35 with the same kT thresholds.
- Jet-3/4 spectra use pT>2 GeV. R32, R43 and third-jet veto scans retain
  2,3,4,5,6,8,10 GeV thresholds and the pT31/pT41 cumulative tails retain
  0.1,0.2,0.3,0.4,0.5 thresholds.

The two independent families differ only by `HardProcessSpin Yes/No`.
Both explicitly retain `SpinCorrelations Yes` and the same polarized hard
generation. Red is full spin, blue is LHE-like: the latter removes inherited
hard spin and polarized backward-ISR conditioning, not ordinary shower-created
spin correlations. This is not a claim of exact LHE equivalence.

The generation floor stays at 3 GeV so events migrating upward under showering
are retained. To improve sampling of harder jets, the common card attaches
`ThePEG::ReweightMinPT` (Power=4, Scale=20 GeV, OnlyColoured=Yes) to
`SubProcess:Preweights`. ThePEG includes this factor during sampling and divides
the event weight by it in `StandardEventHandler::select`. It does **not** attach
a physical `Reweights` factor. Rivet must retain the compensated event weights.
Both families use the identical bias. The saved descriptor/cards and runtime
ReweightMinPT library fingerprints are included in provenance.

Consequently cutflow/accepted-entry fractions are weighted, not counts of
sampled events. Their fraction times requested events is explicitly named as
such in JSON; production sizing uses Neff=(sum w)^2/sum w^2 and independent
moment errors. Do not merge these products with unweighted raw accumulators
or reinterpret them as equal-weight event counts.

## Launch on Odysseus

Use a persistent terminal. The existing thermal guard matches the shared-runner
alias below (`run_validation_campaign.py full`). The thin dedicated wrapper
`scripts/run_mc_poljetshapes_hard_campaign.py` also works, but its name does not
match that existing guard rule. Do not run both controllers for one tag.

```bash
module purge
module use /home/apapaefs/Projects/Herwig/production/hard-process-spin-20260908-v1/modulefiles
module load herwig/pol-hardspin-20260908-v1
cd /home/apapaefs/Projects/HerwigPolarizedPheno

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

python3 campaigns/control/rivet-fidelity-20260905/runtime/run_validation_campaign.py full \
  --measurement MC_POLJETSHAPES_HARD \
  --tag mc_poljetshapes_hard_pilot_500k_20260910_v1 \
  --families nominal,lhe_like_shower \
  --lo-events 500000 --shards 50 --jobs 100 \
  --seed-base 9407000 \
  --plot-comparisons --include-diagnostics --plot-jobs 16
```

This requests 500,000 sampled events **per helicity/family**, 4 million total,
400 shards of 10,000 events, seeds 9407000–9407399. It is a bounded pilot, not
high-statistics production. Fresh run files are generated with the dedicated
hard-spin installation; neither the old `herwig/pol` installation nor ThePEG
needs rebuilding. Resume an interrupted pilot with the identical command and
`--recover-failed`; completed nonempty shard products are reused. Do not change
immutable settings such as the tag's jobs or event counts when resuming.

## Results and validation

Results are under
`campaigns/phenomenology/MC_POLJETSHAPES_HARD/<tag>/`:

- `plots/html/index.html`: complete gallery;
- `plots/html/focus/index.html`: harder-window moments/shapes, sine nulls,
  jet-3/4 spectra and rates, including full/LHE-like ratios or signed differences;
- `postprocess/summary.json`, `central.csv`, YODA predictions and the existing
  sensitivity/moment summaries;
- `postprocess/shard-block-covariance.json`: selected spectra/rates and angular
  moment covariances, estimated from equal-size common-shard blocks;
- `postprocess/statistics-projection.json`: bounded 100M/250M/500M tiers,
  requiring all four standard kT>1 jet/window selections to project at least
  250,000 effective entries in every helicity/family and independent A2UU/A2LL
  differences with uncertainty <=0.002. Missing/masked inputs cannot pass.
  Rare kT>2/hard-sharing/inter-jet limitations are reported separately.

If no tier passes, the largest bounded command is a ceiling suggestion only,
**not** a claim that it meets the targets. No production is started automatically.
The new production seed base is 9507000, disjoint from the pilot and validation.

After completion, audit both namespaces, weighted shape closure, rate
identities, primary acceptance, HTML and sampling normalization against the
existing unbiased sample:

```bash
python3 scripts/check_harder_jet_pilot.py \
  --campaign-dir campaigns/phenomenology/MC_POLJETSHAPES_HARD/mc_poljetshapes_hard_pilot_500k_20260910_v1 \
  --reference-summary campaigns/phenomenology/MC_POLJETSHAPES_LHE/mc_poljetshapes_lhe_1500m_20260908_v1_350/postprocess/summary.json \
  --output campaigns/phenomenology/MC_POLJETSHAPES_HARD/mc_poljetshapes_hard_pilot_500k_20260910_v1/postprocess/validation-audit.json
```

This reports weighted/unbiased ratios for both jets, windows and families,
including helicity spectra and pointwise independent errors. UU and helicities,
or different jet bins, must not be counted as independent tests of one global
normalization hypothesis. A sampling check is not an angular LHE closure test.

To regenerate the gallery without rerunning generation:

```bash
python3 scripts/run_mc_poljetshapes_hard_campaign.py plot \
  --tag mc_poljetshapes_hard_pilot_500k_20260910_v1 \
  --plot-comparisons --include-diagnostics --plot-jobs 16
```

Known polarized-ISR PDFmax/overestimate warnings, predominantly in bottom
channels, remain a separate generator-precision issue. Harder selections and
successful execution do not resolve that caveat. Interpret any isolated excess
against the full scan and sine nulls, with independent-family errors and
appropriate within-sample covariance, before claiming a spin-transfer signal.
