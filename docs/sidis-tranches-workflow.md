# SIDIS first/second tranches and high-z charge ratios

Six schema-6 SIDIS measurements extend the existing COMPASS 2009/2017 and
HERMES 2019 coverage:

| Analysis | Published target | Observable | Central matrix |
|---|---|---|---:|
| `COMPASS_2026_I3096394` | isoscalar | corrected `h`, `pi`, `K` `dM/dz` | P/N x 00 x POS/NEG = 4 |
| `COMPASS_2025_I2840545` | hydrogen | `h`, `pi`, `K` `dM/dz` | P x 00 x POS/NEG = 2 |
| `COMPASS_2010_I862410` | proton | `A1` for identified `pi`, `K` | P x PP/PM/MP/MM x POS/NEG = 8 |
| `HERMES_2013_I1208547` | H/D | five-binning `pi`, `K` multiplicities | P/N x 00 x POS/NEG = 4 |
| `COMPASS_2018_I1624692` | isoscalar | `d2M/(dz dPhT2)` | P/N x 00 x POS/NEG = 4 |
| `COMPASS_2020_I1788430` | isoscalar | high-z antiproton/proton and K-/K+ ratios | P/N x 00 x POS/NEG = 4 |

The versioned sources, normalized snapshots, audits and checksums are under
`data/phenomenology/<analysis>`. Reference YODA is regenerated from the
normalized snapshot by `fetch-data` or automatically during `prepare`; it is a
runtime product and is not committed for the 2020 addition. Rebuild generated
sparse-cell C++ only from the validated snapshots:

```bash
python3 scripts/generate_sidis_tranche_binning.py
python3 scripts/generate_sidis_tranche_cards.py
python3 scripts/generate_sidis_tranche_descriptors.py
```

## Schema-6 target and density contract

`cards.target_components` is a nonempty subset of `P,N`.
`postprocess_config.target_outputs` maps each published target to explicit
component weights. Hydrogen is `P:1`; isoscalar/deuteron outputs use
`P:0.5,N:0.5`. A dataset declares the widths applied after its numerator/DIS
ratio: `z`, `z,pt2`, or `z,phperp`. Polarized `A1` datasets also declare the
longitudinal target scale. The target convention also applies to the
schema-5 multiplicities.

All multiplicity estimators use event-aggregated numerator, denominator and
same-event covariance objects. Shards are summed inside each signed NLO
contribution. POSNLO and NEGNLO normalized bins are then added, target yields
are combined, and the final ratio and density widths are applied. Generic raw
`yodamerge` arithmetic is not a valid final estimator.

Since 2026-09-05, the 2025/2026 multiplicities impose the published
z-bin-dependent nu window on both numerator and denominator, using the
species mass (pion mass for unidentified hadrons). The RICH vertex-angle
surrogate has been removed from generated-level selection. The current
releases have support in every cell at the nominal beam energy; unsupported
cells in the older 2017 releases are explicitly masked. Use new campaign
tags for these corrected selections. See
[the correction record](rivet-fidelity-corrections-20260905.md).

The HERMES full archive defines 6,592 retained cells across five overlapping
3D binnings. Each target/species/binning has its own released statistical
covariance; goodness of fit is reported separately. Readable projections are
formed from integrated numerator and denominator yields rather than averaging
preformed multiplicities. The archive covariance files contain ordinary
transpose-level differences and nine isolated one-sided decimal/exponent
corruptions. The normalized snapshot averages valid transpose pairs and, for
those nine entries, retains the partner that obeys the covariance
Cauchy--Schwarz bound. Raw and normalized matrix hashes and every replacement
are frozen in the reference audit; no diagonal element is changed.

The complete COMPASS 2018 HEPData v1 source contains 4,664 rather than 4,918
rows. The source-count audit in `reference.json` freezes that discrepancy.

No official HEPData record exists for `COMPASS_2020_I1788430` as of
2026-08-28. Its checksum-pinned arXiv TeX/PDF is therefore the numerical
authority. The TeX source of the earlier kaon-ratio paper is pinned as well,
because the 2020 publication explicitly inherits its kaon selection. All 67
paper-table rows are audited. Reconstructed-z bin limits are applied to truth
`z`; published corrected-z means are retained for display. The paper quotes
systematic correlations of about 0.7--0.8, modelled within each separately
reported table with correlation coefficient 0.75. The overlapping x-z and
z-momentum tables are never combined into one goodness of fit.

Charge ratios use the same ordering guarantees as multiplicities, but their
event-aggregated negative- and positive-hadron yields replace the hadron/DIS
pair. Their same-event covariance is retained, P/N yields are combined first,
and the negative/positive ratio is formed last. No density width is applied.

## Campaigns

Use the runtime-locked controller documented in
`campaigns/control/sidis-tranches-20260828/README.md`. Standard interactive
defaults remain 300,000 POSNLO and 30,000 NEGNLO events with one shard;
`--smoke` uses 100 events. These defaults are structural only for the large
sparse-cell analyses. Physics-quality production must explicitly use the
controller-generated event plan or larger overrides.

The statistics gate also treats a populated experimental bin with both zero
prediction and zero Monte Carlo variance as unresolved. Such a zero-entry
pilot bin cannot support `1/sqrt(N)` extrapolation and therefore triggers the
immutable doubled-statistics pilot plan.

The pilot/central commands are:

```bash
CONTROL=/home/apapaefs/Projects/HerwigPolarizedPheno-sidis-tranches-20260828/campaigns/control/sidis-tranches-20260828
"$CONTROL/run.sh" dry-run pilot all
"$CONTROL/run.sh" prepare pilot all
"$CONTROL/run.sh" launch pilot all
"$CONTROL/run.sh" postprocess pilot all
"$CONTROL/run.sh" assess pilot all

"$CONTROL/run.sh" prepare central all
"$CONTROL/run.sh" launch central all
"$CONTROL/run.sh" postprocess central all
"$CONTROL/run.sh" plot central all
"$CONTROL/run.sh" package central all
```

Production is never launched by validation, preparation, assessment or dry
runs. Paper profiles are prepared/launched only per analysis.
