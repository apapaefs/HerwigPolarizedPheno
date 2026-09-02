# SIDIS transverse-momentum and azimuthal diagnostics

Three schema-6 analyses implement the actionable remainder of the SIDIS
assessment. They use the same registry, immutable manifests, sharding,
normalized POSNLO/NEGNLO addition, target-component combination,
postprocessing and plotting path as the other HerwigPolarizedPheno SIDIS
measurements.

| Analysis | Output | Reference status | Central jobs |
|---|---|---|---:|
| `COMPASS_2013_I1236358` | fitted low-$p_T$ $\langle p_T^2\rangle$ for $h^\pm$ in 23 $(x,Q^2)\times8z$ cells | 368 values reconstructed from checksum-pinned paper Tables 1--3 | P/N $\times$ 00 $\times$ POS/NEG = 4 |
| `COMPASS_2014_I1278730` | $A_{UU}^{\cos\phi_h}$ and $A_{UU}^{\cos2\phi_h}$ in the published 1D and 3D bins | 480 values reconstructed from checksum-pinned paper Tables 2--12 | P/N $\times$ 00 $\times$ POS/NEG = 4 |
| `HERMES_2013_I1111237` | direct $\langle\cos\phi_h\rangle_{UU}$ and $\langle\cos2\phi_h\rangle_{UU}$ for H/D, $h^\pm$, $\pi^\pm$ and $K^\pm$ in all 900 four-dimensional cells | exact data-free definition; the official values and covariance endpoint is currently unavailable | P/N $\times$ 00 $\times$ POS/NEG = 4 |

The analyses are generator-level diagnostics. They do not apply detector
acceptance, unfolding, QED radiative corrections or diffractive subtraction,
and they do not imply a controlled Cahn, Boer--Mulders or twist-3
calculation. The COMPASS $A_{LU}^{\sin\phi_h}$ result is deliberately absent:
it requires the polarized muon beam and cannot be formed from the explicit
zero-polarization samples.

## Observable and covariance contracts

`COMPASS_2013_I1236358` fills event-aggregated charged-hadron yields in 36
$p_T^2$ bins over the published fit range $0.1<p_T<0.85$ GeV. The runner adds
normalized signed-NLO contributions, forms the P/N isoscalar yield, and fits
$A\exp[-p_T^2/\langle p_T^2\rangle]$ independently in every $(x,Q^2,z)$ cell.
The inclusive-DIS denominator and $z$ width are constant within one fit and
therefore cancel exactly from its inverse slope.

`COMPASS_2014_I1278730` stores event-aggregated
$\sum2\cos(n\phi_h)$, $\sum\epsilon_n(y)$ and signed same-event covariance
proxies, with

\[
\epsilon_1=\frac{2(2-y)\sqrt{1-y}}{1+(1-y)^2},\qquad
\epsilon_2=\frac{2(1-y)}{1+(1-y)^2}.
\]

The P/N target sum is formed before the ratio. The paper's point-to-point
systematic uncertainty, twice the tabulated statistical uncertainty, is kept
separate in the reference snapshot.

`HERMES_2013_I1111237` uses the exact $5\times5\times6\times6=900$ cell
layout. Its moments use a hadron-count denominator rather than a COMPASS
depolarization denominator. Hydrogen uses the proton sample; deuterium uses
the explicit half-proton plus half-neutron combination. There are 12
target/species/charge groups and two moments, hence 21,600 potential numerical
values. None is synthesized while the official release remains unavailable.

## Commands

Load the same Herwig Polarized environment used for the existing SIDIS
campaigns, then enter the repository:

```bash
source ~/Projects/Herwig/Herwig-REAL-stable-gcc-full/bin/activate
cd ~/Projects/HerwigPolarizedPheno
```

Validate the frozen sources/definitions and rebuild the aggregate Rivet
plugin:

```bash
python3 scripts/run_phenomenology_campaign.py fetch-data \
  --measurement COMPASS_2013_I1236358
python3 scripts/run_phenomenology_campaign.py fetch-data \
  --measurement COMPASS_2014_I1278730
python3 scripts/run_phenomenology_campaign.py fetch-data \
  --measurement HERMES_2013_I1111237

make rivet
make check-rivet
```

Inspect the complete four-job matrices without creating or running a
campaign:

```bash
for analysis in \
  COMPASS_2013_I1236358 \
  COMPASS_2014_I1278730 \
  HERMES_2013_I1111237
do
  python3 scripts/run_phenomenology_campaign.py full \
    --measurement "$analysis" \
    --tag "${analysis}_plan" \
    --profile central \
    --dry-run
done
```

Run one 100-event structural smoke per logical job, including compilation,
generation, postprocessing and plots:

```bash
for analysis in \
  COMPASS_2013_I1236358 \
  COMPASS_2014_I1278730 \
  HERMES_2013_I1111237
do
  python3 scripts/run_phenomenology_campaign.py full \
    --measurement "$analysis" \
    --tag "${analysis}_smoke_20260903" \
    --smoke \
    --jobs 4
done
```

After a successful smoke, the descriptors' starting central sample is three
million POSNLO and 0.3 million NEGNLO events per target. For example:

```bash
python3 scripts/run_phenomenology_campaign.py full \
  --measurement COMPASS_2014_I1278730 \
  --tag compass2014_azimuth_central_v1 \
  --profile central \
  --posnlo-events 3000000 \
  --negnlo-events 300000 \
  --shards 30 \
  --jobs 8
```

Use a new tag when changing any physics or statistics option. To recover only
failed shards with fresh seeds, repeat the identical command and add
`--recover-failed`. The HERMES kaon cells and the finely divided COMPASS slope
cells are expected to need a statistics assessment before any physics claim;
the descriptor defaults are a starting point, not a publication-quality
guarantee.

No event campaign is launched by the implementation or by `fetch-data`,
`make`, or `--dry-run`.
