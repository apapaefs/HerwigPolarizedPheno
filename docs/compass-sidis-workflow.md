# COMPASS SIDIS analyses and campaigns

Three COMPASS SIDIS measurements are implemented by the primary campaign
runner:

These schema-5 analyses remain frozen for reproducibility. The corrected 2026
isoscalar multiplicities are the primary modern comparison; see
`docs/sidis-tranches-workflow.md` for the schema-6 tranche implementations.

| Analysis | Pinned release | Observable | Central jobs |
|---|---|---|---:|
| `COMPASS_2009_I820721` | [HEPData v1](https://www.hepdata.net/record/ins820721?version=1) | $A_{1,d}^{\pi^\pm,K^\pm}(x)$ | 16 |
| `COMPASS_2017_I1444985` | [HEPData v1](https://www.hepdata.net/record/ins1444985?version=1) | $dM^{\pi^\pm,h^\pm}(x,y,z)/dz$ | 4 |
| `COMPASS_2017_I1483098` | [HEPData v1](https://www.hepdata.net/record/ins1483098?version=1) | $dM^{K^\pm}(x,y,z)/dz$ | 4 |

The complete v1 record metadata and every table are vendored under
`data/phenomenology/<analysis>/raw/`. Each source manifest pins the record
DOI, version, table inventory, row counts, headers, and SHA-256 checksums. The
normalized JSON and compressed reference YODA files are generated from these
sources. HEPData v1 is the numerical authority for the 2009 result. The
checksum-pinned TeX paper table is retained only in
`paper-hepdata-discrepancy-audit.json`, which lists all 88 field-level
differences; no paper number replaces a HEPData value.

## Truth definitions

All three analyses reconstruct Lorentz-invariant DIS variables from the
incoming and prompt scattered positive muon and the incoming proton or
neutron. The fixed-target beam energy is 160 GeV. Hadron momentum and polar
angle are evaluated in the laboratory frame.

`COMPASS_2009_I820721` applies
$Q^2>1\,\mathrm{GeV}^2$, $0.1<y<0.9$, $0.004<x<0.3$,
$0.2<z<0.85$, and $10<p_h<50\,\mathrm{GeV}$. Event-aggregated
ordinary, inverse-depolarization, and ordinary--inverse-depolarization
covariance-proxy yields are filled for $\pi^\pm$ and $K^\pm$ in the ten
published $x$ bins. The final estimator uses the finite-muon-mass COMPASS
depolarization factor, R1998, and $\eta A_2=0$. Its target model is

\[
  \sigma_{UU}=\frac{p+n}{2},\qquad
  \sigma_{LL}=0.925\,\frac{p+n}{2}.
\]

The two 2017 analyses apply $Q^2>1\,\mathrm{GeV}^2$,
$W>5\,\mathrm{GeV}$, $0.004<x<0.4$, $0.1<y<0.7$,
$0.2\le z\le0.85$, $12<p_h<40\,\mathrm{GeV}$, and
$10<\theta_h<120$ mrad. They fill only the 311 pion/hadron or 309 kaon cells
released in HEPData, including the restricted $y$ coverage. Every cell has a
hadron numerator, matching inclusive-DIS denominator, and same-event
numerator--denominator covariance proxy. Unidentified $h^\pm$ means stable
charged hadrons, with the experimental pion-mass assumption used in $z$.

The proton and neutron components are summed before the numerator/DIS ratio,
and the result is divided by the cell width in $z$. Flattened-cell objects are
written for machine use and 38 per-$(x,y)$ $z$ slices per species are written
for readable plots. Missing cells and non-positive denominators are masked.
The reference is the final radiative- and diffractive-vector-meson-corrected
multiplicity. Correction columns remain in JSON and the central CSV as
provenance; no detector, radiative, or experimental correction is applied to
Herwig.

## Covariance policy

The 2009 statistical covariance is assembled from the ten published $4\times
4$ SIDIS submatrices in the HEPData species order. One signed vector
$0.08 A_1$ defines a nuisance shared by all 40 points; the residual published
systematic variance is diagonal.

No statistical covariance was released for either 2017 record. Statistical
errors therefore remain diagonal. One record-wide vector
$0.8\,\sigma_\mathrm{syst}$ is correlated across every species and cell, and
$0.6\,\sigma_\mathrm{syst}$ is added on the diagonal. Monte Carlo statistical
variance is added on the diagonal in all goodness-of-fit calculations.

## Runner contract

Schema version 5 adds `unpolarized_sidis` and per-PDF-axis `active` flags. A
missing flag retains the active behavior of older descriptors. Polarized
SIDIS uses `PP`, `PM`, `MP`, and `MM`; unpolarized SIDIS accepts only `00`.
The 2017 descriptors vary the unpolarized NNPDF40 axis and reject any
noncentral `--polarized-pdf-members` selection. The established polarized
matrix element is retained in explicit zero-polarization mode, so the signed
NLO workflow is unchanged while polarized PDFs are disabled.

Inspect central and paper plans without launching events:

```bash
python3 scripts/run_phenomenology_campaign.py prepare \
  --measurement COMPASS_2009_I820721 \
  --tag compass2009_plan --dry-run

python3 scripts/run_phenomenology_campaign.py prepare \
  --measurement COMPASS_2017_I1444985 \
  --tag compass2017_pion_hadron_paper \
  --profile paper --dry-run
```

Run one structural smoke for each analysis after loading the validated
`herwig/pol` environment:

```bash
for analysis in \
  COMPASS_2009_I820721 \
  COMPASS_2017_I1444985 \
  COMPASS_2017_I1483098
do
  python3 scripts/run_phenomenology_campaign.py full \
    --measurement "$analysis" \
    --tag "${analysis}_smoke" \
    --smoke --jobs 4
done
```

The runner combines shards inside POSNLO and NEGNLO first and adds the two
normalized signed contributions before constructing an asymmetry or ratio.
Generic raw-accumulator `yodamerge` arithmetic is not a valid replacement.
Generated cards, `.run` files, logs, YODA, summaries, and plots stay below the
ignored `campaigns/phenomenology/` tree and are not repository products.
