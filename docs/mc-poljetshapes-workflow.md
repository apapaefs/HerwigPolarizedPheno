# MC_POLJETSHAPES shower-spin workflow

`MC_POLJETSHAPES` is an internal 510 GeV particle-level measurement for
comparing the nominal polarized parton shower with the same polarized hard
events showered with spin correlations disabled. Each generated event is sent
to both `MC_POLDIJETS:LEVEL=HADRON` and `MC_POLJETSHAPES`; the original
`MC_POLDIJETS` descriptor, cards, and completed campaigns are unchanged.

The final state contains all stable visible particles. Hadronization, hadron
decays, and remnants are enabled, MPI is disabled, and jets are anti-$k_T$
with $R=0.5$, $|\eta|<1.5$, $p_{T,1}>5$ GeV, and $p_{T,2}>4$ GeV. No shower
history, parton flavour, or hard-process truth label enters the jet-shape
analysis. The `nominal` and `shower_spin_off` cards differ only in
`ShowerHandler:SpinCorrelations`; their seeds are independent.

## Particle-level observables

Jet constituents are reclustered with Cambridge/Aachen and E-scheme
recombination. At a declustering,
$z=p_{T,\mathrm{soft}}/(p_{T,\mathrm{hard}}+p_{T,\mathrm{soft}})$ and
$k_T=p_{T,\mathrm{soft}}\Delta R$. The primary Lund plane follows the hard
branch. The secondary plane follows the hard branch of the selected primary
splitting's soft branch. Signed plane angles are accumulated recursively, so
selected splittings need not be adjacent in the declustering sequence.

The three $\Delta\psi_{12}$ working points are:

- loose: $z_1,z_2>0.1$ and $k_{T,1},k_{T,2}>0.5$ GeV;
- symmetric secondary: $z_1>0.1$, $z_2>0.3$, and both $k_T>0.5$ GeV;
- perturbative: $z_1,z_2>0.1$ and both $k_T>1$ GeV.

Two additional hard-sharing working points require
$0.25<z_1<0.40$ and $z_2>0.35$, with both selected splittings above either
$k_T>1$ GeV or $k_T>2$ GeV. The signed $\Delta\psi_{12}$ and squeezed EEEC
are recorded after these selections. This deliberately concentrates on hard,
reasonably symmetric branchings where a helicity-dependent azimuthal
correlation is less likely to be diluted by strongly ordered soft radiation.

Resolved radiation is analysed with

\[
\Delta\phi_{3|1}=\angle[(\mathrm{beam},j_1),(j_1,j_3)],\qquad
\Delta\psi_{34}=\angle[(j_1,j_3),(j_1,j_4)],
\]

signed about the leading-jet axis. The postprocessor reports cosine and sine
moments in fixed bins of $p_{T,3}/p_{T,1}$ and $p_{T,4}/p_{T,1}$. The sine
moments are null tests.

The hard-sharing $k_T>1$ GeV sample is also divided by two experimentally
constructible gluon-enrichment proxies: $N_\mathrm{const}\geq8$ versus its
complement, and $p_T^D<0.45$ versus its complement, where
$p_T^D=\sqrt{\sum_i p_{T,i}^2}/\sum_i p_{T,i}$. Both use only stable visible
jet constituents; no truth-flavour information enters either category.

The analysis also records the beam--jet-plane to primary-splitting angle,
the two-jet $\Delta\psi_{11'}$, an exact deterministic squeezed EEEC triplet
sum with $p_T$-fraction weights, $Q_{2,\beta}$ and sine-null quadrupoles for
$\beta=1,2$, and the four-jet Bengtsson--Zerwas angle after boosting the four
jets to their centre-of-mass frame and ordering them by energy there.
Inclusive $R_{32}$, $R_{43}$, third-jet veto, and cumulative
$p_{T,3}/p_{T,1}$ and $p_{T,4}/p_{T,1}$ scans use the fixed thresholds in the
checked analysis snapshot.

Every raw distribution yields the four direct helicity cross sections,
$\sigma_{UU}$, $\Delta\sigma_{LL}$, and $A_{LL}$. Angular distributions also
produce unit-normalized `PP`, `PM`, `MP`, `MM`, and `UU` shapes; cosine moments
`A2UU` and `A2LL`; and sine moments `B2UU` and `B2LL`. The binned moment uses
the same midpoint discrete Fourier construction as the plotted histogram and
propagates numerator--denominator covariance within each helicity sample.

Postprocessing writes the red spin-on and blue spin-off prediction YODAs,
an independent-sample on-minus-off YODA, JSON and CSV moment/ranking summaries,
cutflow and accepted-entry summaries, a bounded statistics projection, and a
complete HTML gallery. The comparison is a generator-description test, not a
comparison with experimental points.

For the third- and fourth-jet spectra, $R_{32}$, $R_{43}$, veto and cumulative
tails, and the resolved-angle moments versus jet-$p_T$ fraction, postprocessing
also writes `postprocess/shard-block-covariance.json` and a flat CSV form. The
matrix is a delete-one-common-shard jackknife: one equal-size shard with the
same ordinal is removed from every helicity sample in a shower family. Since
spin on and spin off are independently generated, their covariance matrices
are added for the difference. The focused HTML page summarizes the resulting
covariance-aware $\chi^2$ values and links to the full matrices.

The plot stage also constructs presentation-only spin-on/spin-off ratios for
positive cross sections, normalized shapes, and rate observables. Both the
numerator and denominator Monte Carlo errors are propagated because the two
shower families use independent event samples. Ratios are deliberately not
formed for signed `DeltaSigmaLL`, `ALL`, or cosine/sine moments, where a zero
crossing would make the ratio misleading; those use the existing on-minus-off
plots instead.

In addition to the complete gallery, plotting writes
`plots/html/focus/index.html`. This focused page is linked prominently from the
main `index.html` and places each recommended red/blue overlay beside either
its spin-on/spin-off ratio or its on-minus-off difference. Its fixed hierarchy
covers resolved-radiation moments, hard-sharing declusterings, conditioned
EEEC shapes, particle-level gluon-enriched proxy categories, inclusive
Delta-psi moments and shapes, sine/null tests, the UU, Delta-sigma LL, and
A_LL third- and fourth-jet transverse-momentum spectra, radiation-rate scans,
and inclusive controls. The resolved-jet spectra are kept together so a
change in the amount of additional radiation can be distinguished from a
longitudinal-spin asymmetry. The selection is fixed in code and is not chosen
from the noisy pilot sensitivity ranking.

## Pilot

After synchronizing the repository on Odysseus, run:

```bash
module purge
module load herwig/pol
cd /home/apapaefs/Projects/HerwigPolarizedPheno

python3 scripts/run_mc_poljetshapes_campaign.py full \
  --tag mc_poljetshapes_conditional_pilot_500k_20260830_v1 \
  --families nominal,shower_spin_off \
  --lo-events 500000 \
  --shards 50 --jobs 100 \
  --seed-base 8357000 \
  --plot-comparisons --include-diagnostics
```

This is 500,000 events for each helicity and shower family, split into 400
shards total. The postprocessor evaluates the 250M, 500M, and 1B tiers. It
selects the smallest tier projecting at least 250,000 effective entries for
each hard-sharing $k_T>1$ GeV $\Delta\psi_{12}$ jet/helicity/family sample and
an independent spin-on/off `A2UU` and `A2LL` error no larger than 0.002. The
$k_T>2$ GeV and four-jet limitations are reported but cannot increase the
bounded recommendation above 1B.

The ready-to-copy command is stored in
`postprocess/statistics-projection.json`. It always uses 500,000 events per
shard, 100 jobs, seed base `8407000`, and a fresh tag of the form
`mc_poljetshapes_conditional_spin_<tier>_20260830_v1`. These projections assess
statistical precision; they do not guarantee a physical separation if the
true spin-on/off difference is zero.

For the bounded high-statistics tier, run:

```bash
module purge
module load herwig/pol
cd /home/apapaefs/Projects/HerwigPolarizedPheno

python3 scripts/run_mc_poljetshapes_campaign.py full \
  --tag mc_poljetshapes_conditional_spin_1b_20260830_v1 \
  --families nominal,shower_spin_off \
  --lo-events 1000000000 \
  --shards 2000 --jobs 100 \
  --seed-base 8407000 \
  --plot-comparisons --include-diagnostics
```

This is 1B generated events per helicity and shower family (8B total), with
500,000 events per shard. It uses a new immutable measurement signature and
does not modify `mc_poljetshapes_spin_500m_20260829_v1`.

## Replot the completed conditional production campaign

The focused page and propagated ratios require only the already-postprocessed
family predictions and `postprocess/summary.json`. They do not require new
events or a repeated postprocessing pass. After the running production has
finished and the Odysseus checkout has been updated to the plotting commit,
regenerate the gallery with:

```bash
module purge
module load herwig/pol
cd /home/apapaefs/Projects/HerwigPolarizedPheno

python3 scripts/run_mc_poljetshapes_campaign.py plot \
  --tag mc_poljetshapes_conditional_spin_1b_20260830_v1 \
  --plot-comparisons --include-diagnostics
```

This replaces only the generated `plots/html/` subtree and updates the plot
section of the campaign manifest. Shards, raw YODAs, normalized predictions,
and statistical summaries remain untouched. The blue family is a
spin-averaged shower on polarized hard events, not an unpolarized-beam sample.
The older 500M campaign predates these raw observables and must not be reused
under the new measurement signature.
