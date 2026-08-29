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

## Pilot

After synchronizing the repository on Odysseus, run:

```bash
module purge
module load herwig/pol
cd /home/apapaefs/Projects/HerwigPolarizedPheno

python3 scripts/run_mc_poljetshapes_campaign.py full \
  --tag mc_poljetshapes_pilot_500k_20260829_v1 \
  --families nominal,shower_spin_off \
  --lo-events 500000 \
  --shards 50 --jobs 100 \
  --seed-base 8207000 \
  --plot-comparisons --include-diagnostics
```

This is 500,000 events for each helicity and shower family, split into 400
shards total. The postprocessor evaluates the 100M, 250M, and 500M tiers. It
selects the smallest tier projecting at least 250,000 effective entries for
each baseline $\Delta\psi_{12}$ jet/helicity/family sample and an independent
spin-on/off `A2UU` and `A2LL` error no larger than 0.002. Inter-jet and
four-jet limitations are reported but cannot increase the bounded
recommendation above 500M.

The ready-to-copy command is stored in
`postprocess/statistics-projection.json`. It always uses 500,000 events per
shard, 100 jobs, seed base `8307000`, and a fresh tag of the form
`mc_poljetshapes_spin_<tier>m_20260829_v1`. These projections assess
statistical precision; they do not guarantee a physical separation if the
true spin-on/off difference is zero.
