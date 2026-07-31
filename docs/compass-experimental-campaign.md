# COMPASS inclusive polarized-DIS campaigns

This workflow compares HerwigPol particle-level fixed-target simulations to
the final one-dimensional COMPASS inclusive measurements:

- COMPASS_2016_I1357198: \(A_1^p(x)\) from the 2011 200 GeV
  positive-muon data, HEPData record 72819 v1, Table 2.
- COMPASS_2017_I1501480: \(A_1^d(x)\) from the combined 160 GeV
  positive-muon data, HEPData record 78374 v1, Table 1.

The analyses deliberately do not implement the HEPData two-dimensional
tables. Those tables publish mean \(Q^2\) values but not a complete set of
unambiguous \(Q^2\)-bin boundaries suitable for exact event-level binning.
The HEPData \(g_1\) columns are retained in the raw and normalized snapshots,
but \(g_1\) is not predicted: it requires the external unpolarized \(F_2\)
convention used by COMPASS, and several high-\(x\) deuteron HEPData values
differ from the final journal table.

## Reference-data contract

The official HEPData table-display JSON is vendored byte-for-byte below
data/experimental/<measurement>/raw-hepdata-table.json. The corresponding
reference.json is the normalized, inspectable input used to build the Rivet
reference YODA.

The fetch-data command downloads the same official endpoint and refuses to
proceed unless all of the following agree:

- source SHA256;
- record version, table name, and table DOI;
- column headers, reaction, and beam-momentum qualifier;
- row count, bin edges, means, \(A_1\), and statistical/systematic errors;
- the vendored raw payload and normalized snapshot.

The proton source has a published typo in row 5: the bin is encoded as
0.008-0.001, while its mean is 0.009 and the next bin starts at 0.010. The
normalizer applies the explicitly registered 0.001 -> 0.010 correction only
after verifying the raw value. The official payload remains unchanged.

Refresh and validate both tables with:

~~~bash
python3 scripts/run_experimental_campaign.py fetch-data \
  --measurement COMPASS_2016_I1357198

python3 scripts/run_experimental_campaign.py fetch-data \
  --measurement COMPASS_2017_I1501480
~~~

## Rivet physics

COMPASSInclusiveDIS.hh supplies the shared implementation. It identifies the
highest-energy prompt final-state positive muon, excluding hadron-decay
leptons, and reconstructs

\[
Q^2=-q^2,\qquad x=\frac{Q^2}{2P\cdot q},\qquad
y=\frac{P\cdot q}{P\cdot k},\qquad W^2=(P+q)^2
\]

directly in the fixed-target laboratory frame.

The proton selection uses 200 GeV, \(1<Q^2\le190\,\mathrm{GeV}^2\),
\(0.1<y<0.9\), \(W^2>12\,\mathrm{GeV}^2\), and
\(0.0025<x<0.7\). The deuteron components use 160 GeV,
\(1<Q^2\le100\,\mathrm{GeV}^2\), \(0.1<y<0.9\),
\(W>4\,\mathrm{GeV}\), and \(0.004<x<0.7\). Both analyses also fill a
nested \(Q^2>4\,\mathrm{GeV}^2\) validation view.

The depolarization factor is the exact finite-muon-mass expression printed by
COMPASS. \(R(x,Q^2)\) is the E143 R1998 average of the published
\(R_a,R_b,R_c\) fits. The conversion uses
\(A_{LL}=D(A_1+\eta A_2)\) with \(\eta A_2=0\), following the COMPASS
extraction. The analyses fill ordinary, \(1/D\), and \(1/\sqrt D\)
covariance-proxy histograms. POSNLO and NEGNLO are added only after each
sample and shard has been normalized.

## Deuteron approximation

Herwig generates separate physical proton and neutron targets. The campaign
constructs a per-nucleon deuteron impulse approximation:

\[
\sigma_{UU}^d=\frac{\sigma_{UU}^p+\sigma_{UU}^n}{2},\qquad
\sigma_{LL}^d=0.925\,\frac{\sigma_{LL}^p+\sigma_{LL}^n}{2}.
\]

Here \(0.925=1-1.5\omega_D\) with \(\omega_D=0.05\). The same coefficients
are applied to covariance terms before taking ratios. This is not a nuclear
event generator: binding, Fermi motion, off-shell effects, low-\(x\)
shadowing, and the spin-1 tensor function \(b_1\) are absent. ThePEG applies
its standard neutron isospin map to both the unpolarized and polarized
proton PDF sets.

## Herwig samples

The nominal setup uses NNPDF40 NLO and NNPDFpol2.0 NLO, pure photon
exchange, native DIS-window POWHEG generation, the QCD shower,
hadronization, decays, and remnant handling. MPI and QED shower radiation
are disabled. The nominal sample retains the POWHEG real-emission spin
vertex and shower spin correlations.

The proton matrix contains eight logical jobs:
four physical helicities times POSNLO/NEGNLO. The deuteron matrix contains
sixteen because each helicity/order pair is generated for proton and neutron
targets.

The --comparisons option enables the same labels as the HERMES workflow:

- NLO+PS (polarized; full spin)
- NLO+PS (unpolarized beams)
- LO+PS (polarized)
- NLO+PS (polarized; Born spin only)
- NLO+PS (polarized; shower spin off)

The expanded matrices contain 30 proton and 60 deuteron logical jobs.
Comparison families are opt-in and immutable in the campaign manifest.

NNPDFpol2.0 includes COMPASS information. These comparisons are therefore
primarily implementation and generator-closure tests, not predictions
independent of the measured data.

## Running

Load a clean polarized environment first:

~~~bash
source /opt/homebrew/opt/modules/init/bash
module purge
module load herwig/pol
~~~

Inspect a resolved smoke matrix without writing:

~~~bash
python3 scripts/run_experimental_campaign.py prepare \
  --measurement COMPASS_2017_I1501480 \
  --tag compass-dry \
  --smoke \
  --comparisons \
  --dry-run
~~~

Run nominal smoke campaigns:

~~~bash
python3 scripts/run_experimental_campaign.py full \
  --measurement COMPASS_2016_I1357198 \
  --tag compass-proton-smoke \
  --smoke \
  --jobs 4

python3 scripts/run_experimental_campaign.py full \
  --measurement COMPASS_2017_I1501480 \
  --tag compass-deuteron-smoke \
  --smoke \
  --jobs 4
~~~

Add --comparisons to create the expanded matrices. Production event counts,
sharding, seeds, recovery, and live tracking use the same options as the
HERMES experimental runner. All generated cards, .run files, logs, YODA,
manifests, summaries, and HTML are written below
campaigns/experimental/<measurement>/<tag>/.

The \(Q^2>1\,\mathrm{GeV}^2\) results touch the selected PDFs' lower scale
boundary and must be labelled exploratory. Use the \(Q^2>4\,\mathrm{GeV}^2\)
objects as the conservative generator-validation view.
