# HERMES proton and deuteron preparation, 30 September 2026

The later [Born A_parallel extension](hermes-born-apar-20260930.md) has a
fresh prepared Odysseus tag, `hermes_pd_born_3m_20260930_v4`, with the same
event budget. Use its source area and launch command to obtain the Born
cells and both integrated projections. The v3 preparation below records
the preceding A1/deuteron and normalized-order audit.

The updated `HERMES_2007_I726689` analysis accepts separate fixed proton and
neutron components and derives both published Table 14 target views from one
campaign. The existing proton paths remain available; the deuteron reference
and prediction use `d14-x01-y02`. The `Q2 > 4 GeV2` control is retained for
each target.

## Source and campaign locations

Canonical local analysis source:

```text
/Users/apapaefs/Projects/HerwigPolarizedPheno
```

The isolated development checkout is
`/Users/apapaefs/Projects/HerwigPol/.worktrees/hermes-deuteron-20260930`
on `codex/hermes-deuteron-20260930`. Scoped changes are applied to the
canonical checkout; unrelated existing edits are retained.

The isolated Odysseus source and runtime area is:

```text
/home/apapaefs/Projects/Herwig/validation/hermes-deuteron-20260930/pheno
```

The new nominal production tag is `hermes_pd_3m_20260930_v3`. Its generated
inputs, `.run` files, manifest, logs and outputs reside beneath that isolated
root in:

```text
campaigns/experimental/HERMES_2007_I726689/hermes_pd_3m_20260930_v3/
```

It matches the event budget and shard layout of the completed Odysseus
proton-only `compat-central-3m-20260819-v1` campaign. A new seed range keeps
the preparation separate from that historical sample. Old manifests and
YODA outputs must retain their original proton-only definitions; updating
their labels cannot supply missing neutron events or the new target outputs.

## Production matrix

| Setting | Value |
| --- | --- |
| Target components | `P` (`p+`) and `N` (`n0`) |
| Physical helicities per component | `PP`, `PM`, `MP`, `MM` |
| Contributions per helicity | `POSNLO`, `NEGNLO` |
| Logical jobs | 16 |
| POSNLO events per component and helicity | 3,000,000 |
| NEGNLO events per component and helicity | 300,000 |
| Shards per logical job | 100 |
| Events per POSNLO / NEGNLO shard | 30,000 / 3,000 |
| Total shard jobs | 1,600 |
| Total requested events | 26,400,000 |
| Concurrent workers | 100 |
| Seed base | `560726689` |
| Profile | `central`, nominal full spin |

The opt-in comparison matrix would contain 60 logical jobs, with both
components generated for all five prediction families. It is separate from
this nominal preparation and requires a new tag. The descriptor's generic
100,000/10,000 defaults remain available; the explicit options below select
the user's matched Odysseus budget.

## Preparation command

On Odysseus, load the same module used by the existing HERMES production,
then run from the isolated source root:

```bash
cd /home/apapaefs/Projects/Herwig/validation/hermes-deuteron-20260930/pheno
source /etc/profile.d/modules.sh
module purge
module load herwig/pol
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1

python3 scripts/run_phenomenology_campaign.py prepare \
  --measurement HERMES_2007_I726689 \
  --tag hermes_pd_3m_20260930_v3 \
  --profile central \
  --posnlo-events 3000000 \
  --negnlo-events 300000 \
  --shards 100 \
  --jobs 100 \
  --seed-base 560726689
```

The primary wrapper delegates this fixed-target measurement to
`scripts/run_experimental_campaign.py`. Preparation preflights the resolved
Herwig executable, installed `HwMEDIS`, fixed-target library, Rivet toolchain
and both PDF sets, and compiles the updated Rivet plugin. It regenerates the
16 logical `.run` files before generation. Explicit normalized order
coefficients restore the negative
contribution sign in postprocessing; rebuilding a separate source tree
alone is insufficient evidence for the runtime being used.

## Launching the prepared production

Use the ignored launcher below so Odysseus's existing thermal guard recognizes
the established `run_validation_campaign.py full` command form. It is a
symlink to the primary phenomenology runner; it executes the same updated
fixed-target engine. Keep all prepared generation options unchanged:

```bash
cd /home/apapaefs/Projects/Herwig/validation/hermes-deuteron-20260930/pheno
source /etc/profile.d/modules.sh
module purge
module load herwig/pol
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1

python3 campaigns/launchers/hermes-deuteron/run_validation_campaign.py full \
  --measurement HERMES_2007_I726689 \
  --tag hermes_pd_3m_20260930_v3 \
  --profile central \
  --posnlo-events 3000000 \
  --negnlo-events 300000 \
  --shards 100 \
  --jobs 100 \
  --seed-base 560726689 \
  --progress-interval 5 \
  --max-listed 32
```

`full` reuses the matching prepared configuration, generates pending shards,
postprocesses normalized bins and renders the proton and deuteron overlays.
During execution the atomic tracker is available at
`campaigns/experimental/HERMES_2007_I726689/hermes_pd_3m_20260930_v3/monitor/status.json`
and its `status.txt` sibling. Completed nonempty shard YODAs are recoverable;
do not discard them or change a manifest to resume.

## Verified preparation status

The final `hermes_pd_3m_20260930_v3` manifest is **prepared**: all 1,600
shards are planned, all 16 logical `.run` files are nonempty, and no
production YODA outputs exist. Production events have not been launched.
The earlier v1/v2 preparations remain unchanged. The paper audit found that
the active generator stores NEGNLO as a positive magnitude. v3 records the
corrected `POSNLO-NEGNLO` normalized-bin combination and a new seed range.
v3 has the corrected A1/order-sign definition; v1/v2 are obsolete. The later
[Born A_parallel extension](hermes-born-apar-20260930.md) changes the analysis
signature and requires a fresh v4 preparation to obtain its cell and
x-integrated Q2 histograms. v3 is not a preparation for those new outputs.

Validation completed on Odysseus:

- All **247 regression tests passed**, including the independent depolarization and signed-order tests.
- `hermes_pd_smoke_20260930_v3` completed **60/60 jobs**, requesting 100
  events per job, and postprocessed all five prediction families.
- The primary proton/deuteron references contain 15 published points each;
  both overlays and the four deuteron plot titles were checked. Empty
  low-statistics denominators remain explicitly masked.
- The active `herwig/pol` executable and installed generator-library hashes
  match the completed proton-only HERMES production. The updated Rivet
  plugin is private to the isolated campaign area.
- The thermalguard-compatible launcher passed a preparation dry-run with
  the exact final event counts, shard count, worker count and seed base.

These smoke samples verify execution and output construction; they do not
establish precision or agreement with the data. The nominal production
measurement signature is:

```text
0556fc25305904a1391d535be6b6d0b740dc1d968962e93a9005067f5f3e92dd
```

## Target combination

For each target component and helicity, subtract normalized NEGNLO
magnitudes from POSNLO bins first. Form the helicity averages and differences
separately for each
component, then apply the output-specific target coefficients:

```text
proton:   sigma_UU = sigma_UU^p; sigma_LL = sigma_LL^p
deuteron: sigma_UU = (sigma_UU^p + sigma_UU^n)/2
          sigma_LL = 0.925*(sigma_LL^p + sigma_LL^n)/2
```

The active generator stores both order contributions as positive magnitudes.
For each component and helicity, first form `POSNLO-NEGNLO` at normalized-bin
level. Independent statistical variances and covariance-proxy moments retain
both contributions with squared order coefficients.

Apply the same longitudinal coefficients to the `1/D`-weighted numerator
for `A1`; form both asymmetry ratios only after this combination. Independent
component variances use squared coefficients, and the ordinary/weighted
covariance uses their coefficient products. The proton and deuteron outputs
share the proton events; their individual error bars do not provide a joint
cross-target covariance matrix.

The neutron PDFs use ThePEG's `u <-> d` and `ubar <-> dbar` isospin mapping
of `NNPDF40_nlo_pch_as_01180` and `NNPDFpol20_nlo_as_01180`. Separate neutron
LHAPDF sets are unnecessary. The HERMES 27.6 GeV positron beam, cuts,
R1990 depolarization model and neglected `eta*A2` assumption remain the
same for both components.

## Outputs and interpretation

| Summary selection | Target | A1 path |
| --- | --- | --- |
| `Q2GT1` | proton | `/HERMES_2007_I726689/d14-x01-y01` |
| `Q2GT4` | proton | `/HERMES_2007_I726689/A1_Q2GT4` |
| `D_Q2GT1` | deuteron | `/HERMES_2007_I726689/d14-x01-y02` |
| `D_Q2GT4` | deuteron | `/HERMES_2007_I726689/A1d_Q2GT4` |

Each target also has `A_parallel`, normalized `SigmaUU`/`SigmaLL` and
photon-parity controls. Accepted-kinematics diagnostics retain the proton
convention. Summary JSON records each output's resolved target map in
`output_target_combinations`.

The deuteron curve is a per-nucleon free-proton/free-neutron impulse
approximation. The fixed factor `0.925` models the D-state spin reduction;
binding, Fermi motion, off-shell, shadowing and tensor structure functions
are absent. The published deuteron reference already includes its experimental
tensor-asymmetry correction. No extra correction is applied to those data.
The `Q2 > 1 GeV2` output reaches the PDF validity boundary and remains
exploratory; the `Q2 > 4 GeV2` view is the conservative validation control.

The specific paper-level audit is in
[the HERMES validation report](hermes-paper-validation-20260930.md). It
records the omitted `eta*A2` contribution and common-Q2 averaging difference;
these remain limits on a precision comparison with the published A1.

The complete workflow and scientific definitions are documented in
[the HERMES campaign guide](hermes-experimental-campaign.md).
