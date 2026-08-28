# SIDIS tranche campaign controller

This tracked controller prepares six SIDIS analyses: the five first/second-
tranche measurements plus `COMPASS_2020_I1788430`,
launches only an explicitly selected stage, recovers failed shards without
replacing successful products, and gates central statistics on a checksum-
pinned pilot assessment. It never launches a campaign from `verify`,
`commands`, `dry-run`, `status`, or `assess`.

The initial pilot contains 64.9 million generated events. The starting central
floor is 618.2 million events in 11,800 shards. At 100 Odysseus workers the
generation-only extrapolation is about 7 hours; allow 12--30 hours for the
high-dimensional and rare-high-z Rivet analyses, filesystem traffic,
postprocessing, and recovery. The six paper profiles would contain about 65.0 billion events and
must be staged per analysis only after its central gate.

The COMPASS 2018 HEPData v1 submission contains 4,664 numerical cells, even
though the earlier assessment stated 4,918. The implementation and campaign
gate use all 4,664 released cells and record the discrepancy; no cells are
invented.

## Odysseus workflow

While STAR production is active, use only `verify`, `commands`, and dry runs in
the detached clean side worktree. Do not rebuild the active runtime or launch
SIDIS concurrently.

```bash
CONTROL=/home/apapaefs/Projects/HerwigPolarizedPheno-sidis-tranches-20260828/campaigns/control/sidis-tranches-20260828

"$CONTROL/run.sh" verify
"$CONTROL/run.sh" commands
"$CONTROL/run.sh" dry-run pilot all
"$CONTROL/run.sh" prepare pilot all

screen -DmS sidis-pilot -L -Logfile "$CONTROL/pilot.log" \
  "$CONTROL/run.sh" launch pilot all

"$CONTROL/run.sh" status pilot all
"$CONTROL/run.sh" recover pilot all
"$CONTROL/run.sh" postprocess pilot all
"$CONTROL/run.sh" assess pilot all
```

If masked bins prevent extrapolation, the assessor writes an immutable
doubled-statistics `pilot2` plan. Run it with the corresponding
`prepare pilot2 all` and `launch pilot2 all` commands. After the gate passes:

```bash
"$CONTROL/run.sh" prepare central all

screen -DmS sidis-central -L -Logfile "$CONTROL/central.log" \
  "$CONTROL/run.sh" launch central all

"$CONTROL/run.sh" status central all
"$CONTROL/run.sh" recover central all
"$CONTROL/run.sh" postprocess central all
"$CONTROL/run.sh" plot central all
"$CONTROL/run.sh" package central all

"$CONTROL/run.sh" dry-run paper all
```

`commands` prints the equivalent raw runner invocations, including events,
shards, jobs, seed bases, profiles, PDF selectors, and scales. A selector may
be `all` or one analysis ID. Generated plans, manifests, YODA, logs, plots, and
packages live below ignored campaign/runtime directories and are not source
products.

## Statistical gate

Every primary bin must be finite. At least 90% must satisfy
`sigma_MC <= 0.5 sigma_exp,total`, and every finite bin must satisfy
`sigma_MC <= sigma_exp,total`. Recommendations use
`N_required = N_pilot (sigma_MC/(f sigma_exp))^2`, keep the 10:1 POS/NEG
ratio, round upward to complete shards, and never fall below the central
floors. A masked pilot cannot be extrapolated and instead produces a doubled
pilot plan. A bin with nonzero experimental data but zero prediction and zero
Monte Carlo variance is likewise treated as unresolved, since a zero-entry
sample cannot support a `1/sqrt(N)` extrapolation. Overlapping HERMES binnings
are assessed independently and are never combined into one goodness of fit.
