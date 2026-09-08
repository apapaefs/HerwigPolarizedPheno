# MC_POLJETSHAPES_LHE: hard-process spin information loss

This is a separate internal measurement. It uses the existing 510 GeV polarized
LO pp dijet cards and both exact `MC_POLJETSHAPES` / `MC_POLDIJETS` analysis
instances. All particle-level observables and selections are unchanged.

The two independent production families are:

- `nominal`: `HardProcessSpin Yes`, `SpinCorrelations Yes` (red).
- `lhe_like_shower`: `HardProcessSpin No`, `SpinCorrelations Yes` (blue).

The polarized hard matrix element and PDFs are unchanged. Only the new switch
differs between the families' physics settings. Seeds and output identifiers
differ intentionally. Do not substitute the legacy `shower_spin_off` family:
that test leaves polarized backward-ISR conditioning in place.

The new mode is LHE-like, not generically LHE-equivalent. Release requires the
separate generator validation, including actual ordinary LHE replay. Use only
the validated versioned installation, not an arbitrary `herwig/pol` module.
Preflight tests both values of the new switch and verifies the loaded shower
and hadronic-matrix-element libraries. Manifests and postprocessing summaries
record their fingerprints, analysis instances and both switches.

## Bounded pilot and production

After the runtime has passed closure, a pilot is:

```bash
python3 scripts/run_mc_poljetshapes_lhe_campaign.py full \
  --tag mc_poljetshapes_lhe_pilot_500k_20260908_v1 \
  --families nominal,lhe_like_shower \
  --lo-events 500000 --shards 50 --jobs 100 \
  --seed-base 9107000 --plot-comparisons --include-diagnostics
```

This generates 500,000 events per helicity/family: 4 million in 400 shards.
Never reuse an existing immutable tag with changed source, configuration or
runtime. The dedicated runner retains preparation, execution/recovery,
postprocessing, plotting and status subcommands from the shared engine.

`postprocess/statistics-projection.json` reports the smallest passing tier among
100M, 250M and 500M events per helicity/family. The requirements are at least
250,000 effective entries in **each** loose baseline intra-jet angle and each
helicity/family, and independent full/LHE-like A2UU and A2LL difference errors
at most 0.002. Rare selections do not force an unbounded run. A
`no_bounded_tier_satisfies_all_criteria` result means no tier meets the target;
the capped fallback command is not a claim that it does. Production commands
use 500,000 events per shard, 100 jobs and seed base 9207000. They are not
launched automatically.

## Replotting a completed tag

```bash
python3 scripts/run_mc_poljetshapes_lhe_campaign.py postprocess \
  --tag YOUR_COMPLETED_TAG --include-diagnostics
python3 scripts/run_mc_poljetshapes_lhe_campaign.py plot \
  --tag YOUR_COMPLETED_TAG --plot-comparisons --include-diagnostics
```

The complete HTML gallery includes the focus panel, jet-3 / jet-4 spectra and
radiation rates. Ratios are labelled full spin / LHE-like, and moment differences
full spin minus LHE-like. Existing internal filenames containing `on-off` are
retained for compatibility; their numerator/denominator metadata is explicit.
Independent-sample error propagation is retained. Ratios with poorly determined
denominators and asymmetries/moments with insufficient effective support are
masked, not displayed as apparently precise pulls.

Existing descriptors, manifests, completed campaigns and STAR archives must not
be modified to use this new comparison.
