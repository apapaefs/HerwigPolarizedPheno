# STAR 510 GeV projected-spin production controller

This controller prepares and runs the interim STAR 510 GeV comparison with
the full polarization vector radially projected onto the unit Bloch ball at
construction time. The prescription is explicitly provisional: the strict
negative-ISR guard remains enabled, and no event or branching is vetoed.

The nominal campaign contains 500 million LO events for each of the four
beam-helicity combinations. It is split into 500 shards per helicity (2,000
shards total) and uses the verified 4 GeV generator cut. Inclusive bins 1--4
remain in summaries as `diagnostic_only`; inclusive bins 5--14, whose first
analysis edge is exactly 13.1 GeV, enter the primary plots, pulls, and
59-point covariance fit. All four dijet selections are unchanged.

The completed 3/4/5 GeV, 50-million-event cut scans are reused read-only.
Their common campaign commit is recorded separately from the current checker,
measurement-policy, and projection-runtime hashes. The controller never
prepares, postprocesses, or mutates those scans.

The active `herwig/pol` installation is not modified. `run.sh` prepends the
isolated projected ThePEG build below, and the runtime verifier checks the
effective dynamic links for both the `Herwig` executable and `Herwig.so` core:

```text
/home/apapaefs/Projects/Herwig/production/star510-bloch-ad0c5486/prefix
```

## User commands

On Odysseus:

```bash
CONTROL=/home/apapaefs/Projects/HerwigPolarizedPheno/campaigns/control/star510-pt13p1-bloch-20260822

"$CONTROL/run.sh" verify
"$CONTROL/run.sh" check-star-cut
"$CONTROL/run.sh" status

screen -DmS star510-pt13p1-bloch -L -Logfile "$CONTROL/star510.log" \
  "$CONTROL/run.sh" star510
```

Recover only failed shards, retaining successful nonempty outputs:

```bash
"$CONTROL/run.sh" recover star510
```

After all 2,000 shards are successful:

```bash
"$CONTROL/run.sh" postprocess
"$CONTROL/run.sh" plot
"$CONTROL/run.sh" package
```

`package` exports exactly five vector PDFs and five 600-dpi PNGs into a new
commit-tagged directory below `figures/`, with a JSON provenance inventory.
Neither HERMES/COMPASS nor PHENIX is run or repackaged by this controller.
