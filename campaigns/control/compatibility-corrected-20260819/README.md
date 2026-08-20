# Compatibility-corrected production controller

This controller locks the corrected COMPASS, HERMES, and STAR 510 GeV
comparison campaigns to one Git commit and one Odysseus runtime. It prepares
7,600 production shards and 1,200 STAR generator-cut shards. HERMES low-
(Q^2) is validation-only and withheld from the result package; PHENIX is not
prepared or run.

`runtime-lock.json` pins the Herwig executable, `HerwigDefaults.rpo`,
`HwMEDIS`, `HwMEHadron`, `HwShower`, `FixedTargetLuminosity`, Rivet, and the
complete NNPDF4.0/NNPDFpol2.0 file inventories. Each campaign manifest also
pins its compiled Rivet plugin, generated cards, `.run` files, source commit,
and runtime artifacts.

## Validation and preparation boundary

The maintainer validation action runs the complete Python suite, compiles and
registers every Rivet analysis, and completes commit-tagged smoke campaigns
for all production measurements plus diagnostic-only HERMES legacy. The
production preparation action is separate: it creates fresh cards, plugins,
manifests, and `.run` files, then verifies that every manifest is `prepared`,
every shard is `planned`, every production YODA output is absent, and no
matching `Herwig run` process exists.

The preparation action invokes only the campaign runner's `prepare` stage. It
does not invoke `campaign`, `full`, postprocessing, or plotting. Production is
therefore an explicit user action after handoff.

## User handoff

On Odysseus:

```bash
CONTROL=/home/apapaefs/Projects/HerwigPolarizedPheno/campaigns/control/compatibility-corrected-20260819

"$CONTROL/run.sh" verify
"$CONTROL/run.sh" commands
"$CONTROL/run.sh" status
```

Launch fixed-target and SIDIS production:

```bash
screen -DmS pheno-fixed -L -Logfile "$CONTROL/fixed.log" \
  "$CONTROL/run.sh" fixed
```

Launch and assess the three STAR generator-cut samples:

```bash
screen -DmS star510-cutscan -L -Logfile "$CONTROL/star-cutscan.log" \
  "$CONTROL/run.sh" star-cut-scan

"$CONTROL/run.sh" check-star-cut
```

Only after that gate passes, launch nominal STAR 510 production:

```bash
screen -DmS star510-500m -L -Logfile "$CONTROL/star510.log" \
  "$CONTROL/run.sh" star510
```

Failed shards are recovered with fresh seeds while successful nonempty
outputs are retained:

```bash
"$CONTROL/run.sh" recover fixed
"$CONTROL/run.sh" recover star-cut-scan
"$CONTROL/run.sh" recover star510
```

After every production manifest is complete:

```bash
"$CONTROL/run.sh" postprocess
"$CONTROL/run.sh" plot
"$CONTROL/run.sh" package
```

The package action selects exactly four inclusive fixed-target, 24 HERMES
SIDIS, and five STAR 510 overlays. It retains each vector PDF and renders a
600-dpi PNG sibling in a new commit-tagged directory without overwriting an
existing package.

## STAR gate

The first two finite bins of each primary STAR observable gate the 3-versus-4
GeV comparison. `SigmaUU` must agree within the larger of 2% and three
combined Monte Carlo standard errors; `A_LL` must agree within the larger of
(5\times10^{-4}) and three combined Monte Carlo standard errors. The 5 GeV
comparison is recorded as a non-gating stress test. The report refuses
campaigns with incomplete shards, reused initial seeds, mismatched cuts, or
different source commits. Nominal STAR production also rechecks the report's
manifest and summary hashes before launch.
