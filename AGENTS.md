# HerwigPolarizedPheno working contract

This repository owns Rivet analyses, generator cards, reference inputs,
campaign controllers, postprocessing, plotting tools and tests. The separate
Phenomenological-Investigations-of-Polarized-Collisions-in-Herwig-7 repository
owns the manuscript and selected figures. Do not restore whole-repository
Overleaf mirroring or push this repository's branches to the paper remote.
Herwig/ThePEG implementation source remains an external herwigdispol dependency.

Read README.md, the relevant workflow documentation and
`docs/repository-separation-audit.md` before changing workflows. Existing
campaigns and uncommitted work must be preserved. Do not commit generated
events, logs, ordinary YODA, plots or installed/build products; pinned
compressed reference data and source campaign controllers are intentional.

Keep independent physical-helicity generation distinct from correlated
helicity weights. Combine POSNLO/NEGNLO at normalized-bin level and retain
stored-order signs. Do not infer helicity-dependent shower/Sudakov reweighting
from target weights. Validate the active executable and installed libraries
before running smoke campaigns; do not interrupt or rebuild active campaigns.
Run focused tests for changes and a bounded end-to-end smoke when the runtime
is available. Figure export must retain provenance and published uncertainties.
