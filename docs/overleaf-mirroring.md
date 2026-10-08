# Separate code and Overleaf paper repositories

The July 2026 whole-repository mirroring policy is superseded by the
8 October 2026 split. Keep the existing Git histories; do not force-push,
rebase the published paper lineage or replace it with an unrelated history.

- [HerwigPolarizedPheno](https://github.com/apapaefs/HerwigPolarizedPheno)
  owns analysis source, Herwig cards, reference inputs, campaign scripts,
  plotting tools and tests.
- [Phenomenological-Investigations-of-Polarized-Collisions-in-Herwig-7](https://github.com/apapaefs/Phenomenological-Investigations-of-Polarized-Collisions-in-Herwig-7)
  owns the manuscript, bibliography, styles, selected figures, provenance,
  author notes and paper build. This remains the Overleaf GitHub repository.
- [herwigdispol](https://github.com/apapaefs/herwigdispol) supplies the
  external Herwig/ThePEG implementation; this split does not relocate it.

Publish each repository to its own `origin`. Their `main` commits are
expected to differ. Never run the former `git push overleaf-github main:main`
from the code checkout. An existing `overleaf-github` remote may be retained
for historical reads, but must not be a publication destination for code.

Generate and validate results in the code checkout. Copy selected final
figure assets and their input/source/runtime provenance into the paper
checkout, then run `make paper` there and check the log for undefined
citations and references. Preserve the current manuscript while exporting
figures; do not synchronize whole directory trees. Overleaf can continue to
synchronize with its existing paper repository after the paper-only change.

See [the separation audit](repository-separation-audit.md) for the reviewed
source snapshots and preservation checks.
