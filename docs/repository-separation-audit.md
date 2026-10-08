# Code and paper separation audit — 8 October 2026

Retain the current `HerwigPolarizedPheno` analysis implementation and the
current `Phenomenological-Investigations-of-Polarized-Collisions-in-Herwig-7`
manuscript. No newer analysis implementation exists only in the paper
repository. The paper's unique recent changes are manuscript text and the
four HERMES v5 figure PDFs with their provenance.

This decision supersedes the July 2026 policy of pushing an identical whole
repository to the code and Overleaf remotes. Code and paper now have separate
owners and histories can advance independently without rewriting either one.

## Compared snapshots

| Role | Snapshot |
| --- | --- |
| Code, committed local main | `0f7ed13b74a8f1cfba405891cd56ebcd5130b114` |
| Paper, current main | `c1d55d7f482bd561149083a43ef44eeb7fe2e121` |
| Last common ancestor | `78d1be3962c3e058409634af28dc9cd52dc9b77d` |

The current paper checkout is the sibling directory
`/Users/apapaefs/Projects/Phenomenological-Investigations-of-Polarized-Collisions-in-Herwig-7`.
The same-named checkout inside `HerwigPol/` is an older analysis-audit branch
at `f6af80a` and is not the manuscript baseline.

The code's `41c173a` and paper's `f62d515` both add the pp analysis catalogue.
They are distinct commits with identical stable patch ID
`9bc47d07f21ae76c74df646a9a60b4fa0f6f8320`, not a common ancestor.
The paper's Overleaf import `1193988` deleted most scripts, descriptors,
reference data and workflow documentation and dropped six controller
executable bits. Preserve the complete code versions and executable modes.
Do not interpret those deletions as improved analysis implementations.

The paper's changes after `f62d515` affect only `main.tex`, `main_old.tex`,
`FIGURE_SOURCES.md`, and `figures/hermes-v5-20261007/`.

The paper remote's non-main branches were also checked. The archived original
import (`e37d50d`), experimental-analysis audit (`f6af80a`) and hard-process-spin
branch (`f31b80f`) are ancestors of its current main. The orphan
`codex/phenomenology-paper-restart` branch (`f1072d5`) contains one unique
commit with only seven manuscript, bibliography, style and build files.
No non-main paper branch contains additional unmerged analysis work.

## File-by-file preservation

The audit classified 296 tracked analysis, run-card, controller and
collaborator-material paths in the paper snapshot. The split removes 292
code paths; the four `SIDIS_ByFrank/` author-reference files stay in the
paper repository and are removed from the code repository:

| Classification | Paths | Choice |
| --- | ---: | --- |
| Identical code contents at the same path | 257 | Keep the code copy and its executable modes. |
| Identical collaborator author-reference contents | 4 | Keep the paper copy under `SIDIS_ByFrank/`. |
| HERMES files with newer versions in code | 5 | Keep the corrected code versions. |
| Old proton-only HERMES cards | 30 | Keep their existing explicit `LEGACY` copies in code. |

By directory, 99 of 103 shared analysis files are identical, 142 of 143
shared cards are identical, and all 15 campaign-controller payloads and all
four `SIDIS_ByFrank/` files are identical. The remaining identical path is
`SIDIS_IMPLEMENTATION_ASSESSMENT.md`. Six controller payloads have executable
mode in code and non-executable mode in the paper import.

Every one of the 30 old HERMES cards exactly matches its existing code copy
after replacing `HERMES_2007_I726689` with
`HERMES_2007_I726689_LEGACY` throughout the path and text. They contain no
otherwise missing change and should not be restored under the active names.
The active cards now distinguish proton and neutron components.

The machine-readable [preservation manifest](repository-separation-preservation.json)
records all 296 compared paths, Git blob IDs, modes, normalized legacy-card
comparisons and final owning repositories. The 261 identical-content count
in that snapshot comprises 257 code files and four paper author-reference
files. Generated or untracked campaign products are outside this classification.

## Scientific choices

Keep all current COMPASS implementations. Their shared files are identical
between the repositories, including the September 5 corrections from
`84689f4` and clean-checkout reference handling from `8b86c70`. These retain
species-dependent whole-z-bin nu restrictions, explicit nominal-beam support
masks, signed binned slope fits and covariance-aware azimuthal fits. There is
no conflicting paper-side COMPASS improvement to choose.

The five differing HERMES paths are:

- `analyses/rivet/dis/HERMES_2007_I726689.cc`
- `analyses/rivet/dis/HERMES_2007_I726689.info`
- `analyses/rivet/dis/HERMES_2007_I726689.plot`
- `analyses/rivet/dis/HERMES_2007_I726689.yoda.gz`
- `cards/experimental/HERMES_2007_I726689/HERMES_2007_I726689-Common.in`

Retain these code versions together with their matching descriptors, scripts,
reference snapshots and tests. They carry the validated September 30 and
October 6 evolution: separate proton/neutron components and simulated
deuteron combination; explicit normalized POSNLO-minus-NEGNLO arithmetic;
direct published Born longitudinal cells; distinct model-assisted GD11
projections; full R1990 three-fit average for the A1 proxy; the rectangular
aperture with correlated ring controls; and loose generation cuts that cover
the physical massive-target cuts. Reverting to the paper's old HERMES files
would lose these corrections.

The direct Born asymmetry is ordinary LL/UU. Keep inverse-D A1, direct Born
cells, GD11 projections and acceptance controls distinct. This repository
separation introduces no new estimator or physics correction and does not
justify relabeling existing campaign signatures or regenerating events.

The original local code checkout also contains uncommitted October 7 ratio
renderers and hooks, results browsers and open-charm reference work. Preserve
that work; the paper repository supplies no competing versions. Distinguish
these working files from the committed snapshots above when publishing the
separation.

## Scope and validation

Analysis implementations, generator run cards, reference inputs, scripts,
tests, controller sources and workflow documentation belong to the code
repository. Manuscript source, bibliography, styles, paper catalogues,
selected figure assets and figure provenance belong to the paper repository.
The README and Makefile in each repository now match its role. The distinct
pre-split code manuscript and bibliography are preserved verbatim in the
paper repository at `author-notes/code-draft-20261006/`; the current
manuscript and bibliography remain unchanged.

The final preservation review confirmed all 103 analysis files, 349 cards,
26 configuration files, 363 reference-input files and 15 controller files
match the original local code contents. Every one of the 110 removed code
manuscript/figure/author-reference paths has a paper destination. All but
`FIGURE_SOURCES.md` match exactly; that file retains its historical inventory
and adds the newer HERMES selection and historical-context clarification.
No remaining paper source or figure input changed, and no broken relative
Markdown links were found in the code documentation.

Neither compared snapshot vendors the Herwig or ThePEG generator source.
Those remain external dependencies in the existing generator repositories;
do not describe a cards/workflow separation as a migration of those external
source trees. Any explicit later generator-source consolidation needs its
own version and runtime audit.

Paper compilation, undefined citation/reference checks, all four HERMES
figure provenance hashes and the 296 code-path dispositions are recorded
separately in the separation delivery. Generated campaigns remain untouched.

## Code validation — 8 October 2026

The isolated code checkout preserves 28 local source/reference files,
including the October 7 HERMES ratio pipeline, results browsers and
production-disabled open-charm references. Two intentional repairs affect
the browser and its test: the browser now validates and hashes additional
Born reference snapshots and their pinned raw sources in the runner's
order, and the regression test rejects corruption of either input. The
signature mismatch was reproduced in the original checkout before repair;
it was not introduced by the separation. The original checkout was kept
unchanged during this audit and test phase.

The full Python suite ran 444 tests: 438 passed and six were skipped because
native YODA/Rivet dependencies were unavailable. All 26 focused browser and
open-charm tests passed. Herwig, Rivet and rivet-build were unavailable on
the local PATH, so no native smoke test was run. No campaigns were launched.

The published split includes the existing October 6 v5 commit `0f7ed13`
and the preserved local additions. Initial separation commits `44033a5`
(code) and `95fb220` (paper) were pushed as ordinary fast-forward updates and
verified on both GitHub `main` refs. The two original local `main` branches
were then fast-forwarded to the published split.

Before updating the original code checkout, its incorporated local source
was saved in Git stash `f9a2958361cd3e0c4efde886d93122a5c964a700` (label:
`Preserved local source before reviewed code-paper split 2026-10-08`). This
is a recovery copy; its changes are already incorporated and should not be
blindly reapplied. All 137 unrelated untracked files and the paper editor
lock were verified unchanged. The original code checkout's `overleaf-github`
fetch remote is retained for historical reads, with its push URL disabled.
