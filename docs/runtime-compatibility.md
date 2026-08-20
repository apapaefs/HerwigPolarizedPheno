# Runtime compatibility

The repository is self-contained for the phenomenology workflow, reference
data, Rivet analyses, postprocessing, tests, and paper. Event generation still
requires an installed Herwig Polarized environment.

## Required tools

- Python 3.10 or newer;
- the Python bindings shipped with the active YODA/Rivet installation;
- NumPy, with PrettyTable optional for richer live status output;
- `Herwig`, `rivet`, `rivet-build`, `rivet-mkhtml`, `rivet-config`, `lhapdf`,
  and `lhapdf-config` on `PATH`;
- `NNPDF40_nlo_pch_as_01180` and `NNPDFpol20_nlo_as_01180`; the `paper`
  profile needs all members 0--100 of both ensembles;
- `latexmk` and TeX only when building the manuscript.

`requirements.txt` covers the ordinary pip-installable Python dependencies.
YODA should come from the same installation as Rivet rather than from an
unrelated Python package.

## Herwig Polarized interfaces

The active Herwig installation must provide:

- the polarized DIS matrix-element library (`HwMEDIS`) and physical-helicity
  POWHEG DIS setup;
- the polarized hadron matrix elements (`HwMEHadron`) used for STAR and
  PHENIX samples;
- the active shower library (`HwShower`) used by the generated event samples;
- `FixedTargetLuminosity` and the installed snippets `EPCollider.in`,
  `PolarizedEP.in`, `FixedTarget.in`, `PPCollider.in`, and `PolarizedPP.in`;
- `DISBase::MinimumScale` for the low-scale HERMES projection;
- `ScalePreFactor` on the weak and QCD hard processes;
- incoming-density contraction and hard-spin-vertex support in `MEQCD2to2`;
- the corrected polarized QTilde backward-ISR normalization and parent-density
  arguments for physical-helicity shower samples.

The imported workflow is the audited `apapaefs/HerwigPol` tree at
`822911203d431daa2a026fc697aabf5eeeb060b2`; the published descendant
`04ada4408aa71e0674d0e0c494933badaf6eb688` adds the polarized QTilde ISR
correction. Before interpreting new production, confirm that the executable
and dynamically loaded `HwMEDIS`, `HwMEHadron`, and shower libraries belong to
the intended installation. Rebuilds in a source tree alone do not establish
that `Herwig run` loads them.

After a persistent-interface change, regenerate `.run` files. Never replace
an installed library while a campaign using that installation is running.

## Preflight

From the repository root:

```bash
python3 scripts/run_phenomenology_campaign.py list
make test
make check-rivet
```

The campaign runner performs additional executable, library, PDF, card, and
manifest checks during `prepare` or `full`. Every prepared manifest records
SHA-256 identities for `Herwig`, `HerwigDefaults.rpo`, `HwMEDIS`,
`HwMEHadron`, `HwShower`, `FixedTargetLuminosity`, the Rivet executable, the
campaign-specific Rivet plugin, every generated card and `.run` file, and
every file in both LHAPDF sets. A path without a matching file hash is not a
runtime identity.

The compatibility-corrected Odysseus controller adds a checked-in runtime
lock and refuses validation or production when the checkout is dirty,
detached, on the wrong branch, or different from the corresponding GitHub
branch. Its production mode additionally requires the canonical path
`/home/apapaefs/Projects/HerwigPolarizedPheno`. The validation report is
commit- and runtime-lock-specific, so switching commits invalidates it even
when the executable paths are unchanged.

The Makefile first uses the compiler reported by `rivet-config` when that
executable still exists, then falls back to versioned GNU C++ installations.
Set `RIVET_CXX=/path/to/c++` explicitly if the automatic choice is unsuitable.
