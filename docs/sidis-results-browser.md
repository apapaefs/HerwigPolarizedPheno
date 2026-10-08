# Portable SIDIS results browser

`scripts/build_sidis_results_browser.py` builds a self-contained browsing page,
one detailed page per analysis, and optionally a gzip-compressed tarball. It
uses only the Python standard library (Python 3.10 or newer), reads existing
campaigns, and never starts Herwig, rebuilds plots, changes a manifest, or
contacts a remote machine. Herwig/Rivet modules are not needed for packaging.

The default selection is exactly the 13 `*_central_20260903_v1` campaigns from
the all-SIDIS launch. It never guesses the most recent campaign. Each page
includes the experiment, observable/data type, numerical authority, publication
and data links, implemented cuts, beam/target model, binning, estimator,
limitations, PDF definitions, campaign configuration and source identity.
Generation success is explicitly distinguished from physics validation.

## Build on Odysseus

If only the exporter is new locally, copy this one self-contained script from
the Mac (the metadata and campaign inputs already belong to the source checkout):

```bash
scp /Users/apapaefs/Projects/HerwigPolarizedPheno/scripts/build_sidis_results_browser.py \
  odysseus:/home/apapaefs/Projects/HerwigPolarizedPheno/scripts/
```

Once the exporter script is present in the matching source checkout:

```bash
python3 /home/apapaefs/Projects/HerwigPolarizedPheno/scripts/build_sidis_results_browser.py \
  --output /home/apapaefs/Projects/HerwigPolarizedPheno/campaigns/exports/sidis-results-20260904 \
  --tar
```

This creates the directory above and the adjacent
`sidis-results-20260904.tar.gz`. The command prints the bundle/archive sizes,
plot counts, and archive SHA-256. An existing directory or archive is never
overwritten; use a fresh output name for an updated export.

For a smaller image-only viewing bundle, add `--plot-format png` (no PDFs).
The default includes both PNG and PDF. Use `--dry-run` to check selection,
source compatibility, completeness and plot/CSV payload size without writing.

## Copy to the Mac and open

Run on the Mac, with a new destination directory:

```bash
mkdir -p /Users/apapaefs/Projects/HerwigPolarizedPheno/campaigns/exports
scp odysseus:/home/apapaefs/Projects/HerwigPolarizedPheno/campaigns/exports/sidis-results-20260904.tar.gz \
  /Users/apapaefs/Projects/HerwigPolarizedPheno/campaigns/exports/
tar -xzf /Users/apapaefs/Projects/HerwigPolarizedPheno/campaigns/exports/sidis-results-20260904.tar.gz \
  -C /Users/apapaefs/Projects/HerwigPolarizedPheno/campaigns/exports
open /Users/apapaefs/Projects/HerwigPolarizedPheno/campaigns/exports/sidis-results-20260904/index.html
```

All gallery links and images work directly from disk. No server or internet
connection is required except when following an external publication/data link.

## What is included

- Fresh static HTML with inline CSS and relative links, with no CDN/fonts/JS dependencies.
- Rendered PNG/PDF assets, including available diagnostic plots, preserving
  their relative subdirectory names.
- `central.csv` when that campaign provides it.
- Compact analysis/campaign metadata, numerical-source checksum manifests,
  and `bundle.json` with sizes and SHA-256 checksums for every payload file.

The exporter does **not** include campaign logs, event/shard YODAs, HepMC,
ROOT, `.run` files, grids, cards, work/build directories, libraries, raw paper
or experimental archives, original plot scripts/data helpers, or large dense
postprocessing summaries/covariance matrices. It generates clean gallery pages
rather than retaining Rivet download links to excluded helper files.

This is a viewing package, not a resumable or fully reproducible campaign
backup. Preserve the original campaign, analyzed YODA, covariance, and pinned
reference repository. No original file is removed during export.

## Completeness and source safety

By default, export fails before writing if a selected manifest, successful
shard set, postprocessed prediction/summary, recorded plot index, or PNG/PDF
plot pair is missing. This tests file availability, not numerical quality or
Monte Carlo precision. It does not regenerate missing products.

`--allow-incomplete` creates an explicitly labelled partial/metadata-only
browser, with no invented plots or completion claims. Campaign identity and
generation-signature mismatches remain hard errors even with that option.
Use the source checkout matching the campaign; do not relabel old plots with
new cuts. Input symlinks, path escapes, and output inside the campaign input
tree are rejected. The source-signature calculation is regression-tested
against the campaign runner.

For another campaign or a subset, select analyses and override exact tags:

```bash
python3 /home/apapaefs/Projects/HerwigPolarizedPheno/scripts/build_sidis_results_browser.py \
  --measurement COMPASS_2014_I1278730 \
  --campaign COMPASS_2014_I1278730=compass2014_azimuth_central_20260903_v1 \
  --output /home/apapaefs/Projects/HerwigPolarizedPheno/campaigns/exports/compass2014-view-v2 \
  --tar
```

`--measurement` and `--campaign` can be repeated. `--repo` and
`--campaign-root` allow an explicit source checkout and existing campaign
location, including a mounted copy. Outputs under `campaigns/exports/` remain
ignored by Git.

## Data-authority distinctions

HEPData is authoritative for the 2009, 2017, 2018, 2025 and 2026 COMPASS
measurements. COMPASS 2010, 2013 and 2014 and the high-z 2020 ratios use exact
paper-table transcriptions, not digitized or inferred curves. The associated
HEPData identifiers for 2013/2014 do not change that implemented source.
HERMES 2019 uses the APS supplement. HERMES 2013 multiplicities use the
archived five-binning release as 3D/covariance authority, with HEPData only as
a projection cross-check. HERMES 2013 azimuthal moments remain generator-only:
the exact binning is implemented, but numerical reference values and covariance
are not vendored.
