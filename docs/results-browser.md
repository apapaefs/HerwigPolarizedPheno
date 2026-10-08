# Latest-results browser

`scripts/build_results_browser.py` scans the experimental and phenomenology
campaign directories and creates one portable, searchable plot browser for
all registered analyses. It uses Python 3.10 or newer, the standard library,
and matplotlib to render reference-only data panels. Herwig, Rivet and a web
server are not needed to build or view it.

Run this on the host that holds the campaign outputs. For the current production
results, use the primary Odysseus checkout:

```bash
cd /home/apapaefs/Projects/HerwigPolarizedPheno
python3 scripts/build_results_browser.py --tar
```

Each invocation creates a new `campaigns/exports/results-TIMESTAMP/` directory
and, with `--tar`, a matching `.tar.gz`. The command prints the exact index and
archive paths, selected campaign tags, plot counts and archive checksum. Copy
the archive to your computer, extract it, and open `index.html`. All browsing,
search and plot downloads work offline. Re-run the command to scan for newer
results; an existing browser is a snapshot and does not update itself.

For results already on a Mac, run the same Python command from the local
checkout. The script only finds files on the host where it runs; it does not
automatically contact Odysseus or retrieve event files.

To browse on demand without copying a large archive, start a preview from your
Mac in a separate terminal (use an available port):

```bash
ssh -L 127.0.0.1:8770:127.0.0.1:8770 odysseus \
  'python3 -m http.server 8770 --bind 127.0.0.1 --directory /home/apapaefs/Projects/HerwigPolarizedPheno/campaigns/exports'
```

Open `http://127.0.0.1:8770/` and choose a results directory. Leave the SSH
terminal open while browsing; Ctrl-C closes the preview. The exported archive
remains available for a complete offline copy.

## What the browser includes

- An overview of every registered analysis, including analyses with no current
  results, plus filters for experiment and result availability.
- A searchable, paginated gallery of every exported PNG/PDF plot. Search can
  match the analysis, campaign, plot filename, title or axis labels.
- Analysis information: cuts, beams, target model, PDFs, physics families,
  reference authority, publication links and interpretation limits.
- Displayed and latest campaign identities, shard progress, source signatures,
  run configuration, history and process observations where available.
- Recorded masked-bin counts and fit diagnostics. Large diagnostic arrays are
  explicitly abbreviated; the browser never performs a fit or combines YODA.
- Numerical `central.csv` or legacy `summary.csv` tables, compact analysis metadata, reference source
  manifests and a SHA-256 inventory of the viewing files.
- Reference-only entries from `config/reference/`, currently the 45 COMPASS
  open-charm points in nine panels. These export the published CSV/JSON and
  are labelled **Reference data only**, with simulation requirements visible
  in Analysis info. They are not counted as ready simulation results.

Campaign logs, event files, raw YODA, run files, libraries and dense covariance
arrays are not copied. This export is a viewing bundle, not a campaign backup.

## How automatic selection works

1. Read the analysis registries and manifests in
   `campaigns/phenomenology/<ANALYSIS>/<TAG>/` and
   `campaigns/experimental/<ANALYSIS>/<TAG>/`. Explicitly recorded smoke tests
   are excluded by default.
2. Check the campaign's measurement signature against the current source,
   reference and card definitions. Experimental comparison signatures are
   checked when applicable. Incompatible histories remain visible, but their
   plots are withheld so that old results cannot inherit newer physics labels.
3. Order campaigns by `created_at`, not by tag spelling or a later replot time.
   If creation time is missing, the manifest modification time is used and
   labelled as a fallback in the metadata.
4. Select the newest compatible result whose recorded jobs all succeeded and
   whose prediction, summary, plot index and PNG/PDF pairs exist. A newer
   unfinished attempt is shown separately, with an explicit fallback notice.
   If no complete result exists, available compatible partial plots are shown
   with a warning.

A presentation-only plot-metadata refresh is accepted only if its recorded
historical `.plot` bytes can be read from the local Git history, their checksum
matches, and substituting those bytes reproduces the original generation
signature. The script never fetches Git history or trusts a refresh label alone.

“Ready” describes recorded generation and available products; it does not
establish experimental agreement or statistical convergence. On Linux, the
script observes this user's local campaign controllers and Herwig workers.
An incomplete manifest with no observed controller or worker is labelled
stopped, even if its last recorded state says `running`. On hosts without
`/proc`, process activity is explicitly marked unchecked. These observations
describe the build time, not ongoing monitoring.

The script only reads campaign inputs. It does not resume jobs, rebuild
analyses, regenerate campaign plots, merge histograms or overwrite an existing
export. For reference-only registrations it verifies the pinned data and
renders fresh data-only panels inside the new export. An explicit campaign
selection is not permitted for these entries.

## Options

Preview selections without copying files:

```bash
python3 scripts/build_results_browser.py --dry-run
```

Export a subset, choose a new output directory, or reduce the bundle to PNGs:

```bash
python3 scripts/build_results_browser.py \
  --measurement COMPASS_2013_I1236358 \
  --measurement COMPASS_2014_I1278730 \
  --output campaigns/exports/compass-fit-review-NEW-TAG \
  --plot-format png
```

Repeat `--campaign-root PATH` to scan nondefault campaign roots. A root must
contain `<ANALYSIS>/<TAG>/manifest.json`; specify both roots if replacing the
defaults. Use `--campaign ANALYSIS=TAG` for a deliberate compatible campaign
selection, and `--include-smoke` to consider smoke tests. Explicit campaign
selection still enforces source compatibility. An output directory must be
separate from the input trees and must not already exist.

The older [`build_sidis_results_browser.py`](../scripts/build_sidis_results_browser.py)
and its [13-analysis export instructions](sidis-results-browser.md) remain
available for their pinned SIDIS viewing snapshot.

## Validation

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover \
  -s scripts/tests -p 'test*results_browser.py'
```

Tests cover signature parity with both campaign runners, creation-time ordering,
unfinished and incompatible campaigns, explicit selections, plot-metadata
refresh verification, process observations, portable links and checksums,
safe exports, and preservation of input manifests and existing bundles.
