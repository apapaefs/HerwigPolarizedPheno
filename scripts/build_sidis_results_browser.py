#!/usr/bin/env python3
"""Export a read-only, offline SIDIS plot browser and an optional small tarball.

Python standard library only: no Herwig/Rivet environment, event regeneration,
network access, or campaign mutation. Defaults select the 13 September-2026
central campaigns explicitly, never an arbitrary newest directory.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import hashlib
import html
import json
import os
import re
import shutil
import sys
import tarfile
import tempfile
from pathlib import Path
from urllib.parse import quote, urlsplit


ROOT = Path(__file__).resolve().parents[1]
TOKEN = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*\Z")
# Presentation metadata stays separate from immutable generation descriptors.
# Paper values below are transcribed from tables, not inferred/digitized curves.
CATALOG = {
    "COMPASS_2009_I820721": (
        "compass2009_sidis", "HEPData", "Polarized identified-hadron A1ᵈ",
        "HEPData v1 is the numerical authority; the frozen paper/HEPData discrepancy audit is retained in the source repository."),
    "COMPASS_2010_I862410": (
        "compass2010_sidis", "Paper tables", "Polarized identified-hadron A1ᵖ",
        "Exact paper-TeX table transcription with a frozen extraction audit; not values inferred from a plot."),
    "COMPASS_2013_I1236358": (
        "compass2013_ptslopes", "Paper tables", "Unpolarized low-pT inverse slopes",
        "368 tabulated inverse slopes from paper Tables 1–3. Herwig spectra are fitted separately; tabulated data are not fitted or digitized from figures. Only the published fit errors are available."),
    "COMPASS_2014_I1278730": (
        "compass2014_azimuth", "Paper tables", "Unpolarized cos(φ) and cos(2φ) amplitudes",
        "480 cosine amplitudes from paper Tables 2–12. The implemented numerical source is the paper, despite the associated HEPData record. The beam-spin sin(φ) amplitude is excluded."),
    "COMPASS_2017_I1444985": (
        "compass2017_pion_hadron", "HEPData", "Charged-pion and charged-hadron multiplicities",
        "Exact sparse HEPData v1 cells; published correction columns are provenance, not corrections applied to Herwig events."),
    "COMPASS_2017_I1483098": (
        "compass2017_kaon", "HEPData", "Charged-kaon multiplicities",
        "Exact sparse HEPData v1 cells. Isoscalar hadron and DIS yields are combined before their ratio and z-bin-width normalization."),
    "COMPASS_2018_I1624692": (
        "compass2018_pt2", "HEPData", "pT²-dependent charged-hadron multiplicities",
        "All 4,664 numerical cells actually released in HEPData v1 are retained; missing cells are not invented."),
    "COMPASS_2020_I1788430": (
        "compass2020_charge_ratio", "Paper tables", "High-z antiproton/proton and K−/K+ ratios",
        "67 paper-table rows; no official HEPData submission was found in the pinned 2026-08-28 audit. The three overlapping tables are separate comparison groups."),
    "COMPASS_2025_I2840545": (
        "compass2025_multiplicity", "HEPData", "Hydrogen charged-hadron multiplicities",
        "Complete HEPData v1 submission; final radiative- and diffractive-vector-meson-corrected data compared to uncorrected generator-level hadrons."),
    "COMPASS_2026_I3096394": (
        "compass2026_multiplicity", "HEPData", "Deuteron charged-hadron multiplicities",
        "Complete HEPData v1 submission; the target is modeled by an isoscalar proton/neutron combination, not a nuclear event generator."),
    "HERMES_2013_I1111237": (
        "hermes2013_azimuth", "No numerical data", "Unpolarized azimuthal-moment diagnostic",
        "Paper-defined 5×5×6×6 binning (900 cells per target/species). The 21,600 numerical moments and covariance are not vendored. No experimental values, pseudo-data, or goodness of fit are fabricated."),
    "HERMES_2013_I1208547": (
        "hermes2013_multiplicity", "Official archive + HEPData cross-checks", "Hydrogen/deuterium identified-hadron multiplicities",
        "The archived HERMES five-binning release is the authority for 6,592 usable 3D cells and statistical covariance. The 64 HEPData tables are projection cross-checks, not the 3D source. Overlapping binnings are not independent datasets."),
    "HERMES_2019_I1698889": (
        "hermes2019_sidis", "Official APS supplement", "Longitudinal identified-hadron A∥",
        "Official APS numerical tables and covariance. Overlapping 1D/2D/3D projections are kept separate; spin-asymmetry Fourier amplitudes are not unpolarized cosine moments."),
}

LABELS = {
    "q2_min_gev2": "Q² minimum [GeV²]", "w2_min_gev2": "W² minimum [GeV²]",
    "w_min_gev": "W minimum [GeV]", "pt_gev": "pT [GeV]",
    "pt2_gev2": "pT² [GeV²]", "fit_pt_gev": "Inverse-slope fit pT [GeV]",
    "x": "Bjorken x", "y": "Inelasticity y", "z": "Hadron z",
    "xf_min": "xF minimum", "xF_min": "xF minimum",
    "hadron_lab_angle_mrad": "Hadron laboratory angle [mrad]",
    "virtual_photon_angle_max_rad": "Virtual-photon angle maximum [rad]",
    "hadron_momentum_gev": "Hadron momentum [GeV]",
}
CSS = """
:root{color-scheme:light;--ink:#172b3e;--muted:#536477;--line:#d8e0e8;--blue:#145b90;--wash:#f2f6fa}
*{box-sizing:border-box}body{margin:0;background:white;color:var(--ink);font:16px/1.6 system-ui,-apple-system,sans-serif}
header{border-top:6px solid var(--blue);border-bottom:1px solid var(--line);padding:24px max(24px,calc((100vw - 1200px)/2))}
main{max-width:1248px;margin:auto;padding:24px}h1{font-size:2rem;line-height:1.2;margin:.3em 0}h2{font-size:1.35rem;line-height:1.4}
h3{font-size:1.05rem}a{color:var(--blue);text-underline-offset:3px}a:focus-visible,summary:focus-visible{outline:3px solid #dd8700;outline-offset:4px}
p{margin:.5em 0 1em}.muted{color:var(--muted)}.meta{font-size:.875rem}.mono,code{font-family:ui-monospace,monospace;overflow-wrap:anywhere}
.notice{background:#fff7df;border-left:4px solid #a56800;padding:12px 18px;margin:20px 0}.ok{color:#166044}.warning{color:#854b00}
.tablewrap{overflow-x:auto}table{border-collapse:collapse;width:100%;text-align:left}th,td{padding:14px 12px;border-bottom:1px solid var(--line);vertical-align:top}
th{font-size:.875rem;background:var(--wash)}td small{display:block;font-size:.875rem;color:var(--muted)}.badge{display:inline-block;border:1px solid var(--line);padding:2px 9px;border-radius:4px;font-size:.875rem}
.facts{display:grid;grid-template-columns:1fr 1fr;gap:16px 32px}.panel{border-top:3px solid var(--blue);padding-top:4px;min-width:0}
dl{display:grid;grid-template-columns:minmax(130px,1fr) 2fr;gap:8px 16px}dt{color:var(--muted)}dd{margin:0;overflow-wrap:anywhere}
details{border-bottom:1px solid var(--line);padding:12px 0}summary{cursor:pointer;font-weight:600}pre{white-space:pre-wrap;overflow-wrap:anywhere;font-size:.875rem;background:var(--wash);padding:14px}
.gallery{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,360px),1fr));gap:24px}.plot{margin:0;border:1px solid var(--line);padding:12px;min-width:0}
.plot img{display:block;width:100%;height:auto}.plot figcaption{font-size:.875rem;overflow-wrap:anywhere}.links{display:flex;gap:14px;flex-wrap:wrap;margin:12px 0}
footer{border-top:1px solid var(--line);padding:20px 24px;font-size:.875rem;color:var(--muted)}ul{padding-left:22px}li{margin:5px 0}
@media(max-width:700px){.facts{grid-template-columns:1fr}header,main{padding:18px}h1{font-size:1.65rem}dl{grid-template-columns:1fr}dd{margin-bottom:8px}}
"""


class ExportError(RuntimeError):
    pass


def digest(path: Path) -> str:
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def contained(path: Path, root: Path) -> bool:
    try:
        path.resolve().relative_to(root.resolve())
        return True
    except ValueError:
        return False


def regular_file(path: Path, root: Path) -> bool:
    """Never follow a plot/metadata link into logs or outside the source tree."""
    if not contained(path, root):
        raise ExportError(f"Path escapes source tree: {path}")
    for part in (path, *path.parents):
        if part == root:
            break
        if part.is_symlink():
            raise ExportError(f"Refusing source symlink: {part}")
    return path.is_file() and path.stat().st_size > 0


def read_json(path: Path) -> dict:
    try:
        result = json.loads(path.read_text(encoding="utf-8"))
    except (ValueError, OSError) as exc:
        raise ExportError(f"Cannot read JSON {path}: {exc}") from exc
    if not isinstance(result, dict):
        raise ExportError(f"Expected a JSON object: {path}")
    return result


def measurement_signature(repo: Path, descriptor: dict) -> str:
    """The runner's byte-level contract, without importing its physics runtime.

    Regression-tested against run_phenomenology_campaign._signature. Do not
    include presentation catalog text in the event-generation signature.
    """
    payload = {k: v for k, v in descriptor.items() if not k.startswith("_")}
    result = hashlib.sha256(json.dumps(payload, sort_keys=True,
                                       separators=(",", ":"), ensure_ascii=True).encode())
    analysis = descriptor["analysis"]
    paths = []
    for spec in [analysis, *analysis.get("companions", [])]:
        paths.extend(repo / spec[k] for k in ("source", "info", "plot"))
        paths.extend(repo / p for p in spec.get("support_files", []))
    reference = descriptor["reference"]
    paths.append(repo / reference["snapshot"])
    paths.extend(repo / reference[k] for k in ("source_manifest", "raw_snapshot") if reference.get(k))
    paths.extend(sorted((repo / descriptor["cards"]["directory"]).glob("*.in")))
    for path in paths:
        if not regular_file(path, repo):
            raise ExportError(f"Missing signature input: {path}")
        result.update(path.relative_to(repo).as_posix().encode())
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                result.update(block)
    return result.hexdigest()


def metadata(repo: Path, identifier: str) -> tuple[dict, dict]:
    descriptor_path = repo / "config/phenomenology" / f"{identifier}.json"
    if not regular_file(descriptor_path, repo):
        raise ExportError(f"Missing descriptor: {descriptor_path}")
    descriptor = read_json(descriptor_path)
    if descriptor.get("id") != identifier or descriptor.get("process_kind") not in {"polarized_sidis", "unpolarized_sidis"}:
        raise ExportError(f"Not the requested SIDIS descriptor: {descriptor_path}")
    reference = descriptor["reference"]
    snapshot_path = repo / reference["snapshot"]
    manifest_path = repo / reference["source_manifest"]
    for path in (snapshot_path, manifest_path):
        if not regular_file(path, repo):
            raise ExportError(f"Missing metadata: {path}")
    snapshot, sources = read_json(snapshot_path), read_json(manifest_path)
    if snapshot.get("measurement") != identifier:
        raise ExportError(f"Reference identity mismatch: {snapshot_path}")
    info_path = repo / descriptor["analysis"]["info"]
    if not regular_file(info_path, repo):
        raise ExportError(f"Missing analysis info: {info_path}")
    info = info_path.read_text(encoding="utf-8")
    _, authority, data_type, note = CATALOG[identifier]
    links = []
    def add_link(label: str, url: str | None) -> None:
        if url and urlsplit(url).scheme in {"http", "https"} and url not in [v[1] for v in links]:
            links.append([label, url])
    add_link("Registered source", reference.get("source_url"))
    for value in re.findall(r"arXiv:(\d{4}\.\d{4,5}(?:v\d+)?)", info):
        add_link("Paper (arXiv)", "https://arxiv.org/abs/" + value)
    doi = reference.get("publication_doi") or snapshot.get("provenance", {}).get("doi")
    if not doi:
        match = re.search(r'doi\s*=\s*["{]([^"}]+)', info)
        doi = match.group(1) if match else None
    if doi:
        add_link("Publication DOI", "https://doi.org/" + doi)
    if reference.get("record_doi"):
        add_link("HEPData record", "https://doi.org/" + reference["record_doi"])
    if reference.get("hepdata_record"):
        add_link("Associated HEPData record (not the extraction source)",
                 "https://www.hepdata.net/record/" + str(reference["hepdata_record"]))
    add_link("Archived HERMES numerical release", sources.get("full_archive", {}).get("url"))
    add_link("Original numerical-release endpoint", snapshot.get("data_status", {}).get("official_endpoint"))
    datasets = snapshot.get("datasets", {})
    records = list(datasets.values()) if isinstance(datasets, dict) else datasets
    selected_keys = ("beam", "beam_energy_gev", "selection", "binning", "binnings",
                     "target_outputs", "target_model", "deuteron_impulse_approximation",
                     "depolarization", "longitudinal_target_scale", "systematics",
                     "estimator", "fit_policy", "corrections", "interpretation",
                     "excluded_observable", "data_status", "projection_policy")
    return descriptor, {
        "id": identifier, "experiment": identifier.split("_")[0],
        "title": descriptor["title"], "data_type": data_type,
        "authority": authority, "authority_note": note,
        "observable": snapshot.get("observable"),
        "data_available": reference.get("kind") != "internal_observable_definition",
        "reference_entries": sum(len(d.get("points", [])) for d in records),
        "dataset_count": len(records), "links": links,
        "physics": descriptor.get("physics", {}),
        "pdf_ensembles": descriptor.get("pdf_ensembles", {}),
        "reference_metadata": {k: snapshot[k] for k in selected_keys if k in snapshot},
        "datasets": [{k: d[k] for k in ("id", "observable", "species", "published_target", "projection", "density_widths") if k in d}
                     | {"reference_entries": len(d.get("points", []))} for d in records],
        "metadata_sha256": {"descriptor": digest(descriptor_path), "reference_snapshot": digest(snapshot_path),
                            "source_manifest": digest(manifest_path)},
        "source_manifest": sources,
        "rivet_status": re.search(r"^Status:\s*(.+)$", info, re.M).group(1) if re.search(r"^Status:\s*(.+)$", info, re.M) else "unspecified",
    }


def inspect_campaign(repo: Path, campaign_root: Path, identifier: str, tag: str,
                     formats: str) -> tuple[dict, list[tuple[Path, Path]]]:
    descriptor, item = metadata(repo, identifier)
    campaign = campaign_root / identifier / tag
    manifest_path = campaign / "manifest.json"
    issues, files = [], []
    item.update(tag=tag, issues=issues, plots=[], campaign={})
    if not regular_file(manifest_path, campaign_root):
        issues.append("Campaign manifest is missing; no completion claim is made.")
        return item, files
    manifest = read_json(manifest_path)
    if manifest.get("measurement") != identifier or manifest.get("tag") != tag:
        raise ExportError(f"Campaign identity mismatch: {manifest_path}")
    if manifest.get("configuration", {}).get("measurement_signature") != measurement_signature(repo, descriptor):
        raise ExportError(f"{identifier}: campaign/source signature mismatch; use the matching source checkout. No metadata are silently substituted.")
    jobs = manifest.get("jobs", [])
    counts = Counter(job.get("status", "unknown") for job in jobs)
    if not jobs or counts.get("success", 0) != len(jobs):
        issues.append(f"Generation is not fully successful: {counts.get('success', 0)}/{len(jobs)} shards.")
    runtime = manifest.get("runtime", {})
    provenance = runtime.get("provenance", {})
    item["campaign"] = {
        "status": manifest.get("status"), "created_at": manifest.get("created_at"),
        "configuration": manifest.get("configuration", {}),
        "shards": dict(counts), "total_shards": len(jobs),
        "requested_events_in_successful_shards": sum(j.get("events", 0) for j in jobs if j.get("status") == "success"),
        "source_control": provenance.get("source_control", {}),
        "herwig_version": runtime.get("herwig_version"), "rivet_version": runtime.get("rivet_version"),
        "manifest_sha256": digest(manifest_path), "signature_matches": True,
        "plot_created_at": manifest.get("plots", {}).get("created_at"),
    }
    for name in ("prediction.yoda", "summary.json"):
        if not regular_file(campaign / "postprocess" / name, campaign):
            issues.append(f"Postprocessing output missing: {name}.")
    plot_root = campaign / "plots/html"
    index = manifest.get("plots", {}).get("index")
    if not index or not regular_file(campaign / index, plot_root):
        issues.append("A completed plot index is not recorded or is missing.")
    candidates = {}
    if plot_root.exists():
        # Only inspect the plotting subtree, never campaign work/log/shard trees.
        for path in sorted(plot_root.rglob("*")):
            if path.suffix.lower() not in {".png", ".pdf"}:
                continue
            if not regular_file(path, plot_root):
                continue
            relative = path.relative_to(plot_root)
            candidates.setdefault(relative.with_suffix("").as_posix(), {})[path.suffix.lower()[1:]] = path
    for name, pair in candidates.items():
        if set(pair) != {"png", "pdf"}:
            issues.append(f"Incomplete PNG/PDF plot pair: {name}.")
        plot = {"name": name, "files": {}}
        for fmt, path in pair.items():
            if formats != "both" and fmt != formats:
                continue
            destination = Path("plots") / path.relative_to(plot_root)
            files.append((path, destination))
            plot["files"][fmt] = destination.as_posix()
        if plot["files"]:
            item["plots"].append(plot)
    if not item["plots"]:
        issues.append("No rendered plots available.")
    csv_path = campaign / "postprocess/central.csv"
    if regular_file(csv_path, campaign):
        files.append((csv_path, Path("central.csv")))
        item["central_csv"] = "central.csv"
    return item, files


def esc(value: object) -> str:
    return html.escape(str(value), quote=True)


def link(path: str, label: str) -> str:
    encoded = path if urlsplit(path).scheme in {"http", "https"} else quote(path, safe="/")
    return f'<a href="{esc(encoded)}">{esc(label)}</a>'


def pretty(value: object) -> str:
    if isinstance(value, (dict, list)):
        return json.dumps(value, ensure_ascii=False, indent=2)
    return str(value)


def facts(values: dict) -> str:
    return "<dl>" + "".join(
        f"<dt>{esc(LABELS.get(k, k.replace('_', ' ')))}</dt><dd>{esc(pretty(v))}</dd>"
        for k, v in values.items()) + "</dl>"


def page(title: str, body: str, *, back: bool = False) -> str:
    nav = link("../../index.html", "← All SIDIS analyses") if back else "Herwig Polarized · SIDIS"
    return f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{esc(title)}</title><style>{CSS}</style></head><body>
<header><div class="meta">{nav}</div><h1>{esc(title)}</h1></header><main>{body}</main>
<footer>Offline results bundle · No event files or logs · External publication/data links require internet access.</footer>
</body></html>'''


def analysis_page(item: dict) -> str:
    refs, campaign = item["reference_metadata"], item["campaign"]
    warnings = "".join(f"<li>{esc(v)}</li>" for v in item["issues"])
    banner = f'<div class="notice"><strong>Incomplete export</strong><ul>{warnings}</ul></div>' if warnings else '<p class="ok">Generation, postprocessing and rendered plot products are present; the source signature matches.</p>'
    links = [link(url, label) for label, url in item["links"]]
    links.extend([link("metadata.json", "Metadata JSON"), link("source-manifest.json", "Reference-source checksums")])
    if item.get("central_csv"):
        links.append(link(item["central_csv"], "Numerical results (CSV)"))
    figure_groups = {}
    for plot in item["plots"]:
        name, paths = plot["name"], plot["files"]
        group = str(Path(name).parent)
        image = f'<img loading="lazy" src="{esc(quote(paths["png"]))}" alt="{esc(name)}">' if "png" in paths else '<p>PDF-only plot</p>'
        downloads = " · ".join(link(path, fmt.upper()) for fmt, path in paths.items())
        figure_groups.setdefault(group, []).append(f'<figure class="plot">{image}<figcaption>{esc(Path(name).name)}<div>{downloads}</div></figcaption></figure>')
    gallery = "".join(f'<h3>{esc(group)}</h3><div class="gallery">{"".join(figures)}</div>' for group, figures in figure_groups.items())
    technical = {k: v for k, v in refs.items() if k != "selection"}
    data_count = f'{item["reference_entries"]:,} reference entries in {item["dataset_count"]} groups (projections can overlap)' if item["data_available"] else "No numerical experimental reference; generator-only diagnostic"
    body = f'''<p class="mono meta">{esc(item['id'])} · {esc(item['tag'])}</p>
<p>{esc(item['data_type'])} <span class="badge">{esc(item['authority'])}</span></p>{banner}
<div class="facts"><section class="panel"><h2>Measurement and data provenance</h2>
<p>{esc(item['authority_note'])}</p><p>{esc(data_count)}</p>
<p class="meta">Rivet metadata status: {esc(item['rivet_status'])}. Successful generation is not experimental validation.</p>
<div class="links">{' '.join(links)}</div></section>
<section class="panel"><h2>Implemented cuts</h2>{facts(refs.get('selection', {}))}</section></div>
<details><summary>Beam, targets, estimator, binning and limitations</summary>{facts(item['physics'])}
<pre>{esc(pretty(technical))}</pre></details>
<details><summary>Campaign configuration and provenance</summary><pre>{esc(pretty(campaign))}</pre>
<h3>PDF definitions (active axes only are varied)</h3><pre>{esc(pretty(item['pdf_ensembles']))}</pre></details>
<details><summary>Dataset inventory</summary><pre>{esc(pretty(item['datasets']))}</pre></details>
<h2>Plots · {len(item['plots']):,}</h2>{gallery or '<p class="warning">No plots copied for this analysis.</p>'}'''
    return page(item["title"], body, back=True)


def index_page(items: list[dict], created: str) -> str:
    rows = []
    for item in items:
        status = "Incomplete" if item["issues"] else "Ready"
        rows.append(f'''<tr><td>{esc(item['experiment'])}<small>{esc(item['id'])}</small></td>
<td>{link('analyses/' + item['id'] + '/index.html', item['data_type'])}<small>{esc(item['title'])}</small></td>
<td>{esc(item['authority'])}</td><td class="{'warning' if item['issues'] else 'ok'}">{status}<small>{len(item['plots']):,} plots</small></td></tr>''')
    incomplete = sum(bool(item["issues"]) for item in items)
    return page("SIDIS results", f'''<p>{len(items)} analyses · {sum(len(i['plots']) for i in items):,} plots · {incomplete} incomplete · Exported {esc(created)}</p>
<p>Select an analysis for its cuts, beam and target model, data provenance, numerical results and plots.</p>
<div class="tablewrap"><table><thead><tr><th>Experiment / analysis</th><th>Observable and results</th><th>Numerical source</th><th>Export status</th></tr></thead><tbody>{''.join(rows)}</tbody></table></div>
<div class="notice">These are generator-level comparisons, not full detector-level reproductions. The September central campaigns do not supply PDF/scale uncertainty envelopes. Azimuthal diagnostics are not controlled TMD/twist-3 predictions. HERMES 2013 azimuthal data are unavailable in the pinned inputs.</div>
<p>Paper-table transcription is distinct from inference or digitization of a plotted curve. The detailed pages identify the numerical authority even when a separate HEPData record exists.</p>
<p>{link('bundle.json', 'Bundle inventory and SHA-256 checksums')} · {link('README.txt', 'Offline viewing and packaging notes')}</p>''')


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def write_archive(bundle: Path) -> Path:
    archive = bundle.with_name(bundle.name + ".tar.gz")
    # Publish only a completely closed gzip stream, without replacing a file
    # that appeared since preflight. Temporary and final files share a device.
    with tempfile.TemporaryDirectory(prefix=".sidis-tar-", dir=bundle.parent) as temporary:
        pending = Path(temporary) / "archive.tar.gz"
        with pending.open("xb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with tarfile.open(fileobj=compressed, mode="w") as tar:
                for path in sorted(bundle.rglob("*")):
                    if not path.is_file() or path.is_symlink():
                        continue
                    info = tar.gettarinfo(str(path), arcname=(Path(bundle.name) / path.relative_to(bundle)).as_posix())
                    info.uid = info.gid = info.mtime = 0
                    info.uname = info.gname = ""
                    info.mode = 0o644
                    with path.open("rb") as source:
                        tar.addfile(info, source)
        os.link(pending, archive)
    return archive


def build(repo: Path, campaign_root: Path, output: Path, selections: dict[str, str], *,
          formats: str = "both", allow_incomplete: bool = False,
          archive: bool = False, dry_run: bool = False) -> dict:
    repo, campaign_root, output = repo.resolve(), campaign_root.resolve(), output.absolute()
    if output.exists() or output.is_symlink():
        raise ExportError(f"Refusing to overwrite existing output: {output}")
    if contained(output, campaign_root) or contained(campaign_root, output):
        raise ExportError("Output must be separate from the campaign input tree.")
    tar_path = output.with_name(output.name + ".tar.gz")
    if archive and (tar_path.exists() or tar_path.is_symlink()):
        raise ExportError(f"Refusing to overwrite existing archive: {tar_path}")
    if formats not in {"both", "png", "pdf"}:
        raise ExportError("Plot format must be both, png or pdf.")
    if not selections:
        raise ExportError("No analyses selected.")
    items, payloads = [], {}
    for identifier, tag in sorted(selections.items()):
        if identifier not in CATALOG or not TOKEN.fullmatch(tag):
            raise ExportError(f"Invalid analysis/tag: {identifier}={tag}")
        item, files = inspect_campaign(repo, campaign_root, identifier, tag, formats)
        items.append(item)
        payloads[identifier] = files
    issues = [f"{item['id']}: {issue}" for item in items for issue in item["issues"]]
    if issues and not allow_incomplete:
        raise ExportError("Export refused; use --allow-incomplete only for an explicitly labelled partial browser:\n" + "\n".join(issues))
    report = {"analyses": len(items), "plots": sum(len(i["plots"]) for i in items),
              "payload_bytes": sum(p.stat().st_size for files in payloads.values() for p, _ in files),
              "issues": issues, "output": str(output)}
    if dry_run:
        return report
    created = datetime.now(timezone.utc).isoformat(timespec="seconds")
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".sidis-export-", dir=output.parent) as temporary:
        staging = Path(temporary) / output.name
        staging.mkdir()
        for item in items:
            target = staging / "analyses" / item["id"]
            target.mkdir(parents=True)
            for source, relative in payloads[item["id"]]:
                destination = target / relative
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(source, destination)
            sources = item.pop("source_manifest")
            write_json(target / "source-manifest.json", sources)
            write_json(target / "metadata.json", item)
            (target / "index.html").write_text(analysis_page(item), encoding="utf-8")
        (staging / "index.html").write_text(index_page(items, created), encoding="utf-8")
        (staging / "README.txt").write_text(
            "Open index.html directly in a browser after extracting the archive.\n"
            "All internal links are relative. No Python, Herwig, Rivet, web server or internet is needed to view plots.\n"
            "Only external paper/data links require internet.\n\n"
            "Included: regenerated HTML, PNG/PDF plots (selected format), central.csv when available, compact campaign/analysis metadata and reference-source checksum manifests.\n"
            "Excluded: event/shard YODA, all logs/work/build/runs/cards, original plot scripts and data helpers, dense postprocess summary/covariance, raw experimental archives, installed libraries.\n"
            "This is a viewing bundle, NOT a restart/reproduction archive. Keep the original campaigns and reference repository.\n"
            "Reference file paths inside provenance JSON identify original sources, not bundled downloads.\n"
            "bundle.json records sizes and SHA-256 of every payload file (excluding itself).\n",
            encoding="utf-8")
        inventory = [{"path": p.relative_to(staging).as_posix(), "bytes": p.stat().st_size, "sha256": digest(p)}
                     for p in sorted(staging.rglob("*")) if p.is_file()]
        write_json(staging / "bundle.json", {"created_at": created, "selections": selections,
                   "plot_format": formats, "issues": issues, "files": inventory})
        staging.rename(output)
    report["bundle_bytes"] = sum(p.stat().st_size for p in output.rglob("*") if p.is_file())
    if archive:
        tar_path = write_archive(output)
        report.update(archive=str(tar_path), archive_bytes=tar_path.stat().st_size, archive_sha256=digest(tar_path))
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=ROOT)
    parser.add_argument("--campaign-root", type=Path)
    parser.add_argument("--output", type=Path, required=True, help="New directory outside the campaign input tree")
    parser.add_argument("--measurement", action="append", choices=sorted(CATALOG), help="Repeat to select a subset; default: all 13")
    parser.add_argument("--campaign", action="append", default=[], metavar="ANALYSIS=TAG", help="Override an analysis tag; no automatic latest-tag selection")
    parser.add_argument("--plot-format", choices=("both", "png", "pdf"), default="both")
    parser.add_argument("--allow-incomplete", action="store_true", help="Label absent/partial products explicitly; signature mismatches still fail")
    parser.add_argument("--tar", action="store_true", help="Also create OUTPUT.tar.gz")
    parser.add_argument("--dry-run", action="store_true", help="Check source/products and report payload size without writing")
    args = parser.parse_args(argv)
    selections = {key: CATALOG[key][0] + "_central_20260903_v1" for key in (args.measurement or CATALOG)}
    try:
        for override in args.campaign:
            identifier, separator, tag = override.partition("=")
            if not separator or identifier not in selections or not TOKEN.fullmatch(tag):
                raise ExportError(f"Invalid override (analysis must be selected): {override}")
            selections[identifier] = tag
        result = build(args.repo, args.campaign_root or args.repo / "campaigns/phenomenology",
                       args.output, selections, formats=args.plot_format,
                       allow_incomplete=args.allow_incomplete, archive=args.tar, dry_run=args.dry_run)
        print(json.dumps(result, indent=2))
    except (ExportError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
