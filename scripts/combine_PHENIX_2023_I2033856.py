#!/usr/bin/env python3
"""Pool independent, compatible central PHENIX campaigns without generating events.

Repeat --campaign for each ORIGINAL campaign; --output-dir must be new. Supports
the nominal central expanded PHENIX setup and its finalized YODA v3 estimates.
Requires the campaign's usual Python/Rivet/YODA plotting environment.
Plots use rivet-mkhtml and the repository's existing compatibility/overlay tools.
Use --plot-only --output-dir DIRECTORY to replot an existing pool.
Add --loop-campaign DIRECTORY for independently normalized MEGammaJetLoop-only
samples. Runtime upgrades require a hash-pinned --loop-runtime-audit JSON.
Source campaigns are read only. No generation manifest is created for the pool.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import itertools
import json
import math
import re
from pathlib import Path
import shutil
import struct
import subprocess
import sys

import run_phenomenology_campaign as runner

exp = runner.experimental
MEASUREMENT = "PHENIX_2023_I2033856"
LOOP_MEASUREMENT = MEASUREMENT + "_MEGammaJetLoop-Only"
LOOP_OBJECT = "/Herwig/MatrixElements/PHENIXGammaJetLoop"
OBJECTS = {
    "inclusive_cross_section": "InclusivePhotonCrossSection",
    "isolated_cross_section": "IsolatedPhotonCrossSection",
    "isolated_all": "IsolatedPhotonYield",
}
HELICITIES = ("PP", "PM", "MP", "MM")


def read_json(path):
    return json.loads(path.read_text())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def allocated_elf_sections(path):
    """Compare actual loaded plugin sections, allowing nonloaded symbol metadata."""
    data = path.read_bytes()
    require(data[:6] == b"\x7fELF\x02\x01", f"Expected ELF64 little endian: {path}")
    offset = struct.unpack_from("<Q", data, 40)[0]
    size, count, names_index = struct.unpack_from("<HHH", data, 58)
    sections = [struct.unpack_from("<IIQQQQIIQQ", data, offset + size*i)
                for i in range(count)]
    names = sections[names_index]
    strings = data[names[4]:names[4]+names[5]]
    result = {}
    for section in sections:
        name, kind, flags, address, start, length, _, _, alignment, entry = section
        if flags & 2:  # SHF_ALLOC, including code, constants, dynamic symbols and BSS.
            label = strings[name:].split(b"\0", 1)[0].decode()
            payload = b"" if kind == 8 else data[start:start+length]
            result[label] = (kind, flags, address, length, alignment, entry,
                             hashlib.sha256(payload).hexdigest())
    require(bool(result), f"No allocated ELF sections: {path}")
    return result


def physics_identity(manifest, resolved):
    config = manifest["configuration"]
    require(manifest["measurement"] == MEASUREMENT, "Only expanded PHENIX is supported")
    require(config["families"] == ["nominal"] and
            config["variation_points"] == [[0, 0, 1.0]] and not config["smoke"],
            "Only full nominal central-scale/central-PDF campaigns are supported")
    runtime = manifest["runtime"]["provenance"]
    return {
        "measurement_signature": config["measurement_signature"],
        "resolved_measurement": resolved,
        "slices": config["hard_pt_ranges_gev"],
        "photon_sampling": config.get("photon_shower_sampling"),
        "reference": manifest["reference_snapshot"],
        "runtime_artifacts": {k: v["sha256"] for k, v in runtime["artifacts"].items()},
        "pdfs": {k: v["files"] for k, v in runtime["pdf_inventories"].items()},
        "prepared_cards_and_runs": {
            k: v["sha256"] for k, v in manifest["prepared_artifacts"]["files"].items()
            if k.startswith(("cards/", "runs/"))
        },
    }


def group_key(job):
    return (job["channel"], job["helicity"], *job["hard_pt_range_gev"])


def validate_jobs(manifest, seen_seeds, channels=("photon", "fragmentation")):
    require(manifest["status"] == "complete", "Campaign is not complete")
    config = manifest["configuration"]
    ranges = config["hard_pt_ranges_gev"]
    require(bool(ranges) and all(lo > 0 and (hi is None or hi > lo) for lo, hi in ranges),
            "Invalid hard-pT slices")
    require(all(a[1] is not None and a[1] <= b[0] for a, b in zip(ranges, ranges[1:])),
            "Overlapping or unordered hard-pT slices")
    expected = {(c, h, *bounds) for c, h, bounds in itertools.product(
        channels, HELICITIES, ranges)}
    groups, ids, outputs = {}, set(), set()
    for job in manifest["jobs"]:
        require(job["status"] == "success" and job["events"] > 0,
                f"Unsuccessful/empty job {job['id']}")
        require(job["id"] not in ids and job["output_yoda"] not in outputs,
                "Duplicate job ID or YODA output")
        ids.add(job["id"]); outputs.add(job["output_yoda"])
        require(job["seed"] not in seen_seeds, f"Repeated final seed {job['seed']}")
        seen_seeds.add(job["seed"])
        require((job["family"], job["order"], job["mpi"], job["scale"],
                 job["polarized_pdf_member"], job["unpolarized_pdf_member"],
                 job["analysis_instance"]) ==
                ("nominal", "LO", "on", 1.0, 0, 0, MEASUREMENT),
                f"Unsupported job settings: {job['id']}")
        groups.setdefault(group_key(job), []).append(job)
    require(set(groups) == expected, "Incomplete channel/helicity/slice matrix")
    for key, jobs in groups.items():
        require(len(jobs) == config["shards"] and
                {j["shard"] for j in jobs} == set(range(1, config["shards"]+1)),
                f"Missing/duplicate shards for {key}")
        index = ranges.index(list(key[2:]))
        budgets = config.get("hard_pt_events_by_channel", {}).get(key[0])
        target = budgets[index] if budgets else config["lo_events"]
        require(sum(j["events"] for j in jobs) == target, f"Event budget mismatch: {key}")


def validate_loop_runtime(baseline, loop, audit):
    """Require identical runtimes or an explicit audit of these exact binaries."""
    artifacts = lambda m: {k: v["sha256"] for k, v in
                           m["runtime"]["provenance"]["artifacts"].items()}
    old, new = artifacts(baseline), artifacts(loop)
    if old == new:
        return []
    require(audit is not None, "Loop runtime differs; supply --loop-runtime-audit")
    require(audit.get("schema_version") == 1 and
            audit.get("baseline_artifacts") == old and audit.get("loop_artifacts") == new,
            "Runtime audit does not describe these exact baseline and loop artifacts")
    changed = sorted(k for k in set(old) | set(new) if old.get(k) != new.get(k))
    require(changed == sorted(audit.get("changed_artifacts", [])),
            "Runtime audit has an incorrect changed-artifact list")
    replay = audit.get("evidence", {}).get("partner_replay", {})
    require(replay.get("status") == "PASS" and replay.get("all_compared_outputs_match") is True
            and bool(replay.get("sha256")) and bool(audit["evidence"].get("reason")),
            "Runtime audit lacks recorded successful replay evidence and an explanation")
    return changed


def card_commands(text, loop=False):
    lines = [" ".join(line.split("#", 1)[0].split()) for line in text.splitlines()]
    lines = [line for line in lines if line]
    if loop:
        expected = [f"create Herwig::MEPP2GammaJetLoop {LOOP_OBJECT} HwMEHadron.so",
                    f"set {LOOP_OBJECT}:AmplitudeMode Massive",
                    f"set {LOOP_OBJECT}:ScalePreFactor 1.0"]
        require(all(lines.count(line) == 1 for line in expected),
                "Expected the nominal Massive loop matrix element settings")
        lines = [line.replace(LOOP_OBJECT, "/Herwig/MatrixElements/MEGammaJet")
                 for line in lines if line not in expected]
    return [re.sub(r"^saverun \S+ ", "saverun COMPARISON ", line) for line in lines]


def verify_saved_card(directory, manifest, relative):
    path = (directory / relative).resolve()
    require(path.is_relative_to(directory), "Card path escapes campaign")
    require(digest(path) == manifest["prepared_artifacts"]["files"][relative]["sha256"],
            f"Saved generated card changed: {path}")
    return path.read_text()


def validate_loop_campaign(directory, manifest, resolved, base_dir, baseline, base_resolved, audit):
    """Permit only the intended hard-process replacement and audited runtime change."""
    config = manifest["configuration"]
    require(manifest["measurement"] == resolved["id"] == LOOP_MEASUREMENT,
            "--loop-campaign requires a MEGammaJetLoop-only campaign")
    require(config["families"] == ["nominal"] and
            config["variation_points"] == [[0, 0, 1.0]] and not config["smoke"] and
            not config.get("photon_shower_sampling"), "Unsupported loop settings")
    require(config["hard_pt_ranges_gev"] == baseline["configuration"]["hard_pt_ranges_gev"],
            "Loop hard-pT slices differ from the baseline")
    require(runner._signature(resolved) == config["measurement_signature"],
            "Loop measurement inputs changed since generation")
    # These fields describe the changed process or display/default choices.
    different = {"id", "title", "cards", "campaign", "channels", "families", "physics",
                 "reference", "plotting", "additive_channels"}
    require({k: v for k, v in resolved.items() if k not in different} ==
            {k: v for k, v in base_resolved.items() if k not in different},
            "Loop analysis, PDFs, scales or postprocessing differs")
    process_fields = {"accuracy_scope", "combination_status", "generation_cuts",
                      "hard_process_accuracy", "photon_production"}
    require({k: v for k, v in resolved["physics"].items() if k not in process_fields} ==
            {k: v for k, v in base_resolved["physics"].items() if k not in process_fields},
            "Loop fiducial or event-simulation settings differ")
    require(manifest["reference_snapshot"] == baseline["reference_snapshot"],
            "Loop reference snapshot differs")
    pdfs = lambda m: {k: v["files"] for k, v in
                      m["runtime"]["provenance"]["pdf_inventories"].items()}
    require(pdfs(manifest) == pdfs(baseline), "Loop PDF files differ")
    require(set(resolved["channels"]) == {"photon"} and not resolved.get("additive_channels"),
            "Loop input contains additional hard-process channels")
    for field in ("raw_objects", "observables", "scale_interpretation"):
        require(resolved["channels"]["photon"][field] == base_resolved["channels"]["photon"][field],
                f"Loop observable definition differs: {field}")
    require(resolved["channels"]["photon"]["scale_object"] == LOOP_OBJECT,
            "Unexpected loop matrix-element object")
    require({k: v for k, v in resolved["families"]["nominal"].items() if k != "label"} ==
            {k: v for k, v in base_resolved["families"]["nominal"].items() if k != "label"},
            "Loop helicity or MPI settings differ")
    common = "cards/" + base_resolved["cards"]["common"]
    require(resolved["cards"]["common"] == base_resolved["cards"]["common"] and
            verify_saved_card(directory, manifest, common) ==
            verify_saved_card(base_dir, baseline, common), "Loop common card differs")
    base_jobs = {group_key(j): j for j in baseline["jobs"] if j["channel"] == "photon"}
    loop_jobs = {group_key(j): j for j in manifest["jobs"]}
    require(set(loop_jobs) == set(base_jobs), "Loop helicity/slice matrix differs")
    for key, job in loop_jobs.items():
        original = verify_saved_card(base_dir, baseline, "cards/" + base_jobs[key]["card_input"])
        added = verify_saved_card(directory, manifest, "cards/" + job["card_input"])
        require(card_commands(original) == card_commands(added, loop=True),
                f"Loop generated card has additional physics changes: {key}")
    return validate_loop_runtime(baseline, manifest, audit)


def read_shard(item):
    directory, job = item
    path = (directory / job["output_yoda"]).resolve()
    require(path.is_relative_to(directory), f"YODA path escapes campaign: {path}")
    requested = {f"/{MEASUREMENT}/{obj}": name for name, obj in OBJECTS.items()}
    requested["/RAW/_EVTCOUNT"] = "count"
    blocks, block = {}, None
    for line in path.read_text().splitlines():
        if line.startswith("BEGIN "):
            _, kind, name = line.split()
            block = {"kind": kind, "rows": []} if name in requested else None
            if block is not None:
                require(requested[name] not in blocks, f"Duplicate YODA object: {path}")
                blocks[requested[name]] = block
        elif line.startswith("END "):
            block = None
        elif block is not None:
            if line.startswith("Edges(A1):"):
                block["edges"] = json.loads(line.split(":", 1)[1])
            elif line.startswith("ErrorLabels:"):
                block["labels"] = json.loads(line.split(":", 1)[1])
            elif line.strip() == "---":
                block["data"] = True
            elif block.get("data") and line.strip() and not line.startswith("#") and ":" not in line:
                block["rows"].append([0.0 if x == "---" else float(x) for x in line.split()])
    require(set(blocks) == {*OBJECTS, "count"}, f"Missing YODA objects: {path}")
    count = blocks.pop("count")
    require(count["kind"] == "YODA_COUNTER_V3" and len(count["rows"]) == 1 and
            count["rows"][0][-1] == job["events"], f"Incomplete event count: {path}")
    objects = {}
    for name, block in blocks.items():
        require(block["kind"] == "YODA_ESTIMATE1D_V3" and block.get("labels", []) in ([], ["stats"]),
                f"Expected finalized YODA v3 estimates with only stats errors: {path}")
        edges, rows = block["edges"], block["rows"]
        require(len(rows) == len(edges)+1 and all(a < b for a, b in zip(edges, edges[1:])),
                f"Unexpected binning: {path}")
        rows = rows[1:-1]  # Finalized v3 estimates include underflow/overflow.
        if not block.get("labels"):
            # YODA omits ErrorLabels/columns when the whole histogram is empty.
            require(all(r == [0.0] for r in rows), f"Nonempty estimate without errors: {path}")
            rows = [[0.0, 0.0, 0.0] for _ in rows]
        require(all(len(r) == 3 and all(math.isfinite(v) for v in r) and
                    r[1] <= 0 <= r[2] and math.isclose(-r[1], r[2], rel_tol=1e-10, abs_tol=1e-14)
                    for r in rows), f"Invalid/asymmetric statistical estimate: {path}")
        objects[name] = exp.BinSeries(edges, [r[0] for r in rows], [r[2]**2 for r in rows])
    return group_key(job), job["events"], objects


def pool_groups(reports):
    """Event fractions for overlapping strata; squared fractions for variance."""
    require(bool(reports) and all(set(r) == set(reports[0]) for r in reports),
            "Cannot pool different channel/helicity/slice groups")
    groups = {}
    for key in reports[0]:
        parts = [report[key] for report in reports]
        counts = [part["events"] for part in parts]
        groups[key] = {"events": sum(counts), "objects": {
            name: exp.combine_shard_series([part["objects"][name] for part in parts], counts)
            for name in OBJECTS}}
    return groups


def add_loop_groups(baseline, loop_reports):
    """Pool repeated loop samples separately, then add their physical component."""
    loop = pool_groups(loop_reports)
    require(all(key[0] == "photon" for key in loop), "Unexpected loop channel")
    renamed = {("gammajet_loop", *key[1:]): value for key, value in loop.items()}
    require(not (set(renamed) & set(baseline)), "Loop component already present")
    return {**baseline, **renamed}


def predict(groups):
    # Disjoint hard-pT slices and the established additive channels are summed;
    # the existing PHENIX routine then averages cross sections/forms A_LL.
    samples = {name: {h: exp.add_independent_series([
        group["objects"][name] for key, group in groups.items() if key[1] == h
    ]) for h in HELICITIES} for name in OBJECTS}
    return {k: v for k, v in runner._phenix_prediction(samples).items() if k in OBJECTS}


def collect(directory, manifest):
    shards = {}
    with ThreadPoolExecutor(max_workers=8) as pool:
        for i, (key, events, objects) in enumerate(pool.map(
                read_shard, ((directory, j) for j in manifest["jobs"])), 1):
            shards.setdefault(key, []).append({"events": events, "objects": objects})
            if i % 2000 == 0:
                print(f"{directory.name}: read {i}/{len(manifest['jobs'])} shards", flush=True)
    groups = {key: pool_groups([{key: part} for part in parts])[key]
              for key, parts in shards.items()}
    prediction = predict(groups)
    saved = read_json(directory / "postprocess/summary.json")["variations"][
        "nominal-photon-p000-u000-mu1-mpion"]
    for name, series in prediction.items():
        for field in ("edges", "values", "errors"):
            require(len(series[field]) == len(saved[name][field]) and all(
                math.isclose(x, y, rel_tol=1e-9, abs_tol=1e-12)
                for x, y in zip(series[field], saved[name][field])),
                f"Recomputed {name}/{field} differs from saved results: {directory}")
    return groups, prediction


def write_yoda(path, prediction, reference, annotations=None):
    """Use the campaign's native YODA writer and statistical-error convention."""
    yoda = exp._import_yoda()
    objects = [
        runner._estimate_with_bands(
            yoda, series, reference["channels"][name]["rivet_path"],
            {"Generator": "HerwigPol full event simulation",
             "Combination": "independent campaign pool", **(annotations or {})})
        for name, series in prediction.items()
    ]
    yoda.write(objects, str(path), precision=17)


def configure_dis_font(html, text_rendering):
    """Keep Rivet's DIS serif typography when TeX is unavailable locally."""
    if text_rendering["mode"] != "mathtext":
        return
    from matplotlib import font_manager

    # P052 and TeX Gyre Pagella are Palatino-compatible serif fonts.
    available = {font.name for font in font_manager.fontManager.ttflist}
    family = next((name for name in ("Palatino", "TeX Gyre Pagella", "P052")
                   if name in available), "DejaVu Serif")
    style_path = html / "default.mplstyle"
    style = style_path.read_text()
    for key, value in {
        "font.family": "serif", "font.serif": family,
        "mathtext.rm": family, "mathtext.bf": family + ":bold",
        "mathtext.it": family + ":italic", "mathtext.cal": "STIXGeneral:italic",
    }.items():
        pattern = rf"(?m)^{re.escape(key)}\s*:.*$"
        if re.search(pattern, style):
            style = re.sub(pattern, lambda _: f"{key}: {value}", style)
        else:
            style += f"\n{key}: {value}\n"
    style_path.write_text(style)
    for script in html.rglob("*.py"):
        text = script.read_text()
        # YODA's no-LaTeX fallback otherwise overrides the style with sans serif.
        text = text.replace("plt.rc('font', family='DejaVu Sans')",
                            "plt.rc('font', family='serif')")
        script.write_text(text)
    text_rendering["font_family"] = family


def configure_residual_bounds(script):
    """Keep the standard +/-3 sigma view unless a residual/error needs more room."""
    source = script.read_text()
    marker = "\nplt.savefig("
    require(marker in source and "ratio0_ax" in source, "Missing residual panel or save marker")
    adjustment = r'''
# Include every finite residual and uncertainty, with padding; retain +/-3 sigma.
_residual_extent = [0.0]
for _label, _values in dataf.get('ratio0_yvals', {}).items():
    _values = np.asarray(_values, dtype=float)
    _residual_extent.extend(_values[np.isfinite(_values)].tolist())
    for _errors in dataf.get('ratio0_yerrs', {}).get(_label, []):
        _errors = np.abs(np.asarray(_errors, dtype=float))
        for _sign in (-1, 1):
            _ends = _values + _sign*_errors
            _residual_extent.extend(_ends[np.isfinite(_ends)].tolist())
_residual_min, _residual_max = min(_residual_extent), max(_residual_extent)
_residual_pad = 0.08*max(_residual_max-_residual_min, 1.0)
_residual_limits = (min(-3.0, _residual_min-_residual_pad), max(3.0, _residual_max+_residual_pad))
ratio0_ax.set_ylim(*_residual_limits)
_residual_ticks = range(int(np.ceil(_residual_limits[0])), int(np.floor(_residual_limits[1])))
ratio0_ax.set_yticks(_residual_ticks, [rf'${i}\,\sigma$' for i in _residual_ticks])
'''
    script.write_text(source.replace(marker, adjustment + marker, 1))


def plot_results(output, rivet_mkhtml=None):
    """Plot pooled YODA through the same Rivet path as runner.plot_pp."""
    summary = read_json(output / "summary.json")
    require(summary.get("kind") == "independent_campaign_pool" and
            summary.get("measurement") == MEASUREMENT, "Not a PHENIX campaign pool")
    measurement = read_json(output / "resolved-measurement.json")
    require(runner._signature(measurement) == summary["measurement_signature"],
            "Measurement inputs changed since pooling; refusing a different plot definition")
    reference = read_json(output / "reference.json")
    require(reference == runner._measurement_snapshot(measurement), "Reference data changed")
    tool = str(rivet_mkhtml) if rivet_mkhtml else shutil.which("rivet-mkhtml")
    require(tool and Path(tool).is_file(),
            "rivet-mkhtml not found: activate the Rivet environment or pass --rivet-mkhtml PATH")
    tool = str(Path(tool).resolve())
    prediction = output / "prediction.yoda"
    objects = exp._import_yoda().read(str(prediction))
    expected_paths = {reference["channels"][name]["rivet_path"] for name in OBJECTS}
    require(set(objects) == expected_paths, "Unexpected objects in pooled YODA")
    for name, series in summary["predictions"].items():
        obj = objects[reference["channels"][name]["rivet_path"]]
        require(list(obj.xEdges()) == series["edges"], "Pooled YODA binning differs from summary")
        bins = list(obj.bins())
        for field, actual in (("values", [b.val() for b in bins]),
                              ("errors", [b.errAvg("stat") for b in bins])):
            # YODA 2.1's Estimate1D writer uses seven significant digits even
            # when a higher general writer precision is requested.
            require(len(actual) == len(series[field]) and all(
                math.isclose(x, y, rel_tol=5.1e-7, abs_tol=1e-14)
                for x, y in zip(actual, series[field])),
                f"Pooled YODA {name}/{field} differs from summary")
    html = output / "plots" / "html"
    logs = output / "logs"
    logs.mkdir(exist_ok=True)
    if html.exists():
        shutil.rmtree(html)  # Only this helper's generated Rivet plot subtree.
    html.mkdir(parents=True)
    environment = exp.analysis_environment(measurement, output, {})
    reference_dir = str((runner.DISPOL_ROOT / measurement["analysis"]["reference_yoda"]).parent)
    environment["RIVET_DATA_PATH"] = reference_dir + ":" + environment.get("RIVET_DATA_PATH", "")
    cache = output / "work" / "matplotlib-cache"
    cache.mkdir(parents=True, exist_ok=True)
    environment["MPLCONFIGDIR"] = str(cache)
    environment["DISPOL_FORCE_NO_ERROR_BANDS"] = "1"
    label = ("LO+PS + loop-induced (polarized; full spin)" if summary.get("loop_sources")
             else "LO+PS (polarized; full spin)")
    # A display override leaves the original generation/measurement signatures intact.
    display_config = output / "work" / "combined-display.plot"
    display_config.write_text(
        f"# BEGIN PLOT /{MEASUREMENT}/d.*\n"
        "RatioPlot=1\n"
        "# END PLOT\n")
    command = [
        sys.executable, str(runner.DISPOL_ROOT / "scripts/rivet_mkhtml_safe.py"),
        tool, "--dry-run", "--offline", "--deviation", "--pwd", "-o", str(html),
        "--ratiolabel", measurement["plotting"]["ratio_label"],
        "-M", r".*/DIAGNOSTICS/.*", "-c", str(display_config),
        runner._prediction_plot_argument(prediction, label, measurement["families"]["nominal"]),
    ]
    exp._run_logged(command, output, environment, logs / "rivet-mkhtml-generate.log")
    text_rendering = runner._configure_plot_text_rendering(html, environment)
    configure_dis_font(html, text_rendering)
    scripts = sorted(p for p in html.rglob("*.py") if not p.name.endswith("__data.py"))
    require({p.stem for p in scripts} == {Path(p).name for p in expected_paths},
            "Rivet did not produce the three expected PHENIX plot scripts")
    by_stem = {Path(reference["channels"][name]["rivet_path"]).name: series
               for name, series in summary["predictions"].items()}
    with (logs / "rivet-plot-scripts.log").open("w") as log:
        log.write("text_rendering: " + json.dumps(text_rendering) + "\n")
        for script in scripts:
            # Identical MC-statistics and experimental overlays to plot_pp.
            exp.add_theory_uncertainty_overlay(
                script, {"monte_carlo": by_stem[script.stem]["errors"]}, prediction.name)
            exp.add_experimental_error_overlay(
                script, runner._pp_reference_overlay_points(measurement, reference, script.stem),
                show_statistical=False)
            configure_residual_bounds(script)
            require(exp.plot_script_has_finite_y(script), f"Non-renderable plot: {script}")
            log.write(f"script: {script}\n"); log.flush()
            completed = subprocess.run([sys.executable, str(script)], cwd=script.parent,
                                       env=environment, stdout=log, stderr=subprocess.STDOUT)
            require(completed.returncode == 0,
                    f"Rivet plot script failed: {script}; see {logs / 'rivet-plot-scripts.log'}")
            require(script.with_suffix(".png").is_file() and script.with_suffix(".pdf").is_file(),
                    f"Missing rendered Rivet output: {script}")
    index = exp.write_plot_indexes(html, measurement, scripts)
    summary["plots"] = {
        "renderer": "rivet-mkhtml via rivet_mkhtml_safe.py",
        "index": str(index.relative_to(output)), "created_at": exp.utc_now(),
        "rivet_mkhtml": tool, "text_rendering": text_rendering,
        "plot_metadata_sha256": digest(runner.DISPOL_ROOT / measurement["analysis"]["plot"]),
        "prediction_sha256": digest(prediction),
        "legend_label": label,
        "display_config": str(display_config.relative_to(output)),
        "display_config_sha256": digest(display_config),
        "lower_panel": "normalized residual (MC-Data)/sigma_Data",
        "residual_bounds": "automatic with 8% padding; at least [-3,3]",
        "plotting_script_sha256": digest(Path(__file__)),
    }
    exp.atomic_write_json(output / "summary.json", summary)
    if (output / "provenance.json").exists():
        provenance = read_json(output / "provenance.json")
        provenance["plotting"] = summary["plots"]
        provenance["summary_sha256"] = digest(output / "summary.json")
        exp.atomic_write_json(output / "provenance.json", provenance)
    print(f"Wrote Rivet HTML to {index}")
    return index


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, action="append")
    parser.add_argument("--loop-campaign", type=Path, action="append",
                        help="Add a completed MEGammaJetLoop-only campaign as a separate component")
    parser.add_argument("--loop-runtime-audit", type=Path,
                        help="JSON recording exact runtime hashes and evidence for a loop-era upgrade")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--no-plots", action="store_true", help="Write numerical outputs without invoking Rivet plotting")
    parser.add_argument("--plot-only", action="store_true", help="Replot the saved combined YODA in --output-dir")
    parser.add_argument("--rivet-mkhtml", type=Path, help="Optional path to rivet-mkhtml; default: activated PATH")
    args = parser.parse_args(argv)
    directories = [p.expanduser().resolve() for p in args.campaign or []]
    loop_dirs = [p.expanduser().resolve() for p in args.loop_campaign or []]
    all_dirs = directories + loop_dirs
    output = args.output_dir.expanduser().resolve()
    try:
        if args.plot_only:
            require(not all_dirs and not args.no_plots and not args.loop_runtime_audit,
                    "--plot-only cannot be combined with campaign inputs, an audit or --no-plots")
            plot_results(output, args.rivet_mkhtml)
            return 0
        exp._import_yoda()
        require(len(directories) >= 2 and len(set(directories)) == len(directories),
                "Supply at least two distinct original campaign directories")
        require(len(set(all_dirs)) == len(all_dirs), "Duplicate campaign directory")
        require(loop_dirs or not args.loop_runtime_audit, "Runtime audit requires --loop-campaign")
        require(not output.exists() and all(not output.is_relative_to(p) for p in all_dirs),
                "Output must be a new directory outside all input campaigns")
        manifests = [read_json(p / "manifest.json") for p in directories]
        resolved = [read_json(p / "resolved-measurement.json") for p in directories]
        identities = [physics_identity(m, r) for m, r in zip(manifests, resolved)]
        require(all(i == identities[0] for i in identities),
                "Campaign physics, cards, runtime, PDF, reference or slicing differs")
        require(runner._signature(resolved[0]) == identities[0]["measurement_signature"],
                "Baseline measurement inputs changed since generation")
        loop_manifests = [read_json(p / "manifest.json") for p in loop_dirs]
        loop_resolved = [read_json(p / "resolved-measurement.json") for p in loop_dirs]
        audit = read_json(args.loop_runtime_audit) if args.loop_runtime_audit else None
        runtime_changes = [validate_loop_campaign(p, m, r, directories[0], manifests[0], resolved[0], audit)
                           for p, m, r in zip(loop_dirs, loop_manifests, loop_resolved)]
        # Repeated loop campaigns must also agree with each other, including cards/runtime.
        if loop_manifests:
            def loop_identity(m, r):
                base_view = dict(m, measurement=MEASUREMENT)
                return physics_identity(base_view, r)
            require(all(loop_identity(m, r) == loop_identity(loop_manifests[0], loop_resolved[0])
                        for m, r in zip(loop_manifests, loop_resolved)),
                    "Repeated loop campaign settings differ")
        plugins = [p / "build" / r["analysis"]["plugin"]
                   for p, r in zip(all_dirs, resolved + loop_resolved)]
        for plugin, manifest in zip(plugins, manifests + loop_manifests):
            require(digest(plugin) == manifest["plugin_provenance"]["sha256"], "Saved plugin changed")
        sections = [allocated_elf_sections(p) for p in plugins]
        require(all(s == sections[0] for s in sections), "Loaded Rivet plugin sections differ")
        seen_seeds = set()
        for manifest in manifests:
            validate_jobs(manifest, seen_seeds)
        for manifest in loop_manifests:
            validate_jobs(manifest, seen_seeds, channels=("photon",))
        reference_path = runner.DISPOL_ROOT / manifests[0]["reference_snapshot"]["path"]
        require(digest(reference_path) == manifests[0]["reference_snapshot"]["sha256"],
                "Reference snapshot changed")
        reference = read_json(reference_path)
        reports, loop_reports, sources = [], [], []
        for directory, manifest in zip(all_dirs, manifests + loop_manifests):
            groups, prediction = collect(directory, manifest)
            is_loop = directory in loop_dirs
            (loop_reports if is_loop else reports).append(groups)
            sources.append({"path": str(directory), "tag": manifest["tag"],
                            "role": "additive_loop_component" if is_loop else "repeated_expanded_campaign",
                            "manifest_sha256": digest(directory / "manifest.json"),
                            "resolved_measurement_sha256": digest(directory / "resolved-measurement.json"),
                            "saved_summary_sha256": digest(directory / "postprocess/summary.json"),
                            "measurement_signature": manifest["configuration"]["measurement_signature"],
                            "events": sum(j["events"] for j in manifest["jobs"]),
                            "jobs": len(manifest["jobs"]), "prediction": prediction})
        groups = pool_groups(reports)
        baseline_prediction = predict(groups)
        if loop_reports:
            groups = add_loop_groups(groups, loop_reports)
        prediction = predict(groups)
        summary = {"measurement": MEASUREMENT, "kind": "independent_campaign_pool",
                   "created_at": exp.utc_now(), "sources": sources,
                   "total_events": sum(s["events"] for s in sources),
                   "measurement_signature": identities[0]["measurement_signature"],
                   "photon_sampling": identities[0]["photon_sampling"],
                   "normalization": "Event-count weighted within matching channel/helicity/slice; squared fractions for variances. Add disjoint slices and channels before forming A_LL.",
                   "uncertainty": "MC statistics only, using stored event-level sumW2; no inter-bin covariance",
                   "physics": resolved[0]["physics"], "predictions": prediction,
                   "groups": [{"channel": k[0], "helicity": k[1], "slice": list(k[2:]),
                               "events": g["events"], "objects": {
                                   n: {"edges": s.edges, "values": s.values, "variances": s.variances}
                                   for n, s in g["objects"].items()}} for k, g in groups.items()],
                   "checks": {"independent_final_seeds": True,
                              "physics_and_runtime_match": not any(runtime_changes),
                              "repeated_base_campaigns_physics_and_runtime_match": True,
                              "loaded_plugin_sections_match": True, "all_shards_complete": True,
                              "individual_saved_predictions_reproduced": True}}
        if loop_reports:
            summary.update({
                "loop_sources": [m["tag"] for m in loop_manifests],
                "baseline_prediction": baseline_prediction,
                "loop_prediction": predict(pool_groups(loop_reports)),
                "loop_runtime_changes": runtime_changes,
                "loop_runtime_audit": audit,
                "loop_physics": loop_resolved[0]["physics"],
                "accuracy_scope": "Existing unmerged tree photon+jet and QCD-dijet+PS sum, plus the leading gg -> gamma g loop-induced contribution with QCD/QED PS; not a complete NLO or NNLO calculation.",
            })
            summary["checks"]["loop_settings_and_generated_cards_checked"] = True
            summary["checks"]["runtime_differences_hash_pinned"] = bool(audit) if any(runtime_changes) else True
        serialized = json.dumps(summary, indent=2, allow_nan=False) + "\n"
        output.mkdir(parents=True)
        (output / "summary.json").write_text(serialized)
        (output / "reference.json").write_text(json.dumps(reference, indent=2) + "\n")
        (output / "resolved-measurement.json").write_text(json.dumps(resolved[0], indent=2) + "\n")
        write_yoda(output / "prediction.yoda", prediction, reference, {
            "SourceCampaigns": ",".join(s["tag"] for s in sources),
            "GeneratedEvents": summary["total_events"],
            "PhotonShowerSampling": json.dumps(summary["photon_sampling"], sort_keys=True),
            "HardComponents": "MEGammaJet,MEQCD2to2" + (",MEPP2GammaJetLoop" if loop_reports else ""),
        })
        provenance = {
            "kind": "postprocessing_combination", "created_at": summary["created_at"],
            "sources": sources, "total_events": summary["total_events"],
            "total_jobs": sum(s["jobs"] for s in sources),
            "normalization": summary["normalization"], "uncertainty": summary["uncertainty"],
            "script": {"path": str(Path(__file__).resolve()), "sha256": digest(Path(__file__))},
            "command": [sys.executable, str(Path(__file__).resolve()), *(argv if argv is not None else sys.argv[1:])],
            "runtime_audit": audit,
            "runtime_audit_sha256": digest(args.loop_runtime_audit) if args.loop_runtime_audit else None,
            "prediction_sha256": digest(output / "prediction.yoda"),
            "summary_sha256": digest(output / "summary.json"), "checks": summary["checks"],
        }
        exp.atomic_write_json(output / "provenance.json", provenance)
        if audit:
            exp.atomic_write_json(output / "loop-runtime-audit.json", audit)
        with (output / "central.csv").open("w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["observable", "low", "high", "value", "mc_stat"])
            for name, series in prediction.items():
                for low, high, value, error in zip(series["edges"], series["edges"][1:], series["values"], series["errors"]):
                    writer.writerow([name, low, high, value, error])
        if not args.no_plots:
            plot_results(output, args.rivet_mkhtml)
        print(f"Pooled {summary['total_events']:,} events from {len(sources)} campaigns into {output}")
        return 0
    except (ValueError, KeyError, OSError, ImportError, exp.CampaignError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
