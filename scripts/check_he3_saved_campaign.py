#!/usr/bin/env python3
"""Reprocess and plot saved helium campaigns in new, independent validation copies.

No Herwig, event generation, or plugin build is invoked. Requires the active
Rivet/YODA Python environment for execution; --dry-run needs only Python.
"""
from __future__ import annotations

import argparse
import ast
import copy
import json
import math
import re
import shutil
import sys
from pathlib import Path

import run_experimental_campaign as campaign

IDS = {"E142_1996_I424108", "E154_1997_I443170", "HERMES_1997_I440904"}


def require(condition, message):
    if not condition:
        raise campaign.CampaignError(message)


def contained_file(root, relative):
    path = Path(relative)
    require(not path.is_absolute() and '..' not in path.parts, f"Unsafe campaign path: {relative}")
    result = (root / path).resolve()
    require(result.is_relative_to(root.resolve()), f"Input escapes campaign directory: {relative}")
    require(result.is_file() and result.stat().st_size > 0, f"Missing or empty input: {result}")
    return result


def compare(expected, actual, where="result", *, rtol=1e-9, atol=1e-12):
    """Compare scientific content, ignoring no fields within the supplied tree."""
    if isinstance(expected, dict):
        require(isinstance(actual, dict) and expected.keys() == actual.keys(), f"{where}: different keys")
        for key in expected:
            compare(expected[key], actual[key], f"{where}/{key}", rtol=rtol, atol=atol)
    elif isinstance(expected, (list, tuple)):
        require(isinstance(actual, (list, tuple)) and len(expected) == len(actual), f"{where}: different lengths")
        for index, (left, right) in enumerate(zip(expected, actual)):
            compare(left, right, f"{where}[{index}]", rtol=rtol, atol=atol)
    elif isinstance(expected, (int, float)) and not isinstance(expected, bool):
        require(isinstance(actual, (int, float)), f"{where}: expected numeric value")
        equal = (math.isnan(expected) and math.isnan(actual)) or expected == actual
        require(equal or math.isclose(expected, actual, rel_tol=rtol, abs_tol=atol),
                f"{where}: {actual!r} != baseline {expected!r}")
    else:
        require(expected == actual, f"{where}: {actual!r} != baseline {expected!r}")


def summary_science(summary):
    """Compare old/new helium summary schemas without dropping scientific fields."""
    result = {key: summary[key] for key in
              ('uncertainties', 'masked_bins', 'diagnostics', 'variation_points')}
    rows = []
    for original in summary['bins']:
        row = dict(original)
        require('x_low' in row and 'x_high' in row, 'Helium summary needs x-bin edges')
        # The original helium format implied an x axis and support for every bin.
        # Fill those defaults, but retain the new fields in the comparison so a
        # changed axis, boundary or support flag cannot pass unnoticed.
        row.setdefault('axis', 'x')
        row.setdefault('bin_low', row['x_low'])
        row.setdefault('bin_high', row['x_high'])
        row.setdefault('supported', True)
        require(row['axis'] == 'x', 'Unexpected helium summary axis')
        require(type(row['supported']) is bool, 'Helium bin support must be boolean')
        compare(row['x_low'], row['bin_low'], 'summary lower-edge aliases')
        compare(row['x_high'], row['bin_high'], 'summary upper-edge aliases')
        # This human-readable description is now copied into every row. The
        # acceptance definition is checked by the campaign signature preflight.
        row.pop('acceptance', None)
        rows.append(row)
    result['bins'] = sorted(rows, key=lambda row: (row['selection'], row['bin']))
    return result


def estimate_science(obj, sources=None):
    result = {'edges': list(obj.xEdges()), 'bins': []}
    for index, bin_obj in enumerate(obj.bins(False, True), 1):
        if obj.isMasked(index):
            result['bins'].append({'masked': True})
        else:
            names = set(bin_obj.sources())
            if sources is not None:
                names &= set(sources)
            result['bins'].append({'value': float(bin_obj.val()), 'errors': {
                name: [float(bin_obj.errDown(name)), float(bin_obj.errUp(name))]
                for name in sorted(names)}})
    return result


def numeric_assignments(path):
    """Read generated numeric dictionaries without executing generated Python."""
    class Numbers(ast.NodeTransformer):
        def visit_Name(self, node):
            if node.id in {'nan', 'inf'}:
                return ast.copy_location(ast.Constant(float(node.id)), node)
            return node
    output = {}
    for node in ast.parse(path.read_text()).body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            try:
                output[node.targets[0].id] = ast.literal_eval(Numbers().visit(node.value))
            except (ValueError, TypeError):
                pass
    return output


def check_plot(script, expected, reference, bands):
    require(script.is_file(), f"Missing plot script: {script}")
    for suffix in ('.png', '.pdf'):
        path = script.with_suffix(suffix)
        require(path.is_file() and path.stat().st_size > 0, f"Missing rendered output: {path}")
    assignments = numeric_assignments(script)
    for variable, key in (('_theory_mc', 'monte_carlo'), ('_theory_pdf', 'pdf_68'),
                          ('_theory_scale_dn', 'scale_down'), ('_theory_scale_up', 'scale_up')):
        values = bands.get(key)
        expected_band = ([float('nan') if value is None else value for value in values]
                         if values is not None and any(value is not None for value in values) else None)
        require(variable in assignments, f'Missing uncertainty overlay: {script}/{variable}')
        compare(expected_band, assignments[variable], f'{script.name}/{variable}')
    data = numeric_assignments(script.with_name(script.stem + '__data.py'))
    labels = [key for key in data.get('yvals', {}) if key != 'Data']
    require(len(labels) == 1 and 'Data' in data['yvals'], f"Expected data and one MC curve in {script}")
    label = labels[0]
    compare([point['value'] for point in reference], data['yvals']['Data'], 'plot data', rtol=3e-5, atol=3e-5)
    for index, result in enumerate(expected['bins']):
        if result.get('masked'):
            require(not math.isfinite(data['yvals'][label][index]), f"Masked bin became finite: {script}/{index}")
            continue
        compare(result['value'], data['yvals'][label][index], 'plot MC', rtol=3e-5, atol=3e-5)
        sigma = data['yerrs']['Data'][0][index]
        require(sigma > 0, f"Nonpositive reference uncertainty: {script}/{index}")
        residual = (data['yvals'][label][index] - data['yvals']['Data'][index]) / sigma
        compare(residual, data['ratio0_yvals'][label][index], 'normalized residual', rtol=3e-5, atol=3e-5)
    compare([0.] * len(reference), data['ratio0_yvals']['Data'], 'residual data baseline')
    from rivet_mkhtml_safe import RESIDUAL_LABEL
    labels = [node.args[0].value for node in ast.walk(ast.parse(script.read_text()))
              if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
              and node.func.attr == 'set_ylabel' and node.args
              and isinstance(node.args[0], ast.Constant)]
    require(RESIDUAL_LABEL in labels, f"Missing or incorrectly escaped residual label: {script}")
    require('\x0c' not in script.read_text(), f"Invalid TeX escape in {script}")


def inspect(source, output_root):
    source = source.expanduser().resolve()
    require(not output_root.is_relative_to(source) and not source.is_relative_to(output_root),
            'The output directory must be separate from every source campaign')
    manifest_path = contained_file(source, 'manifest.json')
    manifest = campaign.load_json(manifest_path)
    identifier, tag = manifest['measurement'], manifest['tag']
    require(identifier in IDS, f"This check supports helium DIS only: {identifier}")
    require(campaign.TAG_PATTERN.fullmatch(tag), f"Invalid saved tag: {tag}")
    measurement = campaign.get_measurement(identifier)
    compatibility = campaign._assert_manifest_signatures_current(manifest, measurement, campaign_dir=source)
    groups = campaign._require_complete_matrix(copy.deepcopy(manifest), measurement, source)
    require(not manifest['configuration'].get('comparisons'),
            'Use a nominal-family campaign for this helium integration check')
    post = manifest.get('postprocess', {})
    baseline = campaign.load_json(contained_file(source, post.get('summary_json', 'postprocess/summary.json')))
    require(baseline.get('order_sign_audit', {}).get('order_coefficients'),
            'Baseline has no recorded NLO sign audit; use a corrected campaign')
    inputs = {'manifest.json', 'resolved-measurement.json', post['summary_json'], post['yoda']}
    if post.get('variations_yoda'):
        inputs.add(post['variations_yoda'])
    inputs.update(str(job['output_yoda']) for jobs in groups.values() for job in jobs)
    files = {relative: contained_file(source, relative) for relative in sorted(inputs)}
    require(contained_file(campaign.DISPOL_ROOT, measurement['analysis']['reference_yoda']).is_file(),
            'Pinned reference YODA must exist before validation')
    return {'source': source, 'manifest': manifest, 'measurement': measurement,
            'compatibility': compatibility, 'files': files, 'baseline': baseline,
            'destination': output_root / identifier / tag,
            'copy_bytes': sum(path.stat().st_size for path in files.values())}


def public_plan(plan):
    return {'source': str(plan['source']), 'destination': str(plan['destination']),
            'input_files': len(plan['files']), 'copy_bytes': plan['copy_bytes'],
            'jobs': len(plan['manifest']['jobs']), 'variation_points': plan['manifest']['configuration'].get('variation_points'),
            'compatibility': plan['compatibility'], 'generates_events': False}


def run_check(plan, tool):
    directory, measurement = plan['destination'], plan['measurement']
    directory.mkdir(parents=True, exist_ok=False)
    hashes, copied = {}, []
    manifest = copy.deepcopy(plan['manifest'])
    post = manifest['postprocess']
    baseline_paths = {post['summary_json'], post['yoda'], post.get('variations_yoda')}
    print(f"Copying {len(plan['files'])} files for {measurement['id']} ({plan['copy_bytes']/1e6:.1f} MB)", flush=True)
    for relative, source in plan['files'].items():
        destination = directory / ('baseline' if relative in baseline_paths else '') / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)  # Independent files, never hardlinks/symlinks.
        sha = campaign.sha256_file(destination)
        require(sha == campaign.sha256_file(source), f"Input changed during copy: {source}")
        hashes[relative] = sha
        copied.append(str(destination.relative_to(directory)))
    report = {**public_plan(plan), 'input_sha256': hashes, 'copied_files': copied,
              'status': 'running', 'runner_sha256': campaign.sha256_file(campaign.SCRIPT_PATH),
              'checker_sha256': campaign.sha256_file(Path(__file__)), 'plot_tool': str(tool),
              'plot_tool_sha256': campaign.sha256_file(tool), 'python': sys.executable,
              'python_version': sys.version}
    report_path = directory / 'validation.json'
    campaign.atomic_write_json(report_path, report)
    previous_root = campaign.CAMPAIGN_ROOT
    try:
        # Presentation-only override, limited to validation: exercise residuals
        # even though the committed helium .plot cards disable lower panels.
        controls = directory / 'validation-controls'
        controls.mkdir()
        config = controls / 'residual.plot'
        config.write_text(f'BEGIN PLOT /{measurement["id"]}/.*\nRatioPlot=1\nRatioPlotYLabel=MC/Data\nEND PLOT\n')
        paths = [out['a1'] for out in measurement['outputs'].values() if 'a1' in out and not out.get('summary_only')]
        # Rivet applies --match to MC paths and /REF-prefixed data paths alike.
        pattern = '^(?:/REF)?(?:' + '|'.join(re.escape(path) for path in paths) + ')$'
        shim = controls / 'rivet-mkhtml-validation.py'
        shim.write_text('import runpy, sys\n'
                        + f'sys.argv[1:1] = ["--config", {str(config)!r}, "--match", {pattern!r}]\n'
                        + f'sys.argv[0] = {str(tool)!r}\nrunpy.run_path({str(tool)!r}, run_name="__main__")\n')
        manifest['runtime']['tools']['rivet-mkhtml'] = str(shim)
        manifest['runtime']['rivet_plugin_compiler'] = ''
        manifest['validation_copy'] = {'source': str(plan['source']), 'residual_override': str(config),
                                       'active_plot_tool': str(tool), 'only_configured_a1_plots': True}
        # Do not retain pointers to old rendered products inside this copy.
        for key in list(manifest):
            if key == 'postprocess' or key.startswith('plots'):
                manifest.pop(key)
        campaign.atomic_write_json(directory / campaign.MANIFEST_NAME, manifest)
        (directory / 'logs').mkdir(exist_ok=True)
        campaign.CAMPAIGN_ROOT = directory.parent.parent
        args = campaign.make_parser().parse_args(['postprocess', '--measurement', measurement['id'], '--tag', manifest['tag']])
        print(f"Reprocessing {measurement['id']} with all saved shards", flush=True)
        campaign.postprocess_campaign(args, measurement)
        regenerated = campaign.load_json(directory / 'postprocess/summary.json')
        compare(summary_science(plan['baseline']), summary_science(regenerated), 'summary')
        compare(plan['baseline']['order_sign_audit']['order_coefficients'], regenerated['order_sign_audit']['order_coefficients'], 'NLO signs')
        compare(plan['baseline']['order_sign_audit']['inputs_sha256'], regenerated['order_sign_audit']['inputs_sha256'], 'NLO input digest')
        yoda = campaign._import_yoda()
        report['yoda_module'] = str(getattr(yoda, '__file__', 'unknown'))
        report['yoda_version'] = str(getattr(yoda, '__version__', 'unknown'))
        nominal_path = directory / post['yoda']
        nominal = yoda.read(str(nominal_path))
        for relative in (post['yoda'], post.get('variations_yoda')):
            if not relative:
                continue
            old = yoda.read(str(directory / 'baseline' / relative))
            new = yoda.read(str(directory / relative))
            compare({p: estimate_science(obj) for p, obj in old.items()},
                    {p: estimate_science(obj) for p, obj in new.items()}, relative, rtol=3e-6, atol=1e-12)
        immutable_postprocess = {str(p.relative_to(directory)): campaign.sha256_file(p)
                                 for p in (directory / 'postprocess').iterdir() if p.is_file()}
        report['postprocessing_matches_baseline'] = True
        campaign.atomic_write_json(report_path, report)
        points = campaign.load_reference_snapshot(measurement)
        modes = [('full', []), ('central', ['--hide-scale-uncertainty'])]
        scales = [float(v[2]) for v in manifest['configuration'].get('variation_points', [[0,0,1.]]) if v[0] == v[1] == 0]
        modes += [(f'mu{campaign._scale_token(scale)}', ['--plot-scale', str(scale)]) for scale in sorted(set(scales))]
        report['plots'] = {}
        for mode, flags in modes:
            print(f"Plotting {measurement['id']}: {mode}", flush=True)
            args = campaign.make_parser().parse_args(['plot', '--measurement', measurement['id'], '--tag', manifest['tag'], *flags])
            output = campaign.plot_campaign(args, measurement)
            selected = nominal
            if mode != 'full':
                selected = yoda.read(str(directory / 'work' / f'plot-{mode}' / nominal_path.name))
                expected = {}
                if mode == 'central':
                    for path, obj in nominal.items():
                        names = {name for b in obj.bins(False, True) for name in b.sources()} - {'hard_scale'}
                        expected[path] = estimate_science(obj, names)
                else:
                    scale = float(flags[-1])
                    source = nominal if scale == 1 else yoda.read(str(directory / post['variations_yoda']))
                    for selection, out in measurement['outputs'].items():
                        if out.get('summary_only'): continue
                        for observable in ('a1','apar','parity_pp_mm','parity_pm_mp'):
                            if observable not in out: continue
                            target = out[observable]
                            source_path = target if scale == 1 else f'/{measurement["id"]}/VARIATIONS/p000-u000-mu{campaign._scale_token(scale)}/{observable}_{selection}'
                            expected[target] = estimate_science(source[source_path], {'stat'})
                compare(expected, {p: estimate_science(obj) for p,obj in selected.items()}, f'{mode} errors', rtol=3e-6)
            for path in paths:
                script = output / (path.lstrip('/') + '.py')
                reference = campaign._fixed_reference_overlay_points(measurement, points, script.stem)
                require(reference is not None, f'No reference for {path}')
                science = estimate_science(selected[path])
                if mode.startswith('mu'):
                    bands = {'monte_carlo': [None if b.get('masked') else
                             sum(abs(v) for v in b['errors']['stat']) / 2 for b in science['bins']]}
                else:
                    bands = dict(campaign._fixed_theory_uncertainty_bands(measurement, regenerated, script.stem))
                    if mode == 'central':
                        bands.pop('scale_down', None); bands.pop('scale_up', None)
                check_plot(script, science, reference, bands)
            report['plots'][mode] = {'directory': str(output), 'checked_a1_plots': len(paths), 'residuals_checked': True}
            campaign.atomic_write_json(report_path, report)
        # Plot-only modes must not alter the full postprocessed results.
        compare(immutable_postprocess, {str(p.relative_to(directory)): campaign.sha256_file(p)
                                      for p in (directory / 'postprocess').iterdir() if p.is_file()},
                'postprocessing preserved during plotting')
        report['status'] = 'passed'
    except BaseException as exc:
        report.update(status='failed', error=f'{type(exc).__name__}: {exc}')
        raise
    finally:
        campaign.CAMPAIGN_ROOT = previous_root
        changed = [relative for relative, source in plan['files'].items()
                   if not source.is_file() or campaign.sha256_file(source) != hashes[relative]]
        report['original_inputs_unchanged'] = not changed
        report['changed_source_files'] = changed
        if changed: report['status'] = 'failed'
        campaign.atomic_write_json(report_path, report)
    require(report['original_inputs_unchanged'], 'Source inputs changed during validation; inspect validation.json')
    return report


def make_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--campaign', type=Path, action='append', required=True,
                        help='Existing completed helium campaign directory; repeat for multiple campaigns')
    parser.add_argument('--output-root', type=Path, required=True,
                        help='New validation directory, separate from sources; must not already exist')
    parser.add_argument('--dry-run', action='store_true', help='Inspect inputs and print copy plan; write nothing, no native imports')
    return parser


def main(argv=None):
    args = make_parser().parse_args(argv)
    output = args.output_root.expanduser().resolve()
    try:
        require(not output.exists(), f'Output already exists; choose a new validation directory: {output}')
        plans = [inspect(source, output) for source in args.campaign]
        require(len({plan['destination'] for plan in plans}) == len(plans), 'Duplicate measurement/tag destination')
        if args.dry_run:
            print(json.dumps([public_plan(plan) for plan in plans], indent=2))
            return 0
        campaign._import_yoda()
        tool = shutil.which('rivet-mkhtml')
        require(tool is not None, 'Activate the Rivet environment: rivet-mkhtml is not on PATH')
        import numpy, matplotlib  # noqa: F401 — fail before copying if rendering dependencies are absent.
        reports = [run_check(plan, Path(tool).resolve()) for plan in plans]
        campaign.atomic_write_json(output / 'validation-summary.json', reports)
        print(f'Passed {len(reports)} saved-campaign checks; report: {output / "validation-summary.json"}')
        return 0
    except (campaign.CampaignError, OSError, ValueError, KeyError, ImportError) as exc:
        print(f'error: {exc}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
