#!/usr/bin/env python3
"""Run controlled events through the compiled Rivet plugins (Herwig env needed)."""
import argparse
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile

import numpy as np
import run_experimental_campaign as experimental

ROOT = Path(__file__).resolve().parents[1]
MULTIPLICITIES = ("COMPASS_2017_I1444985", "COMPASS_2017_I1483098",
                  "COMPASS_2025_I2840545", "COMPASS_2026_I3096394")


def check(output):
    output.mkdir(parents=True, exist_ok=True)
    environment = dict(os.environ)
    for name, paths in {
        "RIVET_ANALYSIS_PATH": (ROOT/"build", ROOT/"analyses/rivet/dis", ROOT/"analyses/rivet/pp"),
        "RIVET_DATA_PATH": (ROOT/"analyses/rivet/dis", ROOT/"analyses/rivet/pp"),
    }.items():
        environment[name] = os.pathsep.join(map(str, paths))+os.pathsep+environment.get(name, "")

    def run(command, name):
        with (output/(name+".log")).open("w") as log:
            subprocess.run(command, cwd=ROOT, env=environment, stdout=log,
                           stderr=subprocess.STDOUT, check=True)

    # Match the Makefile's compiler fallback; an inherited generic CXX may
    # use libc++ while this Rivet installation was built with GCC/libstdc++.
    configured = subprocess.check_output(["rivet-config", "--cxx"], text=True).strip()
    executable = shlex.split(configured)[0]
    compiler = os.environ.get("RIVET_CXX") or (
        configured if Path(executable).is_file() and os.access(executable, os.X_OK)
        else shutil.which("g++-16") or shutil.which("g++-15") or shutil.which("g++") or "c++")
    flags = lambda option: shlex.split(subprocess.check_output(["rivet-config", option], text=True))
    run(shlex.split(compiler)+["-std=c++17"]+flags("--cppflags")+
        ["-I"+str(ROOT/"analyses/rivet/dis"), str(ROOT/"scripts/tests/rivet_fidelity_events.cc")]+
        flags("--libs")+["-o", str(output/"generate")], "compile")
    run([str(output/"generate"), str(output)], "generate")
    results = {}

    def rivet(event, analyses):
        destination = output/(event+".yoda")
        command = ["rivet", "--ignore-beams", "-o", str(destination)]
        for analysis in analyses:
            command += ["-a", analysis]
        run(command+[str(output/(event+".hepmc"))], event)
        return destination

    for energy, analysis, observable in (
        (200, "STAR_2021_I1850855", "inclusive_combined_Yield"),
        (510, "STAR_2022_I1949588", "inclusive_Yield"),
    ):
        for configuration in ("three", "outside"):
            event = f"star-{energy}-{configuration}"
            path = rivet(event, [analysis+":LEVEL=HADRON"])
            item = experimental.read_histogram_series(path, f"/{analysis}:LEVEL=HADRON/{observable}")
            integral = float(np.dot(item.values, np.diff(item.edges)))
            np.testing.assert_allclose(integral, 2., atol=1.e-6)
            results[event] = {"inclusive_jet_yield_pb": integral}

    for kind, species in (("pion", "piplus"), ("kaon", "kplus")):
        for selection, expected in (("reject", 0.), ("pass", 1.)):
            event = f"compass-{selection}-{kind}"
            path = rivet(event, MULTIPLICITIES)
            results[event] = {}
            for analysis in MULTIPLICITIES:
                if kind == "pion" and analysis == "COMPASS_2017_I1483098":
                    continue
                if kind == "kaon" and analysis == "COMPASS_2017_I1444985":
                    continue
                snapshot = json.loads((ROOT/"data/phenomenology"/analysis/"reference.json").read_text())
                dataset = snapshot["datasets"][species]
                # Select the actual released cell containing this event.
                cell = next(i for i, p in enumerate(dataset["points"])
                            if p["x_low"] <= .0098 < p["x_high"]
                            and p["y_low"] <= (.35 if selection == "reject" else .45) < p["y_high"]
                            and p["z_low"] <= .23 < p["z_high"])
                series = experimental.read_histogram_series_many(path, {
                    key: f"/{analysis}/{dataset['raw_objects'][key]}"
                    for key in ("numerator", "denominator")})
                for key, item in series.items():
                    np.testing.assert_allclose(item.values[cell], expected, atol=1.e-12)
                    results[event][analysis+"/"+key] = item.values[cell]

    path = rivet("compass-covariance", ["COMPASS_2013_I1236358", "COMPASS_2014_I1278730"])
    # The two positive pions occupy different bins in the same cell.
    for analysis, inputs, cov, dimensions, cell in (
        ("COMPASS_2013_I1236358", "HadronYield_hplus_pt2_cells", "SpectrumCovariance_hplus_cells", 36, 24),
        ("COMPASS_2014_I1278730", "AzimuthalInputs_hplus_cos1_x", "AzimuthalCovariance_hplus_cos1_x", 17, 1),
    ):
        series = experimental.read_histogram_series_many(path, {
            "inputs": f"/{analysis}/{inputs}", "covariance": f"/{analysis}/{cov}"})
        counts = np.asarray(series["inputs"].values[cell*dimensions:(cell+1)*dimensions])
        diagonal = np.asarray(series["inputs"].variances[cell*dimensions:(cell+1)*dimensions])
        bins = np.flatnonzero(counts[:16] if dimensions == 17 else counts)
        assert len(bins) == 2, (analysis, counts)
        np.testing.assert_allclose(counts[bins], [1., 1.])
        # Rivet's default text output rounds values/errors independently.
        np.testing.assert_allclose(diagonal, counts**2, rtol=2.e-6, atol=1.e-8)
        pairs = dimensions*(dimensions-1)//2
        observed = series["covariance"].variances[cell*pairs:(cell+1)*pairs]
        expected = np.outer(counts, counts)[np.triu_indices(dimensions, 1)]
        np.testing.assert_allclose(observed, expected, rtol=2.e-6, atol=1.e-8)
        results[analysis+"/event_covariance"] = "all diagonal/off-diagonal entries close"
    (output/"results.json").write_text(json.dumps(results, indent=2)+"\n")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Retain fixture events, logs and YODA in this directory")
    args = parser.parse_args()
    if args.output:
        check(args.output.resolve())
    else:
        with tempfile.TemporaryDirectory(prefix="rivet-fidelity-") as temporary:
            check(Path(temporary))
