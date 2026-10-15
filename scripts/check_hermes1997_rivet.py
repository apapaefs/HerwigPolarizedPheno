#!/usr/bin/env python3
"""Compile HERMES1997 and check controlled HepMC fixtures; does not run Herwig."""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile

import run_experimental_campaign as campaign
from helium_dis_validation import r1990_fit_b

ROOT = Path(__file__).resolve().parents[1]
ID = "HERMES_1997_I440904"
EDGES = [.023,.04,.055,.075,.1,.14,.2,.3,.4,.6]
MEANS = [.033,.047,.065,.087,.119,.168,.244,.342,.464]
Q2 = [1.22,1.47,1.73,1.99,2.30,2.65,3.07,3.86,5.25]


def check(output: Path) -> dict:
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    environment = dict(os.environ)
    environment["RIVET_ANALYSIS_PATH"] = str(output)+os.pathsep+environment.get("RIVET_ANALYSIS_PATH", "")
    environment["RIVET_DATA_PATH"] = str(ROOT/"analyses/rivet/dis")+os.pathsep+environment.get("RIVET_DATA_PATH", "")

    def run(command, name):
        with (output/(name+".log")).open("w") as log:
            subprocess.run(command, cwd=output, env=environment, stdout=log, stderr=subprocess.STDOUT, check=True)

    flags = lambda option: shlex.split(subprocess.check_output(["rivet-config", option], text=True))
    compiler = shlex.split(os.environ.get("RIVET_CXX", shutil.which("g++") or "c++"))
    include = "-I"+str(ROOT/"analyses/rivet/dis")
    run(compiler+["-std=c++17", "-shared", "-fPIC"]+flags("--cppflags")+
        [include, str(ROOT/f"analyses/rivet/dis/{ID}.cc")]+flags("--libs")+
        ["-o", str(output/f"Rivet{ID}.so")], "build-plugin")
    run(compiler+["-std=c++17"]+flags("--cppflags")+
        [include, str(ROOT/"scripts/tests/hermes1997_fidelity_events.cc")]+flags("--libs")+
        ["-o", str(output/"fixture-events")], "build-fixtures")
    run([str(output/"fixture-events"), str(output)], "create-fixtures")
    run(["rivet", "--show-analysis", ID], "analysis-info")
    info = (output/"analysis-info.log").read_text()
    assert "Parse error" not in info and "ERROR" not in info, info
    assert "27.5" in info and "2112" in info and "e+" in info and "ANY" not in info, info
    report = {}
    for hepmc in sorted(output.glob("*.hepmc")):
        name = hepmc.stem
        command = ["rivet", "-a", ID, "-o", str(output/(name+".yoda")), str(hepmc)]
        if name in {"wrong-lepton", "wrong-target"}:
            command.insert(1, "--ignore-beams")
        run(command, name)
        expected = {
            "nominal-wide": (False, True, True), "wide-only": (False, False, True),
            "outside": (False, False, False), "above-y-cut": (False, False, False),
            "wrong-target": (False, False, False), "below-x-range": (False, False, False),
            "above-x-range": (False, False, False), "below-q2-cut": (False, False, False),
            "wrong-lepton": (False, False, False),
        }.get(name, (True, True, True))
        windows = {}
        for window, accepted in zip(("W025", "W050", "W100"), expected):
            data = campaign.read_histogram_series_many(output/(name+".yoda"), {
                label: f"/{ID}/{hist}_{window}" for label, hist in (
                    ("ordinary", "SigmaX"), ("weighted", "SigmaOverD_X"),
                    ("covariance", "CovarianceProxy_X"), ("xmoment", "XMoment"),
                    ("q2moment", "Q2Moment"))})
            ordinary = data["ordinary"]
            integral = sum(value*(high-low) for value, low, high in zip(ordinary.values, EDGES, EDGES[1:]))
            assert math.isclose(integral, float(accepted), abs_tol=2e-6), (name, window, integral)
            if name.startswith("mean-"):
                _, target, index = name.split("-")
                index = int(index)
                mass = .9382720813 if target == "2212" else .9395654133
                x, q2, energy = MEANS[index], Q2[index], 27.5
                nu = q2/(2*mass*x)
                theta = 2*math.asin(math.sqrt(q2/(4*energy*(energy-nu))))
                epsilon = 1/(1+2*(1+nu*nu/q2)*math.tan(theta/2)**2)
                d = (energy-epsilon*(energy-nu))/(energy*(1+epsilon*r1990_fit_b(x,q2)))
                assert math.isclose(data["weighted"].values[index]/ordinary.values[index], 1/d, rel_tol=3e-6)
                assert math.isclose(data["covariance"].variances[index]/ordinary.variances[index], 1/d, rel_tol=3e-6)
                assert math.isclose(data["xmoment"].values[index]/ordinary.values[index], x, rel_tol=3e-6)
                assert math.isclose(data["q2moment"].values[index]/ordinary.values[index], q2, rel_tol=3e-6)
            windows[window] = integral
        report[name] = windows
    assert len(report) == 20, len(report)
    (output/"results.json").write_text(json.dumps(report, indent=2)+"\n")
    print(f"Passed {len(report)} compiled Rivet fixtures (all 9 bins, three windows, y/beam cuts and 1/D covariance).")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Keep validation fixtures and logs in this directory")
    args = parser.parse_args()
    if args.output:
        check(args.output)
    else:
        with tempfile.TemporaryDirectory(prefix="hermes1997-rivet-check-") as temporary:
            check(Path(temporary))
