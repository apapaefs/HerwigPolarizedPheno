#!/usr/bin/env python3
"""Pinned COMPASS 1211.6849v2 Tables 5–7 and data-only plot export.

No generator, helicity estimator or detector response is supplied by this
module. The separate reference registry must not be dispatched to a DIS run.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]
ID = "COMPASS_2013_I1204782"
DATA = Path("data/reference") / ID
PDF_URL = "https://arxiv.org/pdf/1211.6849v2"
PDF_SHA256 = "2b876dc320d93bdd2815c443fc0bd595eac80655c62de375e9d1a069b12e0acc"
CSV_SHA256 = "7d3cc2c53537f7ef28439bcffde4e6aebbda1574dec14004d4a7c848d2b5e18c"
PT_BINS = [(0.0, 0.3), (0.3, 0.7), (0.7, 1.0), (1.0, 1.5), (1.5, None)]
E_BINS = [(0.0, 30.0), (30.0, 50.0), (50.0, None)]
PT_LABELS = ["0–0.3", "0.3–0.7", "0.7–1", "1–1.5", "> 1.5"]
E_LABELS = ["0–30 GeV", "30–50 GeV", "> 50 GeV"]
SAMPLES = {
    5: ["D0_Kpi", "Dstar_Kpi", "Dstar_Ksubpi"],
    6: ["Dstar_Kpipi0"],
    7: ["Dstar_Kpipipi"],
}
SAMPLE_LABELS = {
    5: "Combined untagged / tagged Kπ samples",
    6: "D* tagged Kππ⁰ sample (π⁰ unobserved)",
    7: "D* tagged Kπππ sample",
}
COLUMNS = ("table,pt_low_GeV,pt_high_GeV,E_low_GeV,E_high_GeV,A_gammaN,"
           "stat,syst,mean_y,mean_Q2_GeV2,mean_pt_GeV,mean_E_GeV,mean_D").split(",")


def encoded(value: object) -> bytes:
    return (json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n").encode()


def source_file(repo: Path, relative: str | Path) -> Path:
    path = repo / relative
    if (not path.is_file() or path.is_symlink()
            or not path.resolve().is_relative_to(repo.resolve())):
        raise ValueError(f"Missing or unsafe open-charm source: {path}")
    return path


def descriptor(repo: Path = ROOT) -> dict:
    path = source_file(repo, Path("config/reference") / f"{ID}.json")
    value = json.loads(path.read_text())
    if (value.get("schema_version") != 1 or value.get("id") != ID
            or value.get("process_kind") != "polarized_open_charm"
            or value.get("analysis", {}).get("status") != "REFERENCE_ONLY"
            or value.get("simulation", {}).get("enabled") is not False
            or value.get("simulation", {}).get("status") != "not_implemented"
            or not value["simulation"].get("reason")
            or any(key in value for key in ("cards", "families", "postprocessor"))):
        raise ValueError("Open-charm registration must remain explicitly reference-only")
    for key, filename in {"snapshot": "reference.json", "source_manifest": "source-manifest.json",
                          "raw_snapshot": "tables-5-7.csv", "publication_pdf": "arXiv-1211.6849v2.pdf"}.items():
        if value.get("reference", {}).get(key) != str(DATA / filename):
            raise ValueError(f"Unexpected open-charm reference path: {key}")
    return value


def discover(repo: Path = ROOT) -> dict:
    paths = sorted((repo / "config/reference").glob("*.json"))
    result = {}
    for path in paths:
        if path.stem != ID:
            raise ValueError(f"Unsupported reference-only measurement: {path.stem}")
        result[ID] = descriptor(repo)
    return result


def parse_csv(payload: bytes) -> list[dict]:
    reader = csv.DictReader(io.StringIO(payload.decode()))
    if reader.fieldnames != COLUMNS:
        raise ValueError("Unexpected open-charm table columns")
    rows = []
    for row in reader:
        if set(row) != set(COLUMNS) or any(value is None for value in row.values()):
            raise ValueError("Malformed open-charm table row")
        point = {key: float(row[key]) if row[key] else None for key in COLUMNS}
        if any(value is not None and not math.isfinite(value) for value in point.values()):
            raise ValueError("Non-finite open-charm table entry")
        if any(point[key] is None for key in COLUMNS if key not in {"pt_high_GeV", "E_high_GeV"}):
            raise ValueError("Missing open-charm table entry")
        if point["stat"] <= 0 or point["syst"] < 0 or not 0 < point["mean_D"] <= 1:
            raise ValueError("Invalid open-charm uncertainty or depolarisation factor")
        rows.append(point)
    expected = [(table, *pt, *energy) for table in SAMPLES for pt in PT_BINS for energy in E_BINS]
    actual = [tuple(p[key] for key in COLUMNS[:5]) for p in rows]
    if actual != expected:
        raise ValueError("Expected all 45 distinct cells in published table order")
    for point in rows:
        for variable in ("pt", "E"):
            low, high, mean = (point[f"{variable}_low_GeV"], point[f"{variable}_high_GeV"], point[f"mean_{variable}_GeV"])
            if mean < low or (high is not None and mean > high):
                raise ValueError("Published mean lies outside its bin")
    return rows


def normalized(rows: list[dict]) -> dict:
    datasets = []
    for table, samples in SAMPLES.items():
        for energy_index, (elow, ehigh) in enumerate(E_BINS):
            points = []
            for point in rows:
                if point["table"] != table or (point["E_low_GeV"], point["E_high_GeV"]) != (elow, ehigh):
                    continue
                points.append({
                    "bin_index": len(points) + 1,
                    "pt_bin_GeV": [point["pt_low_GeV"], point["pt_high_GeV"]],
                    "value": point["A_gammaN"],
                    "stat": point["stat"], "syst": point["syst"],
                    "total": math.hypot(point["stat"], point["syst"]),
                    "means": {key.removeprefix("mean_"): point[key] for key in COLUMNS if key.startswith("mean_")},
                })
            datasets.append({"id": f"table{table}-energy{energy_index+1}",
                             "table": table, "samples": samples, "label": SAMPLE_LABELS[table],
                             "energy_bin_GeV": [elow, ehigh], "energy_label": E_LABELS[energy_index],
                             "points": points})
    return {
        "measurement": ID, "schema_version": 1, "observable": "A_gammaN_to_D0X",
        "provenance": {
            "source": "COMPASS publication, Tables 5–7",
            "publication": "https://arxiv.org/abs/1211.6849v2",
            "publication_pdf_url": PDF_URL,
            "doi": "10.1103/PhysRevD.87.052018", "inspire_id": 1204782,
            "notes": "Exact numerical transcription of arXiv v2 Tables 5–7, visually checked against PDF pages 21, 22 and 25. Experimental data only; no Herwig prediction or extracted Delta g/g fit is included.",
        },
        "binning": {"pt_bins_GeV": PT_BINS, "energy_bins_GeV": E_BINS,
                    "open_upper_bound": "null means open-ended, not an invented finite edge",
                    "plot_axis": "Five categorical published pT intervals; spacings are not physical bin widths"},
        "uncertainties": {"stat": "First published uncertainty", "syst": "Second published uncertainty",
                          "total": "Quadrature of published statistical and systematic errors, for display only",
                          "covariance": "Not supplied in Tables 5–7; no diagonal-fit assumption or chi2 is made",
                          "range": "Measured estimates and errors are retained even outside [-1, 1]"},
        "estimator": {
            "observable": "Photon–nucleon asymmetry, A_muN / D, not the SIDIS A1 estimator",
            "D_equation_2": "y*(2-y-2*y*y*m_mu^2/Q^2)/(y*y*(1-2*m_mu^2/Q^2)+2*(1-y))",
            "signal_weight": "w_S = P_mu*f*D*s/(s+b)",
            "published_means": "Means weighted by w_S^2; all five means per cell, including D, are retained exactly",
            "muon_asymmetry_conversion": "Section 4.2 permits conversion of the published A_gammaN by the tabulated mean D; no conversion is applied here",
            "response": "Signal/background separation, acceptance and signal-purity weighting cannot be reconstructed from Tables 5–7 alone",
        },
        "selection": {
            "scope": "Published reconstruction cuts (Table 4 and Section 3), recorded as documentation; not implemented particle-level cuts",
            "beam": "Positive 160 GeV/c muons; combined proton and deuteron target data from 2002–2007",
            "common": "Incoming and outgoing muon, at least two additional charged tracks, incoming trajectory crossing the target and vertex inside the target; RICH identification",
            "channels": {
                "D0_Kpi": {"mass_offset_MeV": [-400, 400], "abs_cos_theta_star_lt": 0.65, "z": [0.20, 0.85], "pK_GeV": [9.5, 50], "pPi_GeV": [7, 50], "tag": "Exclude reconstructed D* tags"},
                "Dstar_Kpi": {"mass_offset_MeV": [-600, 600], "abs_cos_theta_star_lt": 0.90, "z": [0.20, 0.85], "pK_GeV": [9.5, 50], "pPi_GeV": [2.5, 50], "deltaM_MeV": [3.2, 8.9]},
                "Dstar_Kpipi0": {"mass_offset_MeV": [-600, 600], "abs_cos_theta_star_lt": 0.90, "z": [0.20, 0.85], "pK_GeV": [9.5, 50], "pPi_GeV": [2.5, 50], "deltaM_MeV": [3.2, 8.9], "reconstruction": "Neutral pion is not reconstructed; additional neural-network selection"},
                "Dstar_Ksubpi": {"mass_offset_MeV": [-400, 400], "abs_cos_theta_star_lt": 0.85, "z": [0.25, 0.85], "pK_GeV": [2.5, 9.5], "pPi_GeV": [2.5, 50], "deltaM_MeV": [3.2, 8.9]},
                "Dstar_Kpipipi": {"mass_offset_MeV": [-400, 400], "abs_cos_theta_star_lt": 0.85, "z": [0.30, 0.85], "pK_GeV": [9.5, 50], "pPi_GeV": [2.5, 50], "deltaM_MeV": [4.0, 7.5]},
            },
            "tagged_common": "Slow pion momentum < 8 GeV/c and electron rejection; DeltaM is reconstructed D* mass minus reconstructed D0-candidate mass minus charged-pion mass",
            "candidate_priority": "Dstar_Kpipipi, then Dstar_Kpi or Dstar_Kpipi0, then D0_Kpi, then Dstar_Ksubpi; one of two candidates in the same sample is chosen randomly",
            "charge_conjugates": True,
            "observed_kinematics_not_fiducial_cuts": {"mean_Q2_GeV2": 0.6, "Q2_range_GeV2_approx": [0.001, 30], "y_range_approx": [0.1, 1], "mean_y_approx": 0.63},
        },
        "datasets": datasets,
    }


def expected(repo: Path) -> tuple[dict, dict]:
    raw = source_file(repo, DATA / "tables-5-7.csv").read_bytes()
    pdf = source_file(repo, DATA / "arXiv-1211.6849v2.pdf").read_bytes()
    if hashlib.sha256(raw).hexdigest() != CSV_SHA256 or hashlib.sha256(pdf).hexdigest() != PDF_SHA256:
        raise ValueError("Pinned open-charm CSV/PDF checksum mismatch")
    snapshot = normalized(parse_csv(raw))
    manifest = {
        "measurement": ID, "source_version": "arXiv:1211.6849v2", "retrieved_at": "2026-09-14",
        "sources": [{"path": str(DATA / "arXiv-1211.6849v2.pdf"), "url": PDF_URL, "sha256": PDF_SHA256},
                    {"path": str(DATA / "tables-5-7.csv"), "kind": "manual numerical transcription", "sha256": CSV_SHA256}],
        "snapshot": {"path": str(DATA / "reference.json"), "sha256": hashlib.sha256(encoded(snapshot)).hexdigest()},
        "tables": [5, 6, 7], "pdf_pages_one_based": [21, 22, 25], "printed_pages": [20, 21, 24],
        "points": 45, "verification": "PDF table images plus independent pdftotext row comparison; use check --verify-pdf-text to repeat text comparison",
    }
    return snapshot, manifest


def validate(repo: Path = ROOT, *, verify_pdf_text: bool = False) -> dict:
    descriptor(repo)
    snapshot, manifest = expected(repo)
    for name, expected_value in [("reference.json", snapshot), ("source-manifest.json", manifest)]:
        if source_file(repo, DATA / name).read_bytes() != encoded(expected_value):
            raise ValueError(f"Open-charm {name} differs from the pinned definition; run rebuild after review")
    if verify_pdf_text:
        rows = parse_csv(source_file(repo, DATA / "tables-5-7.csv").read_bytes())
        number = r"([+−-]?\d+\.\d+)"
        # Rotated running section numbers 4/5 can share a text line with a row.
        pattern = re.compile(r"^[ \t]*(?:[45][ \t]{5,})?(0–0\.3|0\.3–0\.7|0\.7–1|1–1\.5|>\s*1\.5)\s+"
                             r"(0–30|30–50|>\s*50)\s+" + number + r"\s*±\s*" + number
                             + r"\s*±\s*" + number + (r"\s+" + number) * 5 + r"\s*$", re.M)
        for table, page in [(5, 21), (6, 22), (7, 25)]:
            text = subprocess.check_output(["pdftotext", "-f", str(page), "-l", str(page), "-layout",
                                            str(repo / DATA / "arXiv-1211.6849v2.pdf"), "-"], text=True)
            found = pattern.findall(text)
            wanted = [row for row in rows if row["table"] == table]
            if len(found) != 15:
                raise ValueError(f"Table {table}: expected 15 PDF rows, got {len(found)}")
            for i, (extracted, point) in enumerate(zip(found, wanted)):
                pt, energy, *values = extracted
                if (re.sub(r"\s+", "", pt) != PT_LABELS[i//3].replace(" ", "")
                        or re.sub(r"\s+", "", energy) != E_LABELS[i%3].replace(" GeV", "").replace(" ", "")
                        or [float(v.replace("−", "-")) for v in values] != [point[key] for key in COLUMNS[5:]]):
                    raise ValueError(f"CSV disagrees with PDF Table {table}, row {i+1}")
    return snapshot


def signature(repo: Path, value: dict) -> str:
    validate(repo)
    digest = hashlib.sha256(encoded({k: v for k, v in value.items() if not k.startswith("_")}))
    for key in ("snapshot", "source_manifest", "raw_snapshot", "publication_pdf"):
        relative = value["reference"][key]
        digest.update(relative.encode())
        digest.update(source_file(repo, relative).read_bytes())
    return digest.hexdigest()


def plot_records(snapshot: dict, formats: str = "both") -> list[dict]:
    return [{"name": dataset["id"],
             "labels": {"Title": f"Table {dataset['table']} · {dataset['label']} · E(D0) {dataset['energy_label']} · data only"},
             "files": {fmt: f"analyses/{ID}/plots/{dataset['id']}.{fmt}"
                       for fmt in ("png", "pdf") if formats in ("both", fmt)}}
            for dataset in snapshot["datasets"]]


def render(snapshot: dict, output: Path, formats: str = "both") -> None:
    if formats not in {"both", "png", "pdf"}:
        raise ValueError("Reference plot format must be both, png or pdf")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 11, "text.usetex": False,
                         "axes.spines.top": False, "axes.spines.right": False}):
        for dataset, record in zip(snapshot["datasets"], plot_records(snapshot, formats)):
            fig, ax = plt.subplots(figsize=(8.2, 5.7))
            fig.subplots_adjust(left=0.12, right=0.97, bottom=0.28, top=0.78)
            points = dataset["points"]
            x = list(range(1, 6))
            y = [point["value"] for point in points]
            ax.axhline(0, color="#a7adb3", linewidth=0.8, zorder=0)
            ax.errorbar(x, y, yerr=[p["total"] for p in points], fmt="none", ecolor="#a7bac8",
                        elinewidth=4, capsize=0, label="Stat. ⊕ syst.", zorder=2)
            ax.errorbar(x, y, yerr=[p["stat"] for p in points], fmt="o", color="#164c70",
                        elinewidth=1.2, markersize=5, capsize=5, label="Data ± stat.", zorder=3)
            ax.set_xticks(x, PT_LABELS)
            ax.set_xlim(0.5, 5.5)
            ax.margins(y=0.15)
            ax.set_xlabel(r"Published $p_T^{D^0}$ interval [GeV/$c$] (categorical)", labelpad=10)
            ax.set_ylabel(r"$A^{\gamma N\to D^0 X}$")
            ax.grid(axis="y", color="#e4e9ed", linewidth=0.6)
            fig.text(0.12, 0.93, "COMPASS · Open-charm spin asymmetries", fontsize=15, weight="bold")
            fig.text(0.12, 0.875, f"Table {dataset['table']}  |  {dataset['label']}", fontsize=11)
            fig.text(0.12, 0.83, r"$E_{D^0}$: " + dataset["energy_label"] + "  ·  160 GeV/c muons", fontsize=11)
            fig.legend(*ax.get_legend_handles_labels(), loc="lower center", bbox_to_anchor=(0.55, 0.08),
                       frameon=False, ncol=2, fontsize=10)
            fig.text(0.12, 0.062, "REFERENCE DATA ONLY · No Herwig prediction", fontsize=9, color="#6f4c25", weight="bold")
            fig.text(0.12, 0.025, "arXiv:1211.6849v2 · Published errors retained; no clipping or fit", fontsize=9, color="#58636c")
            for fmt, relative in record["files"].items():
                destination = output / relative
                destination.parent.mkdir(parents=True, exist_ok=True)
                if destination.exists():
                    raise ValueError(f"Refusing to overwrite reference plot: {destination}")
                fig.savefig(destination, dpi=180, metadata={"Title": record["labels"]["Title"]})
            plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("check", "rebuild", "plot"))
    parser.add_argument("--repo", type=Path, default=ROOT)
    parser.add_argument("--verify-pdf-text", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    try:
        if args.command == "rebuild":
            descriptor(args.repo)
            snapshot, manifest = expected(args.repo)
            for name, value in [("reference.json", snapshot), ("source-manifest.json", manifest)]:
                (args.repo / DATA / name).write_bytes(encoded(value))
        snapshot = validate(args.repo, verify_pdf_text=args.verify_pdf_text)
        if args.command == "plot":
            if args.output is None or args.output.exists():
                raise ValueError("plot requires --output pointing to a new directory")
            render(snapshot, args.output)
        print(f"{ID}: validated 45 published points in 9 panels; reference data only")
        return 0
    except (OSError, ValueError, subprocess.SubprocessError) as exc:
        parser.exit(2, f"error: {exc}\n")


if __name__ == "__main__":
    raise SystemExit(main())
