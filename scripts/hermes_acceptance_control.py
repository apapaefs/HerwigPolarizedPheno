"""Compact shared-event ring/rectangle checks; primary paper galleries stay unique."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path


MEASUREMENT = "HERMES_2007_I726689"
MANIFEST = "cache-manifest.json"
OUTPUTS = {f"AcceptanceControl_{target}.{extension}"
           for target in ("P", "D") for extension in ("png", "pdf")}


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _valid_cache(directory: Path, key: str) -> bool:
    try:
        if directory.is_symlink() or not directory.is_dir():
            return False
        record = json.loads((directory / MANIFEST).read_text(encoding="utf-8"))
        if record["cache_key"] != key or set(record["outputs"]) != OUTPUTS:
            return False
        return all(not (directory / name).is_symlink()
                   and (directory / name).is_file()
                   and (directory / name).stat().st_size > 0
                   and _digest(directory / name) == digest
                   for name, digest in record["outputs"].items())
    except (OSError, ValueError, KeyError, TypeError):
        return False


def _gallery_info(campaign_dir: Path, directory: Path, key: str,
                  source_sha: str, *, created: bool) -> dict:
    return {"directory": directory.relative_to(campaign_dir).as_posix(),
            "cache_key": key, "input_sha256": source_sha, "created": created,
            "difference_uncertainty": "paired same-event covariance",
            "primary": "rectangle_intersect_polar_ring"}


def ensure_gallery(campaign_dir: Path, output_dir: Path) -> dict:
    """Reuse verified plots or create a fresh revision, retaining older outputs."""
    campaign_dir, output_dir = Path(campaign_dir).resolve(), Path(output_dir).resolve()
    source = campaign_dir / "postprocess/summary.json"
    raw, code = source.read_bytes(), Path(__file__).read_bytes()
    summary = json.loads(raw)
    if summary.get("measurement") != MEASUREMENT:
        raise ValueError("Acceptance controls require a HERMES_2007_I726689 summary")
    if summary.get("acceptance") != "rectangle_intersect_polar_ring":
        raise ValueError("Acceptance controls require a v5 rectangle/ring summary")
    key = hashlib.sha256(raw + code).hexdigest()
    source_sha = hashlib.sha256(raw).hexdigest()
    cache_root = output_dir / MEASUREMENT / "acceptance-controls"
    for directory in sorted(cache_root.glob(f"{key}*")):
        if _valid_cache(directory, key):
            if source.read_bytes() != raw:
                raise ValueError("Acceptance-control input changed during cache validation; retry")
            return _gallery_info(campaign_dir, directory, key, source_sha, created=False)
    cache_root.mkdir(parents=True, exist_ok=True)
    attempt = 0
    while True:
        destination = cache_root / (key if attempt == 0 else f"{key}-{attempt:03d}")
        try:
            destination.mkdir()
            break
        except FileExistsError:
            attempt += 1
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = {(row["selection"], row["bin"]): row for row in summary["bins"]}
    for target, label in (("P", "Proton"), ("D", "Deuteron")):
        stem = f"AcceptanceControl_{target}"
        with plt.rc_context({"font.size": 9, "axes.grid": True, "grid.alpha": .2}):
            fig, axes = plt.subplots(2, 2, figsize=(10, 6.5), constrained_layout=True)
            number = 0
            points = []
            for sl in range(1, 20):
                selection = f"Born{target}_X{sl:02d}"
                selected = sorted((v for (s, _), v in rows.items() if s == selection), key=lambda v:v["bin"])
                for primary in selected:
                    number += 1
                    ring = rows[("RingControl_" + selection, primary["bin"])]
                    points.append((number, primary, ring))
            for which, color, title in ((1, "#d55e00", "Rectangle (primary)"), (2, "#0077bb", "Ring control")):
                valid = [(n, p, c) for n,p,c in points if (p,c)[which-1].get("a_parallel") is not None]
                axes[0,0].errorbar([v[0] for v in valid], [v[which]["a_parallel"] for v in valid],
                    yerr=[v[which]["a_parallel_stat"] for v in valid], fmt=".", ms=3, lw=.7, color=color, label=title)
            differences = [(n,c) for n,p,c in points if c.get("ring_minus_primary_a_parallel") is not None]
            axes[1,0].errorbar([n for n,c in differences], [c["ring_minus_primary_a_parallel"] for n,c in differences],
                yerr=[c["ring_minus_primary_a_parallel_stat"] for n,c in differences], fmt=".", color="#333333", lw=.7)
            axes[1,0].axhline(0,color="0.5",lw=.8)
            axes[0,0].set(ylabel=r"Born $A_\parallel$", title="Published cells")
            axes[1,0].set(xlabel="Published cell number", ylabel="Ring − rectangle")
            axes[0,0].legend(fontsize=8)
            selection = "Q2GT1" if target == "P" else "D_Q2GT1"
            for ax, observable, ylabel in ((axes[0,1], "a_parallel", r"Integrated $A_\parallel$"), (axes[1,1], "a1", r"$A_1$ proxy")):
                for prefix,color,title in (("", "#d55e00", "Rectangle"), ("RingControl_", "#0077bb", "Ring")):
                    selected = sorted((v for (s,_),v in rows.items() if s==prefix+selection and v.get(observable) is not None), key=lambda v:v["bin"])
                    ax.errorbar([(v["bin_low"]*v["bin_high"])**.5 for v in selected],
                        [v[observable] for v in selected], yerr=[v[observable+"_stat"] for v in selected],
                        fmt=".-", ms=3, lw=.7, color=color, label=title)
                ax.set(xscale="log", xlabel="x", ylabel=ylabel)
            axes[0,1].set_title(r"$Q^2>1$ GeV$^2$")
            fig.suptitle(label + ": angular acceptance sensitivity (MC statistical errors)")
            for extension in ("pdf", "png"):
                fig.savefig(destination / f"{stem}.{extension}", dpi=180)
            plt.close(fig)
    if source.read_bytes() != raw or Path(__file__).read_bytes() != code:
        raise ValueError("Acceptance-control inputs changed during rendering; retry")
    if any(not (destination / name).is_file() or (destination / name).stat().st_size == 0
           for name in OUTPUTS):
        raise ValueError("Incomplete acceptance-control gallery")
    record = {"schema_version": 1, "cache_key": key, "input_sha256": source_sha,
              "outputs": {name: _digest(destination / name) for name in sorted(OUTPUTS)}}
    (destination / MANIFEST).write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    return _gallery_info(campaign_dir, destination, key, source_sha, created=True)
