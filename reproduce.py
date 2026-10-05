#!/usr/bin/env python3
"""Lightweight publication checks and plotting; never invokes upstream solvers."""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import io
import json
import lzma
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def inside(relative: str | Path) -> Path:
    path = ROOT / relative
    resolved = path.resolve()
    if path.is_symlink() or not resolved.is_relative_to(ROOT):
        raise ValueError(f"Path escapes repository or is a symlink: {relative}")
    return resolved


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_json(relative: str) -> dict:
    return json.loads(inside(relative).read_text())


def input_index() -> dict[str, dict]:
    return {x["destination_path"]: x for x in read_json("manifests/FIGURE_INPUTS.json")["inputs"]}


def check_inputs(paths: list[str]) -> list[dict]:
    index = input_index()
    checks = []
    for rel in paths:
        record = index[rel]
        path = inside(rel)
        if not record["included"] or record["rights_status"] != "APPROVED":
            status = "PENDING_AUTHOR_RIGHTS_REVIEW"
        elif not path.is_file():
            status = "MISSING_INPUT"
        else:
            status = "PASS" if digest(path) == record["sha256"] else "CHECKSUM_MISMATCH"
        checks.append({"path": rel, "status": status})
    return checks


def checksums() -> dict:
    entries = inside("manifests/checksums/SHA256SUMS.txt").read_text().splitlines()
    failures = []
    for line in entries:
        expected, rel = line.split("  ", 1)
        path = inside(rel)
        if not path.is_file() or digest(path) != expected:
            failures.append(rel)
    return {"status": "FAIL" if failures else "PASS", "files_checked": len(entries), "failures": failures,
            "note": "Verifies included release files only; does not establish missing-input coverage."}


def parquet(relative: str):
    import pandas as pd
    with lzma.open(inside(relative), "rb") as stream:
        return pd.read_parquet(io.BytesIO(stream.read()), engine="pyarrow")


def validate_frozen() -> dict:
    index = input_index()
    checks = check_inputs(list(index))
    result = {"status": "PASS", "checks": checks, "scientific_checks": []}
    if any(x["status"] != "PASS" for x in checks):
        result["status"] = "INCOMPLETE"
        return result
    import numpy as np
    import pandas as pd
    from scipy.stats import spearmanr
    points = pd.read_csv(inside("data/figure_inputs/fig2/figure2_DG_points.csv.xz"), float_precision="round_trip")
    assert not points.duplicated(["panel", "reaction_id"]).any()
    expected = {"D": (2454, .8571160663367106), "E": (2485, .8339440056039364),
                "F": (838, .8351435823183291), "G": (843, .8493657186460156)}
    for panel, (count, rho) in expected.items():
        block = points.loc[points.panel.eq(panel)]
        x = block.x_direction_explained_magnitude_variance.to_numpy(float)
        y = block.y_directional_usefulness_fraction.to_numpy(float)
        assert np.isfinite(x).all() and np.isfinite(y).all()
        assert len(block) == count
        assert np.isclose(spearmanr(x, y).statistic, rho, rtol=0, atol=1e-6)
    result["scientific_checks"].append("Figure 2 usefulness panel counts and frozen correlations")
    pairs = parquet("data/figure_inputs/fig3/fig3_gain_reaction_evaluation.parquet.xz")
    counts = pairs.groupby("anchor_setting").size().to_dict()
    assert counts == {"A1": 660237, "A2-L": 660237}  # Historical storage label.
    assert np.isfinite(pairs[["correct_gain_mean", "wrong_gain_mean", "directional_advantage"]]).all().all()
    assert np.allclose(pairs.directional_advantage, pairs.correct_gain_mean - pairs.wrong_gain_mean,
                       rtol=0, atol=5e-12)
    # The check describes existing endpoints; it does not redefine production tolerances.
    result["scientific_checks"].append("Figure 3 pair counts and full correct-minus-wrong endpoint")
    frames = []
    for method in ["CORDA", "GIMME", "iMAT", "RIPTiDe"]:
        rel = f"data/figure_inputs/fig4/fig4_geometry_paired/algorithm={method}/part.parquet.xz"
        frames.append(parquet(rel).assign(algorithm=method))
    geometry = pd.concat(frames, ignore_index=True)
    assert len(geometry) == 1671600 and geometry.evaluation_id.nunique() == 400
    assert geometry.reaction_id.nunique() == 4179
    assert not geometry.duplicated(["evaluation_id", "reaction_id"]).any()
    for method, block in geometry.groupby("algorithm"):
        assert len(block) == 417900
    result["scientific_checks"].append("Figure 4 400 evaluations × 4,179 reactions, one-to-one paired keys")
    return result


def state() -> dict:
    registry = read_json("manifests/FIGURE_REGISTRY.json")
    checks = check_inputs(list(input_index()))
    return {"source_commit": registry["source_commit"], "inputs": checks,
            "figures": [{"id": x["id"], "status": x["status"], "issue": x.get("issue", "")}
                        for x in registry["figures"]]}


def write_new(path: Path, payload: str) -> None:
    path = inside(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as stream:
        stream.write(payload)


def plot(selected: str | None, output_root: str) -> dict:
    registry = read_json("manifests/FIGURE_REGISTRY.json")["figures"]
    names = set(selected.split(",")) if selected else {x["id"] for x in registry}
    if names - {x["id"] for x in registry}:
        raise ValueError("Unknown figure IDs: " + ", ".join(sorted(names - {x["id"] for x in registry})))
    output = inside(output_root)
    if output == ROOT or not output.is_relative_to(ROOT / "reproduced"):
        raise ValueError("Plot outputs must be under reproduced/.")
    jobs = []
    for figure in registry:
        if figure["id"] not in names:
            continue
        job = {"id": figure["id"]}
        if figure["status"] != "READY_WHEN_INPUTS_CLEARED":
            job.update(status="PENDING_AUTHOR_LINEAGE_REVIEW", issue=figure["issue"])
            jobs.append(job)
            continue
        checks = check_inputs(figure["inputs"])
        if any(x["status"] != "PASS" for x in checks):
            job.update(status="PENDING_INPUTS", checks=checks)
            jobs.append(job)
            continue
        path = inside(figure["plotting_script"])
        versions = {name: importlib.metadata.version(name) for name in ["numpy", "pandas", "scipy", "pyarrow", "matplotlib"]}
        identity = {"plotting_sha256": digest(path), "input_sha256": {p: input_index()[p]["sha256"] for p in figure["inputs"]},
                    "versions": versions, "python": sys.version.split()[0]}
        job_dir = inside(output / figure["id"])
        marker = job_dir / "JOB_MANIFEST.json"
        if job_dir.exists():
            old = json.loads(marker.read_text()) if marker.is_file() else {}
            outputs_valid = bool(old.get("output_sha256")) and all(
                inside(job_dir / name).is_file() and digest(inside(job_dir / name)) == sha
                for name, sha in old.get("output_sha256", {}).items())
            job["status"] = "SKIPPED_VALID" if old.get("identity") == identity and outputs_valid else "FAILED_OUTPUT_EXISTS"
            jobs.append(job)
            continue
        job_dir.mkdir(parents=True)
        tempdir = inside(".cache/tmp")
        tempdir.mkdir(parents=True, exist_ok=True)
        mpl = inside(".cache/matplotlib")
        mpl.mkdir(parents=True, exist_ok=True)
        env = dict(os.environ, RECOMB_OUTPUT_DIR=str(job_dir), MPLCONFIGDIR=str(mpl),
                   TMPDIR=str(tempdir), PYTHONDONTWRITEBYTECODE="1", PYTHONPYCACHEPREFIX=str(inside(".cache/pycache")))
        with (job_dir / "plot.log").open("w") as log:
            proc = subprocess.run([sys.executable, "-B", str(path)], cwd=ROOT, env=env, stdout=log, stderr=log)
        hashes = {str(p.relative_to(job_dir)): digest(p) for p in job_dir.iterdir() if p.suffix in [".svg", ".png", ".pdf"]}
        if proc.returncode or not hashes:
            job.update(status="FAILED", log=str((job_dir / "plot.log").relative_to(ROOT)))
        else:
            write_new(marker, json.dumps({"identity": identity, "output_sha256": hashes}, indent=2) + "\n")
            job.update(status="COMPLETED", outputs=list(hashes))
        jobs.append(job)
    result = {"status": "PASS" if all(x["status"] in ["COMPLETED", "SKIPPED_VALID"] for x in jobs) else "INCOMPLETE", "jobs": jobs}
    if jobs:
        output.mkdir(parents=True, exist_ok=True)
        # Append-only run history; existing results and logs are not overwritten.
        with (output / "RUN_HISTORY.jsonl").open("a") as stream:
            stream.write(json.dumps(result, sort_keys=True) + "\n")
    return result


def tables(output_root: str) -> dict:
    paths = ["data/figure_inputs/fig2/figure2_DG_points.csv.xz", "data/figure_inputs/fig3/fig3_gain_summary.parquet.xz",
             "data/figure_inputs/fig3/fig3_directional_advantage_summary.parquet.xz"]
    checks = check_inputs(paths)
    if any(x["status"] != "PASS" for x in checks):
        return {"status": "INCOMPLETE", "checks": checks}
    output = inside(output_root)
    if not output.is_relative_to(ROOT / "reproduced"):
        raise ValueError("Derived tables must be under reproduced/.")
    import pandas as pd
    from scipy.stats import spearmanr
    points = pd.read_csv(inside(paths[0]), float_precision="round_trip")
    rows = []
    for panel, block in points.groupby("panel", sort=True):
        rows.append({"panel": panel, "n_plotted_reactions": len(block),
                     "spearman_rho": float(spearmanr(block.x_direction_explained_magnitude_variance,
                                                     block.y_directional_usefulness_fraction).statistic)})
    derived = {"fig2_panel_statistics.tsv": pd.DataFrame(rows),
               "fig3_pair_gain_summary.tsv": parquet(paths[1]),
               "fig3_full_endpoint_separation.tsv": parquet(paths[2])}
    if any(inside(output / name).exists() for name in derived):
        raise FileExistsError("Derived table exists; choose a new reproduced/ output directory.")
    for name, frame in derived.items():
        write_new(output / name, frame.to_csv(sep="\t", index=False))
    return {"status": "PASS", "outputs": [str((output / name).relative_to(ROOT)) for name in derived]}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ["status", "checksums", "validate"]:
        sub.add_parser(name)
    p = sub.add_parser("figures")
    p.add_argument("--only", help="Comma-separated IDs from status")
    p.add_argument("--output-root", default="reproduced/figures")
    p = sub.add_parser("tables")
    p.add_argument("--output-root", default="reproduced/tables")
    args = parser.parse_args()
    try:
        boundary = Path(subprocess.check_output(
            ["git", "rev-parse", "--show-toplevel"], cwd=ROOT, text=True,
            env=dict(os.environ, GIT_OPTIONAL_LOCKS="0")).strip()).resolve()
        if boundary != ROOT:
            raise ValueError("Run from a Git checkout whose root contains reproduce.py.")
        if args.command == "status": result = state()
        elif args.command == "checksums": result = checksums()
        elif args.command == "validate": result = validate_frozen()
        elif args.command == "figures": result = plot(args.only, args.output_root)
        else: result = tables(args.output_root)
    except (AssertionError, OSError, ValueError, KeyError, subprocess.CalledProcessError) as exc:
        print(json.dumps({"status": "FAIL", "error": str(exc)}, indent=2))
        return 1
    print(json.dumps(result, indent=2))
    return 0 if args.command == "status" or result["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
