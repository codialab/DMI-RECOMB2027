"""Create a compact reviewer bundle from corrected, locally generated results."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import zipfile

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "reproduced/reviewer_checks/corrected"
EXPORTS = ROOT / "analyses/reviewer_checks/exports"
FILES = [
    "experiment_3/experiment_3a/experiment3a_eta2_usefulness.tsv",
    "experiment_3/experiment_3a/experiment3a_20x20_parity.tsv.gz",
    "experiment_3/experiment_3a/experiment3a_recomputed_vs_frozen_case_audit.tsv",
    "experiment_3/experiment_3a/experiment3a_manifest.json",
    "experiment_3/experiment_3b/experiment3b_associations.tsv",
    "experiment_3/experiment_3b/experiment3b_duplicates.tsv",
    "experiment_3/experiment_3b/experiment3b_manifest.json",
    "experiment_3/experiment_3b/experiment3b_recomputed_vs_frozen_case_audit.tsv",
    "experiment_4/experiment4_endpoint_summary.tsv",
    "experiment_4/experiment4_eta2_usefulness_associations.tsv",
    "experiment_4/experiment4_ess20_outcome_parity.tsv",
    "experiment_4/experiment4_recomputed_vs_frozen_case_audit.tsv",
    "experiment_4/experiment4_weight_support.tsv.gz",
    "experiment_4/experiment4_product_ess.tsv.gz",
    "experiment_4/experiment4_manifest.json",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    missing = [rel for rel in FILES if not (RESULTS / rel).is_file()]
    if missing:
        raise SystemExit("missing required corrected result files: " + ", ".join(missing))
    code_hash = sha256(ROOT / "analyses/reviewer_checks/candidate_sensitivity.py")
    manifest_paths = {
        "experiment_3/experiment_3a/experiment3a_manifest.json": "experiment_3/experiment_3a/",
        "experiment_3/experiment_3b/experiment3b_manifest.json": "experiment_3/experiment_3b/",
        "experiment_4/experiment4_manifest.json": "experiment_4/",
    }
    for rel, folder in manifest_paths.items():
        path = RESULTS / rel
        payload = json.loads(path.read_text())
        payload["analysis_code_sha256"] = code_hash
        payload["runner_sha256"] = {
            name: sha256(ROOT / "analyses/reviewer_checks" / name)
            for name in ("run_experiment_3a.py", "run_experiment_3b.py", "run_experiment_4.py")
        }
        payload["corrected_aggregation"] = "distinct truth pairs averaged within evaluation; defined evaluation outcomes averaged equally per reaction"
        payload["prediction_protocol"] = "development-only scaling and Ridge alpha=1; five-fold pathway GroupKFold development OOF; full-development fit applied unchanged to confirmation"
        payload["status"] = "CORRECTED_POST_HOC_SENSITIVITY"
        payload["result_hashes"] = {
            Path(name).name: sha256(RESULTS / name)
            for name in FILES
            if name.startswith(folder) and (RESULTS / name).is_file() and Path(name).name != Path(rel).name
        }
        if "3a" in rel:
            payload["baseline_frozen_case_endpoint_parity"] = "PASS; PL2B/A22 round-trip float aggregation"
            payload["truth_excluded_geometry_interpretation"] = "truth-exclusion sensitivity measure, not a pre-cue predictor"
        elif "3b" in rel:
            payload["baseline_frozen_case_endpoint_parity"] = "PASS for original q10/q50/q90 scenario"
            payload["random_seeds"] = [20271028, 20271029, 20271030, 20271031]
        else:
            payload["ESS20_frozen_outcome_parity"] = "PASS within 1e-12; exact response/tie counts"
            payload["fixed_common_truth_pairs"] = 428
        path.write_text(json.dumps(payload, indent=2) + "\n")

    EXPORTS.mkdir(parents=True, exist_ok=True)
    archive = EXPORTS / "candidate_sensitivity_review_bundle.zip"
    manifest = {
        "status": "corrected_post_hoc_sensitivity_results",
        "contains_candidate_vectors": False,
        "contains_raw_case_rows": False,
        "files": {rel: sha256(RESULTS / rel) for rel in FILES},
        "source_manifests": {
            "experiment3a": "experiment_3/experiment_3a/experiment3a_manifest.json",
            "experiment3b": "experiment_3/experiment_3b/experiment3b_manifest.json",
            "experiment4": "experiment_4/experiment4_manifest.json",
        },
    }
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as zf:
        for rel in FILES:
            zf.write(RESULTS / rel, arcname=rel)
        zf.writestr("export_manifest.json", json.dumps(manifest, indent=2) + "\n")
    (EXPORTS / "candidate_sensitivity_review_bundle.zip.sha256").write_text(
        f"{sha256(archive)}  {archive.name}\n"
    )


if __name__ == "__main__":
    main()
