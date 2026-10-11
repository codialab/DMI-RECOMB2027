import csv
import io
import json
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import dmi_bridge_uncertainty_reduction_v1 as runner


def test_feature_means_derives_coverage_from_polarity_masses(tmp_path):
    path = tmp_path / "features.tsv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["reaction_id", "sign_magnitude_eta2", "directional_entropy3", "dominant_sign_mass", "n_supported_sign_states", "p_negative", "p_positive"], delimiter="\t")
        writer.writeheader()
        writer.writerow(dict(reaction_id="R", sign_magnitude_eta2="0.5", directional_entropy3="0.2", dominant_sign_mass="0.8", n_supported_sign_states="2", p_negative="0.3", p_positive="0.6"))
        writer.writerow(dict(reaction_id="R", sign_magnitude_eta2="0.7", directional_entropy3="0.4", dominant_sign_mass="0.6", n_supported_sign_states="3", p_negative="0.2", p_positive="0.5"))
    means = runner._feature_means([path], ["R"], ("sign_magnitude_eta2", "directional_entropy3", "dominant_sign_mass", "n_supported_sign_states", "p_negative", "p_positive"))
    assert means["R"]["non_tie_coverage"] == pytest.approx(0.8)
    assert means["R"]["sign_magnitude_eta2"] == pytest.approx(0.6)


def test_support_cache_position_accepts_frozen_a1_and_a2_schemas():
    index = {("alg", "CT2A", "ens"): 7}
    assert runner._support_cache_position({"cache_index": 3}, index) == 3
    assert runner._support_cache_position({"cache_key": ("alg", "CT2A", "ens")}, index) == 7


def test_accumulator_treats_csv_nan_as_missing():
    acc = runner.CaseAccumulator(["E"], ["R"])
    row = {"truth_status": "NON_TIE", "relative_variance_eligible": "true", "directionally_useful": "false"}
    for cue in ("correct", "incorrect"):
        row[f"prediction_gain_{cue}"] = "nan"
        for lam in (0.0, 0.25, 0.5, 1.0):
            row[f"variance_change_status_{cue}_lambda_{lam:g}"] = "NEGLIGIBLE"
            row[f"prediction_gain_{cue}_lambda_{lam:g}"] = "nan"
    for field in ("relative_variance_reduction_correct_lambda_0.25", "relative_variance_reduction_incorrect_lambda_0.25", "absolute_variance_reduction_correct_lambda_0.25", "absolute_variance_reduction_incorrect_lambda_0.25", "prediction_gain_correct", "prediction_gain_incorrect", "crps_improvement_correct_lambda_0.25", "crps_improvement_incorrect_lambda_0.25", "baseline_sd", "baseline_iqr", "baseline_central90_width", "sd_correct_lambda_0.25", "sd_incorrect_lambda_0.25", "iqr_correct_lambda_0.25", "iqr_incorrect_lambda_0.25", "central90_width_correct_lambda_0.25", "central90_width_incorrect_lambda_0.25"):
        row[field] = "nan"
    acc.add({"evaluation_id": "E", "algorithm": "M", "rna_context_key": "RNA"}, "R", row)
    assert not np.isnan(acc.sums["gain_correct"][0, 0])
    assert acc.counts["gain_correct"][0, 0] == 0


def test_sensitivity_stratum_uses_matching_lambda_gain_category():
    acc = runner.CaseAccumulator(["E"], ["R"])
    row = {"truth_status": "NON_TIE", "relative_variance_eligible": True, "directionally_useful": False}
    for cue in ("correct", "incorrect"):
        row[f"prediction_gain_{cue}"] = 1.0
        for lam in (0.0, 0.25, 0.5, 1.0):
            row[f"variance_change_status_{cue}_lambda_{lam:g}"] = "DECREASED"
            row[f"prediction_gain_{cue}_lambda_{lam:g}"] = -1.0 if lam == 0.5 else 1.0
    for field in ("relative_variance_reduction_correct_lambda_0.25", "relative_variance_reduction_incorrect_lambda_0.25", "absolute_variance_reduction_correct_lambda_0.25", "absolute_variance_reduction_incorrect_lambda_0.25", "prediction_gain_correct", "prediction_gain_incorrect", "crps_improvement_correct_lambda_0.25", "crps_improvement_incorrect_lambda_0.25", "baseline_sd", "baseline_iqr", "baseline_central90_width", "sd_correct_lambda_0.25", "sd_incorrect_lambda_0.25", "iqr_correct_lambda_0.25", "iqr_incorrect_lambda_0.25", "central90_width_correct_lambda_0.25", "central90_width_incorrect_lambda_0.25"):
        row[field] = 0.0
    acc.add({"evaluation_id": "E", "algorithm": "M", "rna_context_key": "RNA"}, "R", row)
    assert acc.stratum_summaries[("M", "RNA", "correct", "lambda_0.5", "DECREASED", "WORSENED")] == 1


def test_direct_runner_imports_and_run_gate_without_external_pythonpath():
    import os
    import subprocess
    env = dict(os.environ)
    env.pop("PYTHONPATH", None)
    result = subprocess.run([sys.executable, "scripts/dmi_bridge_uncertainty_reduction_v1.py", "--run"], cwd=ROOT, env=env, capture_output=True, text=True)
    assert "ModuleNotFoundError" not in result.stderr
    assert "frozen input" in result.stderr.lower() or "BLOCKED_MISSING_FROZEN_INPUTS" in result.stdout or result.returncode == 0


def test_crossed_bootstrap_reports_paired_a2_minus_a1(tmp_path):
    pd = pytest.importorskip("pandas")
    paths = {}
    for anchor, shift in (("A1", 0.0), ("A2", 2.0)):
        rows = []
        for i in range(5):
            for j in range(5):
                row = {"ct2a_mouse": f"c{i}", "gl261_mouse": f"g{j}"}
                for column in ("relative_correct", "relative_incorrect", "absolute_correct", "absolute_incorrect", "gain_correct", "gain_incorrect", "crps_correct", "crps_incorrect"):
                    row[column] = shift + i + j
                rows.append(row)
        path = tmp_path / f"{anchor}.tsv.xz"
        pd.DataFrame(rows).to_csv(path, sep="\t", index=False, compression="xz")
        paths[anchor] = path
    records = runner._mouse_cluster_bootstrap(paths, replicates=20, seed=9)
    paired = next(r for r in records if r["anchor"] == "A2_minus_A1" and r["outcome"] == "relative_correct")
    assert paired["mean"] == pytest.approx(2.0)
    assert paired["ci95_low"] == pytest.approx(2.0)
    assert paired["ci95_high"] == pytest.approx(2.0)


def test_analyze_anchor_interrupted_resume_matches_fresh_run(tmp_path, monkeypatch):
    from scripts import dmi_bridge_pl2_sign_only_core_v1 as pl2
    from scripts import dmi_bridge_pl2c_geometry_utility_core_v1 as pl2c
    from scripts import dmi_bridge_uncertainty_reduction_core_v1 as core

    mats = np.array([[[0.0], [1.0], [2.0], [3.0]], [[0.0], [0.5], [1.0], [1.5]]])
    specs, cache_rows = [], []
    for i in range(2):
        pair = {"truth_pair_id": f"P{i}", "evaluation_id": "E", "algorithm": "M", "rna_context_key": "RNA", "ct2a_mouse": "c0", "gl261_mouse": "g0"}
        specs.append({"pair": pair,
                      "ct2a": {"cache_index": 0, "held_samples": np.array([0, 1]), "held_weights": np.array([0.5, 0.5]), "truth_sample": 2 + i},
                      "gl261": {"cache_index": 1, "held_samples": np.array([0, 1]), "held_weights": np.array([0.5, 0.5]), "truth_sample": 2 + i}})
        delta = (mats[0, [0, 1], 0][:, None] - mats[1, [0, 1], 0][None, :]).reshape(-1)
        joint = np.full(4, 0.25)
        truth = float(mats[0, 2 + i, 0] - mats[1, 2 + i, 0])
        pred = pl2.evaluate_truth_case(delta_b=delta, baseline_weights=joint, truth_delta_b=truth)
        lam_rows = {format(float(x["lambda_total"]), ".17g"): x for x in pred["lambda_results"]}
        row = {"truth_pair_id": f"P{i}", "reaction_id": "R", "truth_status": "NON_TIE", "baseline_magnitude_estimate": core.analyze_distribution(delta, joint, truth)["baseline"]["mean"], "lambda_0_25_random_sign_expected_gain": 0.0}
        for lam in (0.25, 0.5, 1.0):
            pr = lam_rows[format(lam, ".17g")]
            row[f"lambda_{str(lam).replace('.', '_')}_correct_absolute_error_gain"] = pr["correct_abs_error_gain"]
            row[f"lambda_{str(lam).replace('.', '_')}_wrong_absolute_error_gain"] = pr["wrong_abs_error_gain"]
        cache_rows.append(row)

    manifest = tmp_path / "parts.json"
    manifest.write_text(json.dumps({"parts": [{"filename": "p0"}, {"filename": "p1"}]}))
    monkeypatch.setattr(runner, "_case_part_reader", lambda _paths, idx: (io.StringIO(""), iter([cache_rows[idx]])))
    monkeypatch.setattr(pl2c, "case_responses", lambda **_kwargs: {"directionally_useful": False})
    cache_index = {("M", "CT2A", "e0"): 0, ("M", "GL261", "e1"): 1}

    def execute(path, *, interrupt=False):
        path.mkdir()
        (path / "RUN_STATE.json").write_text(json.dumps({"completed_parts": {}}))
        if interrupt:
            original = runner._case_part_reader
            fired = {"value": False}
            def fail_once(paths, idx):
                if idx == 1 and not fired["value"]:
                    fired["value"] = True
                    raise RuntimeError("synthetic interruption")
                return original(paths, idx)
            monkeypatch.setattr(runner, "_case_part_reader", fail_once)
            with pytest.raises(RuntimeError, match="synthetic interruption"):
                runner._analyze_anchor("A1", specs, ["R"], {"R": 0}, mats, cache_index,
                                       {"manifest_path": manifest, "artifacts": {"p0": manifest, "p1": manifest}}, path, ["E"])
            monkeypatch.setattr(runner, "_case_part_reader", original)
        return runner._analyze_anchor("A1", specs, ["R"], {"R": 0}, mats, cache_index,
                                      {"manifest_path": manifest, "artifacts": {"p0": manifest, "p1": manifest}}, path, ["E"])

    resumed = execute(tmp_path / "resumed", interrupt=True)[0]
    fresh = execute(tmp_path / "fresh")[0]
    assert np.array_equal(resumed.truth_counts, fresh.truth_counts)
    assert resumed.counts.keys() == fresh.counts.keys()
    for key in resumed.sums:
        assert np.allclose(resumed.sums[key], fresh.sums[key], equal_nan=True)
