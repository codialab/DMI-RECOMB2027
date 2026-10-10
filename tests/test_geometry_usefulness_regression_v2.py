import numpy as np
import pandas as pd

from analyses import geometry_usefulness_regression_v2 as v2


def test_top_indices_uses_ceil_and_stable_reaction_id_tiebreak():
    score = np.array([0.9, 0.9, 0.2, 0.1])
    ids = np.array(["R2", "R1", "R3", "R4"])
    selected = v2.top_indices(score, ids, 0.26)
    assert selected.tolist() == [1, 0]


def test_paired_comparison_reports_m4_minus_m3_topk_changes():
    ids = np.array([f"R{i}" for i in range(10)])
    y = np.array([1, 1, 1, 0, 0, 0, 0, 0, 0, 0], dtype=float)
    gain = np.array([0.2, 0.1, 0.3, -0.1, -0.2, -0.3, -0.1, -0.2, -0.1, -0.1])
    pred3 = np.arange(10, dtype=float)
    pred4 = y.copy()
    result = v2.comparison_metrics(ids, y, gain, pred3, pred4)
    assert result["mae_improvement_m3_minus_m4"] > 0
    assert result["spearman_change_m4_minus_m3"] > 0
    assert result["top10_usefulness_change_m4_minus_m3"] > 0
    assert result["top20_usefulness_change_m4_minus_m3"] > 0


def test_gain_topk_reports_positive_gain_fraction_for_both_ranking_targets():
    gain_pred = pd.DataFrame({
        "reaction_id": ["R1", "R2", "R3", "R4"],
        "subsystem": ["S1", "S1", "S2", "S2"],
        "observed_gain": [0.3, -0.2, 0.1, -0.1],
        "observed_usefulness": [1.0, 0.0, 1.0, 0.0],
        "gain_M1": [0.4, 0.3, 0.2, 0.1],
        "gain_M2": [0.4, 0.3, 0.2, 0.1],
        "gain_M3": [0.4, 0.3, 0.2, 0.1],
        "gain_M4": [0.4, 0.3, 0.2, 0.1],
    })
    freq = pd.DataFrame([
        {"anchor": "A1", "reaction_id": rid, "model": model, "predicted_usefulness_frequency": score}
        for model in ("M1", "M2", "M3", "M4")
        for rid, score in zip(["R1", "R2", "R3", "R4"], [0.1, 0.2, 0.3, 0.4])
    ])
    result = v2.make_gain_topk("A1", gain_pred, freq)
    top10 = result.loc[result.top_fraction.eq(0.1)]
    assert set(top10.ranking_target) == {"gain", "usefulness_frequency"}
    assert top10.positive_mean_gain_fraction.notna().all()
    assert np.allclose(top10.positive_mean_gain_fraction, [1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0])


def test_calibration_bins_and_topk_share_the_same_finite_population():
    frame = pd.DataFrame({
        "reaction_id": [f"R{i}" for i in range(10)],
        "model": ["M1"] * 10,
        "predicted_usefulness_frequency": np.linspace(0.05, 0.5, 10),
        "observed_usefulness_frequency": np.linspace(0.0, 0.9, 10),
    })
    rows = pd.DataFrame(v2.calibration_rows("A1", frame))
    bins = rows.loc[rows.metric_scope.eq("bin")]
    assert bins.n.sum() == len(frame)
    assert bins.n.tolist() == [2, 2, 2, 2, 2]
    assert rows.loc[rows.metric_scope.eq("overall"), "n"].iloc[0] == len(frame)
    assert set(rows.loc[rows.metric_scope.eq("top_k"), "n"]) == {1, 2}
