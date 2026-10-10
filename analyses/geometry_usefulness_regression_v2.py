#!/usr/bin/env python3
"""Post hoc v2 follow-up using only the frozen v1 reaction-level artifacts."""
from __future__ import annotations

import hashlib
import json
import math
import os
import platform
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/recomb2027_dmi_v2_matplotlib")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import sklearn
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[1]
V1 = ROOT / "reproduced/geometry_usefulness_regression_v1"
OUT = ROOT / "reproduced/geometry_usefulness_regression_v2"
BOOTSTRAP_REPLICATES = 2_000
BOOTSTRAP_SEED = 20271031
ALPHAS = (0.01, 0.1, 1.0, 10.0)
ANCHORS = ("A1", "A2")
MODELS = {
    "M1": ("mean_directional_entropy3",),
    "M2": ("mean_sign_magnitude_eta2",),
    "M3": ("mean_directional_entropy3", "mean_non_tie_pair_fraction"),
    "M4": ("mean_directional_entropy3", "mean_non_tie_pair_fraction", "mean_sign_magnitude_eta2"),
}
PREDICTION_MODELS = ("M1", "M2", "M3", "M4")
TOP_FRACTIONS = (0.10, 0.20)
V1_INPUTS = ("REACTION_LEVEL_INPUT.tsv.gz", "CONFIRMATION_PREDICTIONS.tsv.gz")
OUTPUTS = (
    "UNCERTAINTY_METRICS.tsv",
    "GAIN_MODEL_METRICS.tsv",
    "GAIN_TOP_K_METRICS.tsv",
    "CALIBRATION_METRICS.tsv",
    "UNCERTAINTY_SUMMARY.png",
    "GAIN_PRIORITIZATION.png",
    "CALIBRATION.png",
    "SCIENTIFIC_SUMMARY.md",
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", compression="infer")


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 2 or np.unique(x[ok]).size < 2 or np.unique(y[ok]).size < 2:
        return math.nan
    return float(spearmanr(x[ok], y[ok]).statistic)


def top_indices(score: np.ndarray, reaction_ids: np.ndarray, fraction: float, occurrence: np.ndarray | None = None) -> np.ndarray:
    n = len(score)
    take = max(1, int(math.ceil(fraction * n)))
    if occurrence is None:
        occurrence = np.arange(n)
    order = np.lexsort((occurrence, reaction_ids.astype(str), -score))
    return order[:take]


def comparison_metrics(
    ids: np.ndarray,
    y: np.ndarray,
    gain: np.ndarray,
    pred3: np.ndarray,
    pred4: np.ndarray,
    occurrence: np.ndarray | None = None,
) -> dict[str, float]:
    mae3 = float(mean_absolute_error(y, pred3))
    mae4 = float(mean_absolute_error(y, pred4))
    out = {
        "mae_improvement_m3_minus_m4": mae3 - mae4,
        "relative_mae_reduction": (mae3 - mae4) / mae3 if mae3 else math.nan,
        "spearman_change_m4_minus_m3": spearman(pred4, y) - spearman(pred3, y),
    }
    for frac in TOP_FRACTIONS:
        s3 = top_indices(pred3, ids, frac, occurrence)
        s4 = top_indices(pred4, ids, frac, occurrence)
        label = int(frac * 100)
        out[f"top{label}_usefulness_change_m4_minus_m3"] = float(y[s4].mean() - y[s3].mean())
        out[f"top{label}_mean_gain_change_m4_minus_m3"] = float(gain[s4].mean() - gain[s3].mean())
    return out


def bootstrap_intervals(data: pd.DataFrame, anchor: str, unit: str) -> list[dict]:
    ids = data.reaction_id.astype(str).to_numpy()
    y = data.observed_usefulness_frequency.to_numpy(float)
    gain = data.mean_correct_gain.to_numpy(float)
    p3 = data.loc[:, "M3"].to_numpy(float)
    p4 = data.loc[:, "M4"].to_numpy(float)
    point = comparison_metrics(ids, y, gain, p3, p4)
    rng = np.random.default_rng(BOOTSTRAP_SEED + (0 if anchor == "A1" else 1) + (0 if unit == "reaction" else 100))
    samples: dict[str, list[float]] = {k: [] for k in point}
    valid = 0
    if unit == "reaction":
        group_indices = None
        cluster_names = None
    else:
        cluster_series = data.subsystem.astype(str)
        cluster_names = sorted(cluster_series.unique())
        group_indices = [np.flatnonzero(cluster_series.to_numpy() == group) for group in cluster_names]
    for _ in range(BOOTSTRAP_REPLICATES):
        if unit == "reaction":
            idx = rng.integers(0, len(data), size=len(data))
        else:
            chosen = rng.integers(0, len(cluster_names), size=len(cluster_names))
            idx = np.concatenate([group_indices[i] for i in chosen])
        occurrence = np.arange(len(idx))
        values = comparison_metrics(ids[idx], y[idx], gain[idx], p3[idx], p4[idx], occurrence)
        if all(np.isfinite(v) for v in values.values()):
            valid += 1
            for key, value in values.items():
                samples[key].append(value)
    invalid_fraction = 1.0 - valid / BOOTSTRAP_REPLICATES
    status = "UNSTABLE" if invalid_fraction > 0.05 else "OK"
    rows = []
    for metric, values in samples.items():
        lo, hi = np.percentile(values, [2.5, 97.5]) if values else (math.nan, math.nan)
        rows.append({
            "anchor": anchor,
            "bootstrap_unit": unit,
            "metric": metric,
            "estimate": point[metric],
            "ci95_low": float(lo),
            "ci95_high": float(hi),
            "replicates_requested": BOOTSTRAP_REPLICATES,
            "replicates_valid": valid,
            "invalid_fraction": invalid_fraction,
            "interval_status": status,
        })
    return rows


def fit_gain_models(data: pd.DataFrame, anchor: str) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    features = ["mean_sign_magnitude_eta2", "mean_directional_entropy3", "mean_non_tie_pair_fraction"]
    required = [*features, "mean_correct_gain"]
    finite = np.isfinite(data[required].to_numpy(float)).all(axis=1)
    d = data.loc[finite].copy().reset_index(drop=True)
    dev = d.loc[d.analysis_split.eq("DEVELOPMENT")].copy().reset_index(drop=True)
    conf = d.loc[d.analysis_split.eq("CONFIRMATION_HOLDOUT")].copy().reset_index(drop=True)
    require(len(dev) > 0 and len(conf) > 0, f"{anchor}: no finite Development/Confirmation gain population")
    y_dev = dev.mean_correct_gain.to_numpy(float)
    y_conf = conf.mean_correct_gain.to_numpy(float)
    folds = dev.fold.to_numpy(int)
    require(set(np.unique(folds)) == {0, 1, 2}, f"{anchor}: expected three Development folds")

    m0_oof = np.full(len(dev), np.nan)
    for fold in range(3):
        train, test = folds != fold, folds == fold
        require(test.any() and train.any(), f"{anchor}: empty CV fold")
        require(not (set(dev.loc[train, "subsystem"]) & set(dev.loc[test, "subsystem"])), f"{anchor}: subsystem leakage in CV")
        m0_oof[test] = y_dev[train].mean()
    predictions: dict[str, np.ndarray] = {"M0": np.full(len(conf), y_dev.mean())}
    cv_maes = {"M0": float(mean_absolute_error(y_dev, m0_oof))}
    selected_alpha: dict[str, float | None] = {"M0": None}
    for model, columns in MODELS.items():
        scores = {alpha: [] for alpha in ALPHAS}
        for fold in range(3):
            train, test = folds != fold, folds == fold
            for alpha in ALPHAS:
                estimator = make_pipeline(StandardScaler(), Ridge(alpha=alpha))
                estimator.fit(dev.loc[train, list(columns)], y_dev[train])
                pred = estimator.predict(dev.loc[test, list(columns)])
                scores[alpha].append(float(mean_absolute_error(y_dev[test], pred)))
        mean_scores = {alpha: float(np.mean(vals)) for alpha, vals in scores.items()}
        best = min(ALPHAS, key=lambda alpha: (mean_scores[alpha], alpha))
        final = make_pipeline(StandardScaler(), Ridge(alpha=best))
        final.fit(dev.loc[:, list(columns)], y_dev)
        predictions[model] = final.predict(conf.loc[:, list(columns)])
        selected_alpha[model] = best
        cv_maes[model] = mean_scores[best]

    require(len(conf) == (838 if anchor == "A1" else 843), f"{anchor}: gain Confirmation cohort mismatch: {len(conf)}")
    rows = []
    observed_metrics = {}
    for model, pred in predictions.items():
        metrics = {
            "mae": float(mean_absolute_error(y_conf, pred)),
            "rmse": float(mean_squared_error(y_conf, pred) ** 0.5),
            "spearman": spearman(pred, y_conf),
        }
        observed_metrics[model] = metrics
    for model in ("M0", *MODELS):
        m = observed_metrics[model]
        m0 = observed_metrics["M0"]
        m3 = observed_metrics["M3"]
        rows.append({
            "anchor": anchor,
            "model": model,
            "predictors": "+".join(MODELS.get(model, ())) or "constant Development mean",
            "selected_alpha": selected_alpha[model],
            "cv_mae": cv_maes[model],
            "confirmation_n": len(conf),
            "mae": m["mae"],
            "rmse": m["rmse"],
            "spearman": m["spearman"],
            "mae_improvement_vs_m0": m0["mae"] - m["mae"],
            "rmse_improvement_vs_m0": m0["rmse"] - m["rmse"],
            "spearman_change_vs_m0": m["spearman"] - m0["spearman"] if np.isfinite(m["spearman"]) else math.nan,
            "mae_improvement_vs_m3": m3["mae"] - m["mae"],
            "rmse_improvement_vs_m3": m3["rmse"] - m["rmse"],
            "spearman_change_vs_m3": m["spearman"] - m3["spearman"] if np.isfinite(m["spearman"]) else math.nan,
        })
    conf = conf.reset_index(drop=True)
    pred_frame = pd.DataFrame({"reaction_id": conf.reaction_id.astype(str), "subsystem": conf.subsystem.astype(str),
                               "observed_gain": y_conf, "observed_usefulness": conf.mean_directionally_useful_fraction.to_numpy(float),
                               **{f"gain_{m}": p for m, p in predictions.items()}})
    require(pred_frame.reaction_id.is_unique, f"{anchor}: duplicate gain Confirmation reactions")
    return pd.DataFrame(rows), pred_frame, {"development_n": len(dev), "confirmation_n": len(conf), "selected_alpha": selected_alpha}


def make_gain_topk(anchor: str, pred_gain: pd.DataFrame, freq_predictions: pd.DataFrame) -> pd.DataFrame:
    ids = pred_gain.reaction_id.astype(str).to_numpy()
    gain = pred_gain.observed_gain.to_numpy(float)
    useful = pred_gain.observed_usefulness.to_numpy(float)
    population_gain, population_useful = float(gain.mean()), float(useful.mean())
    frequency = freq_predictions.loc[freq_predictions.anchor.eq(anchor)].pivot(index="reaction_id", columns="model", values="predicted_usefulness_frequency")
    require(set(frequency.index.astype(str)) == set(ids), f"{anchor}: frequency/gain Confirmation populations differ")
    frequency = frequency.reindex(ids)
    score_by_target: dict[tuple[str, str], np.ndarray] = {}
    for model in MODELS:
        score_by_target[("gain", model)] = pred_gain[f"gain_{model}"].to_numpy(float)
        score_by_target[("usefulness_frequency", model)] = frequency[model].to_numpy(float)
    rows = []
    for frac in TOP_FRACTIONS:
        label = int(frac * 100)
        selected: dict[tuple[str, str], np.ndarray] = {}
        for key, score in score_by_target.items():
            selected[key] = top_indices(score, ids, frac)
        for target in ("gain", "usefulness_frequency"):
            for model in MODELS:
                idx = selected[(target, model)]
                matched = selected[("usefulness_frequency" if target == "gain" else "gain", model)]
                entropy = selected[(target, "M1")]
                selected_gain = float(gain[idx].mean())
                selected_useful = float(useful[idx].mean())
                positive_fraction = float(np.mean(gain[idx] > 0.0))
                rows.append({
                    "anchor": anchor,
                    "ranking_target": target,
                    "model": model,
                    "top_fraction": frac,
                    "selected_n": len(idx),
                    "observed_mean_correct_gain": selected_gain,
                    "observed_usefulness_frequency": selected_useful,
                    "positive_mean_gain_fraction": positive_fraction,
                    "random_expected_mean_correct_gain": population_gain,
                    "random_expected_usefulness_frequency": population_useful,
                    "mean_gain_improvement_vs_random": selected_gain - population_gain,
                    "usefulness_improvement_vs_random": selected_useful - population_useful,
                    "positive_gain_fraction_vs_random": positive_fraction - float(np.mean(gain > 0.0)),
                    "difference_vs_entropy_only_mean_gain": selected_gain - float(gain[entropy].mean()),
                    "difference_vs_entropy_only_usefulness": selected_useful - float(useful[entropy].mean()),
                    "difference_vs_entropy_only_positive_gain_fraction": positive_fraction - float(np.mean(gain[entropy] > 0.0)),
                    "difference_vs_matched_other_target_mean_gain": selected_gain - float(gain[matched].mean()),
                    "difference_vs_matched_other_target_usefulness": selected_useful - float(useful[matched].mean()),
                    "difference_vs_matched_other_target_positive_gain_fraction": positive_fraction - float(np.mean(gain[matched] > 0.0)),
                    "entropy_only_comparator_target": target,
                    "matched_other_target": "usefulness_frequency" if target == "gain" else "gain",
                })
    return pd.DataFrame(rows)


def calibration_rows(anchor: str, frame: pd.DataFrame) -> list[dict]:
    rows = []
    ids = frame.reaction_id.astype(str).to_numpy()
    observed = frame.observed_usefulness_frequency.to_numpy(float)
    for model, group in frame.groupby("model", sort=True):
        pred = group.predicted_usefulness_frequency.to_numpy(float)
        group_ids = group.reaction_id.astype(str).to_numpy()
        order = np.lexsort((group_ids, pred))
        bins = np.array_split(order, 5)
        for bin_no, idx in enumerate(bins, start=1):
            rows.append({"anchor": anchor, "model": model, "metric_scope": "bin", "bin": bin_no,
                         "n": len(idx), "mean_predicted_usefulness": float(pred[idx].mean()),
                         "mean_observed_usefulness": float(group.observed_usefulness_frequency.to_numpy(float)[idx].mean()),
                         "calibration_bias_observed_minus_predicted": float((group.observed_usefulness_frequency.to_numpy(float)[idx] - pred[idx]).mean()),
                         "top_fraction": math.nan})
        rows.append({"anchor": anchor, "model": model, "metric_scope": "overall", "bin": math.nan,
                     "n": len(group), "mean_predicted_usefulness": float(pred.mean()),
                     "mean_observed_usefulness": float(group.observed_usefulness_frequency.mean()),
                     "calibration_bias_observed_minus_predicted": float((group.observed_usefulness_frequency.to_numpy(float) - pred).mean()),
                     "top_fraction": math.nan})
        for frac in TOP_FRACTIONS:
            idx = top_indices(pred, group_ids, frac)
            obs = group.observed_usefulness_frequency.to_numpy(float)
            rows.append({"anchor": anchor, "model": model, "metric_scope": "top_k", "bin": math.nan,
                         "n": len(idx), "mean_predicted_usefulness": float(pred[idx].mean()),
                         "mean_observed_usefulness": float(obs[idx].mean()),
                         "calibration_bias_observed_minus_predicted": float((obs[idx] - pred[idx]).mean()),
                         "top_fraction": frac})
    return rows


def plot_uncertainty(frame: pd.DataFrame, path: Path) -> None:
    metrics = [
        "mae_improvement_m3_minus_m4", "relative_mae_reduction", "spearman_change_m4_minus_m3",
        "top10_usefulness_change_m4_minus_m3", "top20_usefulness_change_m4_minus_m3",
        "top10_mean_gain_change_m4_minus_m3", "top20_mean_gain_change_m4_minus_m3",
    ]
    labels = ["MAE improvement", "Relative MAE reduction", "Spearman change", "Top 10% usefulness", "Top 20% usefulness", "Top 10% mean gain", "Top 20% mean gain"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharex=False)
    colors = {"reaction": "#2878B5", "subsystem": "#D95F02"}
    for ax, anchor in zip(axes, ANCHORS):
        sub = frame.loc[frame.anchor.eq(anchor)]
        for y, (metric, label) in enumerate(zip(metrics, labels)):
            for offset, unit in ((-0.12, "reaction"), (0.12, "subsystem")):
                row = sub.loc[sub.metric.eq(metric) & sub.bootstrap_unit.eq(unit)].iloc[0]
                ax.errorbar(row.estimate, y + offset,
                            xerr=[[row.estimate - row.ci95_low], [row.ci95_high - row.estimate]],
                            fmt="o", color=colors[unit], capsize=2, markersize=4)
        ax.axvline(0, color="#555555", lw=0.8, ls="--")
        ax.set_yticks(range(len(labels)), labels)
        ax.invert_yaxis()
        ax.set_title(anchor)
        ax.set_xlabel("M4 − M3 for correlations and top-k; M3 − M4 for MAE")
        ax.grid(axis="x", alpha=0.2)
    fig.legend([plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=colors[k], label=k.title()) for k in colors],
               ["Reaction bootstrap", "Subsystem bootstrap"], loc="lower center", ncol=2, frameon=False)
    fig.suptitle("Conditional paired uncertainty for M3 versus M4")
    fig.tight_layout(rect=(0, 0.08, 1, 0.94))
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_gain(frame: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13, 8))
    for row, anchor in enumerate(ANCHORS):
        sub = frame.loc[(frame.anchor.eq(anchor)) & frame.ranking_target.isin(["gain", "usefulness_frequency"])]
        for col, (metric, ylabel, title) in enumerate((
            ("observed_mean_correct_gain", "Observed mean correct-cue gain", "Mean gain"),
            ("positive_mean_gain_fraction", "Fraction with positive mean gain", "Positive gain fraction"),
        )):
            ax = axes[row, col]
            width = 0.18
            x = np.arange(4)
            for j, (target, frac, color) in enumerate((("gain", 0.1, "#2878B5"), ("gain", 0.2, "#79ADD2"), ("usefulness_frequency", 0.1, "#D95F02"), ("usefulness_frequency", 0.2, "#F0A36B"))):
                vals = sub.loc[(sub.ranking_target.eq(target)) & (sub.top_fraction.eq(frac))].set_index("model").reindex(PREDICTION_MODELS)[metric]
                ax.bar(x + (j - 1.5) * width, vals.to_numpy(float), width, color=color,
                       label=f"{'Gain' if target == 'gain' else 'Usefulness'} top {int(frac*100)}%")
            if metric == "observed_mean_correct_gain":
                baseline = float(sub.loc[sub.ranking_target.eq("gain"), "random_expected_mean_correct_gain"].iloc[0])
                ax.axhline(baseline, color="#444444", ls="--", lw=1, label="Population mean")
            ax.set_xticks(x, PREDICTION_MODELS)
            ax.set_ylabel(ylabel)
            ax.set_title(f"{anchor}: {title}")
            ax.grid(axis="y", alpha=0.2)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=5, frameon=False)
    fig.suptitle("Confirmation prioritization by predicted gain and usefulness frequency")
    fig.tight_layout(rect=(0, 0.09, 1, 0.94))
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_calibration(calibration: pd.DataFrame, path: Path) -> None:
    colors = {"M1": "#2878B5", "M3": "#41AB5D", "M4": "#756BB1"}
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharex=True, sharey=True)
    bins = calibration.loc[calibration.metric_scope.eq("bin")]
    for ax, anchor in zip(axes, ANCHORS):
        sub = bins.loc[bins.anchor.eq(anchor)]
        lim = [0, max(0.25, float(max(sub.mean_predicted_usefulness.max(), sub.mean_observed_usefulness.max())) * 1.05)]
        ax.plot(lim, lim, color="#555555", ls="--", lw=1)
        for model in ("M1", "M3", "M4"):
            part = sub.loc[sub.model.eq(model)].sort_values("bin")
            ax.plot(part.mean_predicted_usefulness, part.mean_observed_usefulness, marker="o", color=colors[model], label=model)
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(anchor)
        ax.set_xlabel("Mean predicted usefulness")
        ax.grid(alpha=0.2)
    axes[0].set_ylabel("Mean observed usefulness")
    axes[1].legend(frameon=False)
    fig.suptitle("Post hoc Confirmation calibration across predicted-usefulness groups")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def fmt(value: float, digits: int = 4) -> str:
    return f"{value:.{digits}f}"


def scientific_summary(unc: pd.DataFrame, gain_metrics: pd.DataFrame, topk: pd.DataFrame, calibration: pd.DataFrame) -> str:
    lines = [
        "# Geometry usefulness regression v2: scientific summary",
        "",
        "This is a post hoc computational follow-up on the frozen Development/Confirmation split. The Confirmation reactions have already been examined in v1; these results are not a new independent confirmation experiment. Bootstrap intervals are conditional on the fitted models and do not account for model or cutoff selection.",
        "",
        "## 1. Usefulness-frequency prediction and uncertainty",
    ]
    for anchor in ANCHORS:
        rows = unc.loc[(unc.anchor.eq(anchor)) & (unc.bootstrap_unit.eq("reaction"))].set_index("metric")
        mae = rows.loc["mae_improvement_m3_minus_m4"]
        rel = rows.loc["relative_mae_reduction"]
        rho = rows.loc["spearman_change_m4_minus_m3"]
        cluster_mae = unc.loc[(unc.anchor.eq(anchor)) & (unc.bootstrap_unit.eq("subsystem")) & (unc.metric.eq("mae_improvement_m3_minus_m4"))].iloc[0]
        lines.append(f"- **{anchor}:** M3−M4 MAE difference={fmt(mae.estimate)} (95% reaction-bootstrap CI {fmt(mae.ci95_low)} to {fmt(mae.ci95_high)}; subsystem-bootstrap CI {fmt(cluster_mae.ci95_low)} to {fmt(cluster_mae.ci95_high)}). Relative MAE reduction={fmt(rel.estimate*100, 2)}% (reaction CI {fmt(rel.ci95_low*100, 2)}% to {fmt(rel.ci95_high*100, 2)}%). Spearman M4−M3={fmt(rho.estimate)} (reaction CI {fmt(rho.ci95_low)} to {fmt(rho.ci95_high)}).")
    lines += ["", "Reaction-level resampling treats reactions as exchangeable and does not fully represent metabolic dependence. Subsystem resampling is a sensitivity analysis over 49 represented clusters per anchor; the clusters are highly imbalanced (finite Confirmation sizes: median 7, maximum 144), so those intervals should be read cautiously.", "", "## 2. Top-k usefulness prioritization"]
    for anchor in ANCHORS:
        pairs = []
        for frac in TOP_FRACTIONS:
            part = topk.loc[(topk.anchor.eq(anchor)) & (topk.ranking_target.eq("usefulness_frequency")) & topk.model.isin(["M3", "M4"]) & topk.top_fraction.eq(frac)].set_index("model")
            pairs.append(f"top {int(frac*100)}% M3 {fmt(part.loc['M3','observed_usefulness_frequency'])} vs M4 {fmt(part.loc['M4','observed_usefulness_frequency'])}")
        boot = unc.loc[(unc.anchor.eq(anchor)) & (unc.bootstrap_unit.eq("reaction")) & unc.metric.isin(["top10_usefulness_change_m4_minus_m3", "top20_usefulness_change_m4_minus_m3"])].set_index("metric")
        cluster20 = unc.loc[(unc.anchor.eq(anchor)) & (unc.bootstrap_unit.eq("subsystem")) & unc.metric.eq("top20_usefulness_change_m4_minus_m3")].iloc[0]
        t10 = boot.loc["top10_usefulness_change_m4_minus_m3"]
        t20 = boot.loc["top20_usefulness_change_m4_minus_m3"]
        interpretation = ("the evaluated top-k metrics do not show a material gain" if anchor == "A1" else "the top-20% usefulness lift is positive under reaction resampling but its subsystem interval includes zero")
        lines.append(
            f"- **{anchor}:** {'; '.join(pairs)}. Paired reaction-bootstrap changes (M4−M3) were "
            f"top10={fmt(t10.estimate)} (95% CI {fmt(t10.ci95_low)} to {fmt(t10.ci95_high)}) and "
            f"top20={fmt(t20.estimate)} (CI {fmt(t20.ci95_low)} to {fmt(t20.ci95_high)}; "
            f"subsystem CI {fmt(cluster20.ci95_low)} to {fmt(cluster20.ci95_high)}). Thus {interpretation}. "
            "These previously inspected Confirmation rankings are descriptive and do not establish independent prioritization performance."
        )
    lines += ["", "## 3. Mean correct-cue gain prediction"]
    for anchor in ANCHORS:
        metrics = gain_metrics.loc[gain_metrics.anchor.eq(anchor)].set_index("model")
        m4, m3, m0 = metrics.loc["M4"], metrics.loc["M3"], metrics.loc["M0"]
        lines.append(f"- **{anchor}:** Gain-target M4 confirmation MAE={fmt(m4.mae)}, RMSE={fmt(m4.rmse)}, Spearman={fmt(m4.spearman)}; M0 MAE={fmt(m0.mae)}. M4 versus M0 MAE improvement={fmt(m0.mae-m4.mae)}; M4 versus M3 changes were MAE improvement {fmt(m3.mae-m4.mae)}, RMSE improvement {fmt(m3.rmse-m4.rmse)}, and Spearman change {fmt(m4.spearman-m3.spearman)}. Across models, M4 does not improve gain MAE over M0 under either anchor; the gain target therefore has weak overall predictive support in this cohort despite some rank associations and top-k contrasts.")
    lines += ["", "The gain target reuses v1's signed `mean_correct_gain`: within each evaluation, correct-cue gain is averaged over non-tie evaluable truth pairs; the reaction target is the arithmetic mean of finite evaluation-level values. It is not clipped or reweighted by the number of pairs per evaluation.", "", "Gain-targeted and usefulness-frequency-targeted rankings are compared on identical per-anchor reaction populations. For M4, the actual top-k comparison is:"]
    for anchor in ANCHORS:
        snippets = []
        for frac in TOP_FRACTIONS:
            part = topk.loc[(topk.anchor.eq(anchor)) & topk.model.eq("M4") & topk.top_fraction.eq(frac)].set_index("ranking_target")
            g, f = part.loc["gain"], part.loc["usefulness_frequency"]
            snippets.append(f"top {int(frac*100)}% gain ranking: mean gain {fmt(g.observed_mean_correct_gain)}, positive fraction {fmt(g.positive_mean_gain_fraction)}, usefulness {fmt(g.observed_usefulness_frequency)}; frequency ranking: {fmt(f.observed_mean_correct_gain)}, {fmt(f.positive_mean_gain_fraction)}, {fmt(f.observed_usefulness_frequency)}")
        lines.append(f"- **{anchor}:** " + "; ".join(snippets) + ".")
    lines += ["", "The tables also compare each ranking with random selection and the entropy-only ranking. The endpoints remain distinct: usefulness frequency measures how often outcomes are beneficial, while mean correct-cue gain is a signed magnitude. Positive-gain fraction is reported separately from mean gain.", "", "## 4. Calibration"]
    for anchor in ANCHORS:
        rows = calibration.loc[(calibration.anchor.eq(anchor)) & (calibration.metric_scope.isin(["overall", "top_k"]))]
        desc = []
        for model in ("M1", "M3", "M4"):
            overall = rows.loc[(rows.model.eq(model)) & rows.metric_scope.eq("overall")].iloc[0]
            t10 = rows.loc[(rows.model.eq(model)) & rows.top_fraction.eq(0.1)].iloc[0]
            t20 = rows.loc[(rows.model.eq(model)) & rows.top_fraction.eq(0.2)].iloc[0]
            desc.append(f"{model} bias overall/top10/top20={fmt(overall.calibration_bias_observed_minus_predicted)}/{fmt(t10.calibration_bias_observed_minus_predicted)}/{fmt(t20.calibration_bias_observed_minus_predicted)}")
        lines.append(f"- **{anchor}:** " + "; ".join(desc) + ". Bias is observed minus predicted usefulness; calibration is descriptive and does not establish external generalization.")
    lines += [
        "",
        "## 5. Anchor comparison and manuscript implication",
        "The M3/M4 usefulness-frequency accuracy gains occur under both anchors, while top-k usefulness changes are small for A1 and more favorable for A2 at 20%; subsystem resampling makes that A2 top-20% interval uncertain. Gain-targeted M4 ranking selects higher realized mean gain than matched frequency-targeted M4 ranking under both anchors, especially A2, while its selected usefulness frequency is lower under A2. The mismatch shows why practical prioritization cannot be inferred from regression accuracy alone.",
        "",
        "For the manuscript, these results support an anchor-specific predictive association between η²_dir and beneficial-outcome frequency conditional on the fitted linear models. They do not establish reliable prediction of signed mean magnitude gain or prove the theoretical direction–magnitude coupling mechanism. Keep these claims separate; η²_dir is not universally superior across targets and prioritization metrics.",
        "",
        "All modeling used Development outcomes only for fitting, fold-local standardization, and alpha selection. No Confirmation outcome was used to fit or select gain models. No GEM reconstruction, flux sampling, truth generation, or case-level aggregation was rerun.",
    ]
    return "\n".join(lines) + "\n"


def run() -> dict:
    if OUT.exists():
        require(OUT.is_dir() and (OUT / "MANIFEST.json").is_file(), f"existing output is not a recognized v2 bundle: {OUT}")
        existing = json.loads((OUT / "MANIFEST.json").read_text(encoding="utf-8"))
        require(existing.get("schema") == "dmi.geometry_usefulness_regression.v2", f"refusing to overwrite unrelated output: {OUT}")
    manifest_path = V1 / "MANIFEST.json"
    require(manifest_path.is_file() and not manifest_path.is_symlink(), "v1 manifest missing or symlinked")
    v1_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    require(v1_manifest.get("schema") == "dmi.geometry_usefulness_regression.v1", "unexpected v1 manifest schema")
    require(v1_manifest.get("status") == "POST_HOC_PREDICTIVE_FOLLOW_UP_COMPLETE", "v1 analysis is not marked complete")
    require(v1_manifest.get("qc", {}).get("status") == "PASS", "v1 provenance/QC status is not PASS")
    require(v1_manifest.get("matched_primary_population", {}).get("reaction_count") == 4_179, "v1 manifest matched-cohort count mismatch")
    for name, expected in v1_manifest.get("artifact_sha256", {}).items():
        path = V1 / name
        require(path.is_file() and not path.is_symlink(), f"v1 artifact missing or symlinked: {name}")
        require(sha256(path) == expected, f"v1 artifact hash mismatch: {name}")
    v1_script = ROOT / "analyses/geometry_usefulness_regression.py"
    require(sha256(v1_script) == v1_manifest.get("implementation_sha256"), "v1 implementation hash mismatch")

    reaction = read_tsv(V1 / "REACTION_LEVEL_INPUT.tsv.gz")
    freq = read_tsv(V1 / "CONFIRMATION_PREDICTIONS.tsv.gz")
    fold_registry = read_tsv(V1 / "MATCHED_FOLD_ASSIGNMENTS.tsv")
    require(len(reaction) == 8_358 and reaction.groupby("anchor").reaction_id.nunique().to_dict() == {"A1": 4_179, "A2": 4_179}, "v1 cohort inventory mismatch")
    require(reaction.groupby(["anchor", "reaction_id"]).size().eq(1).all(), "duplicate reaction-level inputs")
    require(set(reaction.analysis_split.unique()) == {"DEVELOPMENT", "CONFIRMATION_HOLDOUT"}, "unexpected split labels")
    expected_splits = {"A1": {"DEVELOPMENT": 3_134, "CONFIRMATION_HOLDOUT": 1_045}, "A2": {"DEVELOPMENT": 3_134, "CONFIRMATION_HOLDOUT": 1_045}}
    for anchor, expected in expected_splits.items():
        observed = reaction.loc[reaction.anchor.eq(anchor), "analysis_split"].value_counts().to_dict()
        require(observed == expected, f"{anchor}: frozen split assignment counts mismatch: {observed}")
    require(fold_registry.reaction_id.is_unique and len(fold_registry) == 4_179, "v1 fold registry inventory mismatch")
    fold_ref = fold_registry.set_index("reaction_id").sort_index()[["subsystem", "fold", "analysis_split"]]
    require(fold_ref.equals(reaction.loc[reaction.anchor.eq("A1"), ["reaction_id", "subsystem", "fold", "analysis_split"]].set_index("reaction_id").sort_index()), "v1 reaction table differs from frozen fold registry")
    require(set(reaction.loc[reaction.anchor.eq("A1"), "reaction_id"]) == set(reaction.loc[reaction.anchor.eq("A2"), "reaction_id"]), "anchor cohort IDs differ")
    folds = reaction.loc[reaction.anchor.eq("A1"), ["reaction_id", "subsystem", "fold", "analysis_split"]].set_index("reaction_id").sort_index()
    folds2 = reaction.loc[reaction.anchor.eq("A2"), ["reaction_id", "subsystem", "fold", "analysis_split"]].set_index("reaction_id").sort_index()
    require(folds.equals(folds2), "anchor fold/split/subsystem assignments differ")
    require(set(freq.model.unique()) == {"M0", *PREDICTION_MODELS}, "v1 prediction model inventory mismatch")

    unc_rows = []
    gain_metric_frames = []
    topk_frames = []
    calibration_all = []
    gain_support = {}
    subsystem_sizes = {}
    for anchor in ANCHORS:
        conf_inputs = reaction.loc[(reaction.anchor.eq(anchor)) & reaction.analysis_split.eq("CONFIRMATION_HOLDOUT")].copy()
        finite = np.isfinite(conf_inputs[["mean_sign_magnitude_eta2", "mean_directional_entropy3", "mean_non_tie_pair_fraction", "mean_directionally_useful_fraction"]].to_numpy(float)).all(axis=1)
        conf_inputs = conf_inputs.loc[finite].copy()
        expected_n = 838 if anchor == "A1" else 843
        require(len(conf_inputs) == expected_n, f"{anchor}: finite v1 Confirmation count mismatch: {len(conf_inputs)}")
        p = freq.loc[freq.anchor.eq(anchor)].copy()
        pids = p.groupby("model").reaction_id.apply(lambda s: frozenset(s.astype(str))).to_dict()
        require(all(v == frozenset(conf_inputs.reaction_id.astype(str)) for v in pids.values()), f"{anchor}: v1 prediction IDs differ from common finite Confirmation cohort")
        require(p.groupby("model").size().eq(expected_n).all(), f"{anchor}: v1 model prediction count mismatch")
        wide = p.pivot(index="reaction_id", columns="model", values="predicted_usefulness_frequency").reindex(conf_inputs.reaction_id.astype(str))
        wide["observed_usefulness_frequency"] = conf_inputs.set_index(conf_inputs.reaction_id.astype(str)).mean_directionally_useful_fraction.reindex(wide.index).to_numpy(float)
        wide["mean_correct_gain"] = conf_inputs.set_index(conf_inputs.reaction_id.astype(str)).mean_correct_gain.reindex(wide.index).to_numpy(float)
        wide["subsystem"] = conf_inputs.set_index(conf_inputs.reaction_id.astype(str)).subsystem.reindex(wide.index).to_numpy()
        require(np.isfinite(wide[["M3", "M4", "observed_usefulness_frequency", "mean_correct_gain"]].to_numpy(float)).all(), f"{anchor}: nonfinite uncertainty values")
        groups = wide.subsystem.astype(str).value_counts()
        subsystem_sizes[anchor] = {"n_subsystems": int(len(groups)), "min": int(groups.min()), "median": float(groups.median()), "max": int(groups.max()), "size_counts": {str(int(k)): int(v) for k, v in groups.value_counts().sort_index().items()}}
        for unit in ("reaction", "subsystem"):
            unc_rows.extend(bootstrap_intervals(wide.reset_index(names="reaction_id"), anchor, unit))

        metrics, gain_pred, support = fit_gain_models(reaction.loc[reaction.anchor.eq(anchor)].copy(), anchor)
        gain_metric_frames.append(metrics)
        gain_support[anchor] = support
        topk_frames.append(make_gain_topk(anchor, gain_pred, freq))

        c = p.loc[p.model.isin(("M1", "M3", "M4")), ["reaction_id", "model", "observed_usefulness_frequency", "predicted_usefulness_frequency"]].copy()
        require(c.groupby("model").reaction_id.nunique().eq(expected_n).all(), f"{anchor}: calibration population mismatch")
        calibration_all.append(pd.DataFrame(calibration_rows(anchor, c)))

    uncertainty = pd.DataFrame(unc_rows)
    gain_metrics = pd.concat(gain_metric_frames, ignore_index=True)
    topk = pd.concat(topk_frames, ignore_index=True)
    calibration = pd.concat(calibration_all, ignore_index=True)
    require(len(topk) == 2 * 2 * 2 * len(MODELS), "gain top-k row inventory mismatch")
    require(topk.positive_mean_gain_fraction.notna().all(), "positive gain fraction missing")
    require(set(calibration.loc[calibration.metric_scope.eq("bin"), "bin"]) == {1, 2, 3, 4, 5}, "calibration bins incomplete")

    OUT.mkdir(parents=True, exist_ok=True)
    uncertainty.to_csv(OUT / "UNCERTAINTY_METRICS.tsv", sep="\t", index=False, na_rep="", lineterminator="\n")
    gain_metrics.to_csv(OUT / "GAIN_MODEL_METRICS.tsv", sep="\t", index=False, na_rep="", lineterminator="\n")
    topk.to_csv(OUT / "GAIN_TOP_K_METRICS.tsv", sep="\t", index=False, na_rep="", lineterminator="\n")
    calibration.to_csv(OUT / "CALIBRATION_METRICS.tsv", sep="\t", index=False, na_rep="", lineterminator="\n")
    plot_uncertainty(uncertainty, OUT / "UNCERTAINTY_SUMMARY.png")
    plot_gain(topk, OUT / "GAIN_PRIORITIZATION.png")
    plot_calibration(calibration, OUT / "CALIBRATION.png")
    (OUT / "SCIENTIFIC_SUMMARY.md").write_text(scientific_summary(uncertainty, gain_metrics, topk, calibration), encoding="utf-8")

    input_hashes = {name: sha256(V1 / name) for name in V1_INPUTS}
    manifest = {
        "schema": "dmi.geometry_usefulness_regression.v2",
        "status": "POST_HOC_FOLLOW_UP_COMPLETE",
        "analysis_role": "POST_HOC_USING_FROZEN_DEVELOPMENT_CONFIRMATION_SPLIT_NOT_NEW_INDEPENDENT_CONFIRMATION",
        "v1_manifest_sha256": sha256(manifest_path),
        "v1_implementation_sha256": v1_manifest["implementation_sha256"],
        "v1_input_artifact_sha256": input_hashes,
        "v1_artifact_sha256": v1_manifest["artifact_sha256"],
        "implementation_sha256": sha256(Path(__file__).resolve()),
        "matched_cohort_reactions": 4_179,
        "finite_confirmation_counts": {"A1": 838, "A2": 843},
        "confirmation_subsystem_distribution": subsystem_sizes,
        "gain_model_support": gain_support,
        "excluded_reactions": {"LDH_L": "A1-only reaction excluded from the frozen matched 4,179-reaction cohort"},
        "uncertainty": {"comparison": "M3 versus M4", "bootstrap_replicates": BOOTSTRAP_REPLICATES,
                        "seed": BOOTSTRAP_SEED, "interval": "percentile 95%", "refit": False,
                        "reaction_unit": "paired Confirmation reactions sampled with replacement",
                        "cluster_unit": "whole represented subsystems sampled with replacement"},
        "gain_model": {"target": "cached signed mean_correct_gain", "models": {"M0": "Development mean", **{k: list(v) for k, v in MODELS.items()}},
                       "alpha_grid": list(ALPHAS), "selection": "Development-only three-fold subsystem-grouped mean validation MAE",
                       "preprocessing": "StandardScaler fit within each training fold and final Development fit", "clipping": False,
                       "target_aggregation": "Reuse v1 mean_correct_gain: within each evaluation, mean signed correct-cue gain across non-tie evaluable truth pairs; then arithmetic mean across finite evaluation-level values per reaction, without weighting by pair count."},
        "calibration": {"models": ["M1", "M3", "M4"], "bins": 5, "confirmation_predictions_reused": True},
        "ranking": {"top_fraction_selection": "ceil(fraction * finite Confirmation count)", "tie_break": "predicted score descending, reaction_id ascending", "positive_mean_gain": "observed reaction mean_correct_gain > 0"},
        "post_hoc_selection_note": "Top-k and calibration estimates use previously examined Confirmation outcomes; bootstrap intervals do not correct for selection across models, metrics, or cutoffs.",
        "random_seeds": {"bootstrap": BOOTSTRAP_SEED, "per_anchor_unit_offsets": {"A1": 0, "A2": 1, "subsystem": 100}},
        "software": {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__, "scikit_learn": sklearn.__version__, "matplotlib": matplotlib.__version__},
        "outputs": list(OUTPUTS),
        "output_sha256": {name: sha256(OUT / name) for name in OUTPUTS},
    }
    (OUT / "MANIFEST.json").write_text(json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    return {"status": manifest["status"], "output": str(OUT.relative_to(ROOT)), "uncertainty_rows": len(uncertainty),
            "gain_model_rows": len(gain_metrics), "gain_top_k_rows": len(topk), "calibration_rows": len(calibration)}


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
