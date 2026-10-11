"""Numerical helpers for the additive weak-cue uncertainty analysis.

This module is deliberately independent of repository I/O and frozen PL2
implementations.  It operates on a finite weighted distribution of signed
contrasts and leaves frozen prediction-gain definitions untouched.
"""
from __future__ import annotations

from collections.abc import Mapping, Sequence
import math

import numpy as np

TIE_TOL = 1e-12
PRIMARY_LAMBDA = 0.25
LAMBDA_GRID = (0.0, 0.25, 0.5, 1.0)
RELATIVE_VARIANCE_EPS_MULTIPLIER = 100.0
VARIANCE_CHANGE_EPS_MULTIPLIER = 100.0
PREDICTORS = {
    "eta2": ("sign_magnitude_eta2",),
    "eta2_plus_entropy": ("sign_magnitude_eta2", "directional_entropy3"),
    "eta2_plus_dominant_mass": ("sign_magnitude_eta2", "dominant_sign_mass"),
    "eta2_plus_coverage": ("sign_magnitude_eta2", "non_tie_coverage"),
    "eta2_plus_sign_states": ("sign_magnitude_eta2", "n_supported_sign_states"),
    "all_geometry": (
        "sign_magnitude_eta2", "directional_entropy3", "dominant_sign_mass",
        "non_tie_coverage", "n_supported_sign_states",
    ),
}


def _vector(values: Sequence[float], name: str) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if array.ndim != 1 or array.size == 0 or not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be a non-empty finite 1-D vector")
    return array


def normalize_weights(weights: Sequence[float]) -> np.ndarray:
    weight = _vector(weights, "weights")
    if np.any(weight < 0.0):
        raise ValueError("weights must be non-negative")
    total = math.fsum(float(value) for value in weight)
    if not math.isfinite(total) or total <= 0.0:
        raise ValueError("weights must have positive finite mass")
    return weight / total


def _sorted_magnitudes(delta: Sequence[float], weights: Sequence[float], tie_tol: float):
    signed = _vector(delta, "delta")
    weight = normalize_weights(weights)
    if signed.size != weight.size:
        raise ValueError("delta and weights must have equal length")
    if not math.isfinite(tie_tol) or tie_tol < 0:
        raise ValueError("tie_tol must be finite and non-negative")
    magnitude = np.abs(signed)
    sign = np.where(signed > tie_tol, 1, np.where(signed < -tie_tol, -1, 0))
    order = np.argsort(magnitude, kind="mergesort")
    return signed, magnitude, sign, weight, order


def _group_moments(magnitude: np.ndarray, sign: np.ndarray, weight: np.ndarray):
    masses = np.zeros(3, dtype=float)
    means = np.zeros(3, dtype=float)
    variances = np.zeros(3, dtype=float)
    second = np.zeros(3, dtype=float)
    for index, state in enumerate((-1, 0, 1)):
        mask = sign == state
        mass = float(np.sum(weight[mask]))
        masses[index] = mass
        if mass <= 0.0:
            continue
        conditional = weight[mask] / mass
        values = magnitude[mask]
        means[index] = float(np.dot(conditional, values))
        variances[index] = float(np.dot(conditional, (values - means[index]) ** 2))
        second[index] = variances[index] + means[index] ** 2
    return masses, means, variances, second


def _mixture_variance(masses: np.ndarray, means: np.ndarray, variances: np.ndarray) -> tuple[float, float]:
    mean = float(np.dot(masses, means))
    variance = float(np.dot(masses, variances + (means - mean) ** 2))
    return mean, max(0.0, variance)


def weighted_quantile_sorted(
    sorted_values: np.ndarray, sorted_weights: np.ndarray, q: float
) -> float:
    """Inverse-CDF weighted quantile matching PL1's convention."""
    if not 0.0 <= q <= 1.0:
        raise ValueError("q must lie in [0,1]")
    if sorted_values.ndim != 1 or sorted_weights.shape != sorted_values.shape:
        raise ValueError("sorted values and weights must be same-length vectors")
    if sorted_values.size == 0 or np.any(sorted_weights < 0.0):
        raise ValueError("sorted support and weights are invalid")
    total = float(np.sum(sorted_weights))
    if not math.isfinite(total) or total <= 0.0:
        raise ValueError("sorted weights must have positive finite mass")
    if q <= 0.0:
        return float(sorted_values[0])
    cumulative = np.cumsum(sorted_weights / total)
    position = int(np.searchsorted(cumulative, q, side="left"))
    return float(sorted_values[min(position, sorted_values.size - 1)])


def crps_sorted(
    sorted_values: Sequence[float], sorted_weights: Sequence[float], truth_magnitude: float
) -> float:
    """Weighted discrete CRPS in O(n) given values sorted once in O(n log n)."""
    values = _vector(sorted_values, "sorted_values")
    weights = normalize_weights(sorted_weights)
    if values.size != weights.size:
        raise ValueError("values and weights must have equal length")
    if np.any(np.diff(values) < 0.0):
        raise ValueError("values must be sorted")
    truth = float(truth_magnitude)
    if not math.isfinite(truth) or truth < 0.0:
        raise ValueError("truth_magnitude must be finite and non-negative")
    first_term = float(np.dot(weights, np.abs(values - truth)))
    before = np.cumsum(weights) - weights
    half_pairwise = float(np.dot(weights * values, 2.0 * before + weights - 1.0))
    return max(0.0, first_term - half_pairwise)


def distribution_summary(
    sorted_magnitude: np.ndarray,
    sorted_weights: np.ndarray,
    truth_magnitude: float,
    *,
    compute_crps: bool = True,
) -> dict[str, float]:
    """Summarize one weighted magnitude distribution using a shared sort order."""
    weights = normalize_weights(sorted_weights)
    if sorted_magnitude.shape != weights.shape:
        raise ValueError("magnitude and weight lengths differ")
    mean = float(np.dot(weights, sorted_magnitude))
    variance = float(np.dot(weights, (sorted_magnitude - mean) ** 2))
    variance = max(0.0, variance)
    q05 = weighted_quantile_sorted(sorted_magnitude, weights, 0.05)
    q25 = weighted_quantile_sorted(sorted_magnitude, weights, 0.25)
    q75 = weighted_quantile_sorted(sorted_magnitude, weights, 0.75)
    q95 = weighted_quantile_sorted(sorted_magnitude, weights, 0.95)
    return {
        "mean": mean,
        "variance": variance,
        "sd": math.sqrt(variance),
        "q05": q05,
        "q25": q25,
        "q75": q75,
        "q95": q95,
        "iqr": q75 - q25,
        "central90_width": q95 - q05,
        "crps": crps_sorted(sorted_magnitude, weights, truth_magnitude) if compute_crps else math.nan,
    }


def _group_crps_components(
    magnitude: np.ndarray, sign: np.ndarray, baseline_weights: np.ndarray,
    truth_magnitude: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Precompute conditional CRPS terms reused by every cue and lambda.

    The first component is E_s|M-truth|. The second is E_st|M_s-M_t|.
    Both are calculated in O(n log n) total using sorted conditional supports;
    the post-cue CRPS is then a three-state mixture calculation in O(1).
    """
    states = (-1, 0, 1)
    first = np.zeros(3, dtype=float)
    cross = np.zeros((3, 3), dtype=float)
    conditional = []
    for index, state in enumerate(states):
        mask = sign == state
        mass = float(np.sum(baseline_weights[mask]))
        if mass <= 0:
            conditional.append(None)
            continue
        conditional_weights = baseline_weights[mask] / mass
        values = magnitude[mask]
        order = np.argsort(values, kind="mergesort")
        values, conditional_weights = values[order], conditional_weights[order]
        first[index] = float(np.dot(conditional_weights, np.abs(values - truth_magnitude)))
        conditional.append((values, conditional_weights))
    for i, left in enumerate(conditional):
        if left is None:
            continue
        x, wx = left
        for j, right in enumerate(conditional):
            if right is None:
                continue
            y, wy = right
            # E|X-Y| with X and Y sorted, using prefix moments of Y.
            prefix_w = np.cumsum(wy)
            prefix_wy = np.cumsum(wy * y)
            left_mass = np.searchsorted(y, x, side="left")
            right_mass = np.searchsorted(y, x, side="right")
            sum_left_w = np.where(left_mass > 0, prefix_w[np.maximum(left_mass - 1, 0)], 0.0)
            sum_left_y = np.where(left_mass > 0, prefix_wy[np.maximum(left_mass - 1, 0)], 0.0)
            sum_tie_w = np.where(right_mass > 0, prefix_w[np.maximum(right_mass - 1, 0)], 0.0) - sum_left_w
            sum_right_w = 1.0 - sum_left_w - sum_tie_w
            sum_right_y = float(prefix_wy[-1]) - sum_left_y - x * sum_tie_w
            expectation = x * sum_left_w - sum_left_y + sum_right_y - x * sum_right_w
            cross[i, j] = float(np.dot(wx, np.maximum(expectation, 0.0)))
    return first, cross


def relative_variance_status(
    baseline_variance: float,
    baseline_second_moment: float,
    baseline_mean: float,
    *,
    multiplier: float = RELATIVE_VARIANCE_EPS_MULTIPLIER,
) -> tuple[bool, float]:
    """Return relative-statistic eligibility and its scale-aware threshold."""
    scale = max(float(baseline_second_moment), float(baseline_mean) ** 2)
    threshold = multiplier * np.finfo(float).eps * scale
    return bool(baseline_variance > threshold and scale > 0.0), threshold


def variance_change_status(
    baseline_variance: float,
    updated_variance: float,
    baseline_second_moment: float,
    updated_second_moment: float,
    baseline_mean: float,
    updated_mean: float,
    *,
    multiplier: float = VARIANCE_CHANGE_EPS_MULTIPLIER,
) -> tuple[str, float]:
    """Classify V0 - V1 as DECREASED, INCREASED, or NEGLIGIBLE."""
    scale = max(
        baseline_second_moment, updated_second_moment,
        baseline_mean ** 2, updated_mean ** 2,
    )
    tolerance = multiplier * np.finfo(float).eps * scale
    change = baseline_variance - updated_variance
    if change > tolerance:
        return "DECREASED", tolerance
    if change < -tolerance:
        return "INCREASED", tolerance
    return "NEGLIGIBLE", tolerance


def analyze_distribution(
    delta: Sequence[float],
    baseline_weights: Sequence[float],
    truth_delta: float,
    *,
    lambdas: Sequence[float] = LAMBDA_GRID,
    tie_tol: float = TIE_TOL,
) -> dict[str, object]:
    """Compute baseline and correct/wrong cue distributions without refitting."""
    truth = float(truth_delta)
    if not math.isfinite(truth):
        raise ValueError("truth_delta must be finite")
    signed, magnitude, sign, weights, order = _sorted_magnitudes(delta, baseline_weights, tie_tol)
    sorted_magnitude = magnitude[order]
    sorted_base_weights = weights[order]
    baseline = distribution_summary(sorted_magnitude, sorted_base_weights, abs(truth))
    masses, means, within_variances, seconds = _group_moments(magnitude, sign, weights)
    baseline_mean, baseline_variance = _mixture_variance(masses, means, within_variances)
    baseline["mean"] = baseline_mean
    baseline["variance"] = baseline_variance
    baseline["sd"] = math.sqrt(baseline_variance)
    baseline_second_moment = float(np.dot(masses, seconds))
    crps_first_by_sign, crps_cross_by_sign = _group_crps_components(
        magnitude, sign, weights, abs(truth)
    )
    baseline_crps_mixture = float(np.dot(masses, crps_first_by_sign) - 0.5 * masses @ crps_cross_by_sign @ masses)
    if not math.isclose(baseline_crps_mixture, baseline["crps"], rel_tol=1e-11, abs_tol=1e-12):
        raise RuntimeError("baseline CRPS mixture identity failed")
    eligible_relative, relative_threshold = relative_variance_status(
        baseline_variance, baseline_second_moment, baseline_mean
    )
    truth_direction = 1 if truth > tie_tol else (-1 if truth < -tie_tol else 0)
    updates: dict[str, dict[str, dict[str, object]]] = {}
    support_mask = weights > 0.0
    for lam_value in lambdas:
        lam = float(lam_value)
        if not math.isfinite(lam) or lam < 0.0:
            raise ValueError("lambdas must be finite and non-negative")
        key = format(lam, ".17g")
        updates[key] = {}
        if truth_direction == 0:
            continue
        for cue_name, cue_direction in (("correct", truth_direction), ("incorrect", -truth_direction)):
            if lam == 0.0:
                current_weights = weights
                summary = dict(baseline)
                current_masses, current_means, current_vars = masses, means, within_variances
            else:
                factor = np.exp(lam * cue_direction * sign)
                raw = weights * factor
                current_weights = raw / float(np.sum(raw))
                factors = np.asarray([math.exp(lam * cue_direction * state) for state in (-1, 0, 1)])
                unnormalized_masses = masses * factors
                current_masses = unnormalized_masses / float(np.sum(unnormalized_masses))
                current_means, current_vars = means, within_variances
                current_mean, current_variance = _mixture_variance(current_masses, current_means, current_vars)
                summary = distribution_summary(
                    sorted_magnitude, current_weights[order], abs(truth), compute_crps=False
                )
                summary["mean"] = current_mean
                summary["variance"] = current_variance
                summary["sd"] = math.sqrt(current_variance)
            summary["crps"] = baseline["crps"] if lam == 0.0 else float(
                np.dot(current_masses, crps_first_by_sign)
                - 0.5 * current_masses @ crps_cross_by_sign @ current_masses
            )
            if not np.array_equal(current_weights > 0.0, support_mask):
                raise RuntimeError("positive exponential reweighting changed weighted support")
            current_second_moment = float(np.dot(current_masses, seconds))
            status, tolerance = variance_change_status(
                baseline_variance, summary["variance"], baseline_second_moment,
                current_second_moment, baseline_mean, summary["mean"],
            )
            summary.update({
                "absolute_variance_reduction": baseline_variance - summary["variance"],
                "relative_variance_reduction": (
                    1.0 - summary["variance"] / baseline_variance if eligible_relative else math.nan
                ),
                "relative_variance_eligible": eligible_relative,
                "relative_variance_threshold": relative_threshold,
                "variance_change_status": status,
                "variance_change_tolerance": tolerance,
                "crps_improvement": baseline["crps"] - summary["crps"],
                "support_count_unchanged": int(np.count_nonzero(current_weights > 0.0)) == int(np.count_nonzero(support_mask)),
            })
            updates[key][cue_name] = summary
    return {
        "truth_direction": truth_direction,
        "truth_magnitude": abs(truth),
        "baseline": baseline,
        "baseline_second_moment": baseline_second_moment,
        "relative_variance_eligible": eligible_relative,
        "relative_variance_threshold": relative_threshold,
        "updates": updates,
        "support_count": int(np.count_nonzero(support_mask)),
    }


def exact_key_match(left: Sequence[tuple], right: Sequence[tuple]) -> list[tuple]:
    left_keys, right_keys = list(left), list(right)
    if len(left_keys) != len(set(left_keys)) or len(right_keys) != len(set(right_keys)):
        raise ValueError("comparison keys must be unique")
    if set(left_keys) != set(right_keys):
        raise ValueError("comparison populations differ")
    return sorted(left_keys)


def crossed_mouse_multiplicities(
    ct2a_mice: Sequence[str], gl261_mice: Sequence[str], rng: np.random.Generator
) -> dict[tuple[str, str], int]:
    """Cross independently sampled mouse identities into pair multiplicities."""
    left, right = tuple(map(str, ct2a_mice)), tuple(map(str, gl261_mice))
    if not left or not right or len(set(left)) != len(left) or len(set(right)) != len(right):
        raise ValueError("mouse identities must be non-empty and unique within cohort")
    sampled_left = rng.choice(left, size=len(left), replace=True)
    sampled_right = rng.choice(right, size=len(right), replace=True)
    left_counts = {mouse: int(np.count_nonzero(sampled_left == mouse)) for mouse in left}
    right_counts = {mouse: int(np.count_nonzero(sampled_right == mouse)) for mouse in right}
    return {
        (mouse_left, mouse_right): left_counts[mouse_left] * right_counts[mouse_right]
        for mouse_left in left for mouse_right in right
        if left_counts[mouse_left] and right_counts[mouse_right]
    }


def prospective_model_specs() -> Mapping[str, tuple[str, ...]]:
    """Frozen predictor combinations; outcomes/usefulness are not predictors."""
    return dict(PREDICTORS)


def _prediction_metrics(observed: np.ndarray, predicted: np.ndarray) -> dict[str, float]:
    residual = observed - predicted
    total = float(np.sum((observed - float(np.mean(observed))) ** 2))
    return {
        "mae": float(np.mean(np.abs(residual))),
        "rmse": float(np.sqrt(np.mean(residual ** 2))),
        "r2": 1.0 - float(np.sum(residual ** 2)) / total if total > 0.0 else math.nan,
    }


def fit_prospective_models(
    rows: Sequence[Mapping[str, object]],
    *,
    target: str,
    split_field: str = "analysis_split",
) -> list[dict[str, object]]:
    """Fit fixed ridge models using only baseline geometry predictors.

    Rows must already be one per reaction and restricted to a common finite
    population across all predictors and both cue targets. Scaling is learned
    from Development only; Confirmation outcomes are read solely for scoring.
    """
    if target not in {"correct_relative_variance_reduction", "incorrect_relative_variance_reduction"}:
        raise ValueError("target must distinguish correct and incorrect cue outcomes")
    import pandas as pd
    from sklearn.linear_model import Ridge
    from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    frame = pd.DataFrame(list(rows))
    required = set(PREDICTORS["all_geometry"]) | {
        target, "correct_relative_variance_reduction",
        "incorrect_relative_variance_reduction", split_field, "reaction_id",
    }
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"prospective model input lacks columns: {sorted(missing)}")
    numeric = list(set(PREDICTORS["all_geometry"]) | {
        "correct_relative_variance_reduction", "incorrect_relative_variance_reduction"
    })
    matrix = frame[numeric].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    if not np.isfinite(matrix).all():
        raise ValueError("prospective rows must use the shared finite population")
    if frame["reaction_id"].duplicated().any():
        raise ValueError("prospective input must contain one row per reaction")
    development = frame[split_field].astype(str).eq("DEVELOPMENT").to_numpy()
    confirmation = frame[split_field].astype(str).eq("CONFIRMATION_HOLDOUT").to_numpy()
    if not development.any() or not confirmation.any() or np.any(development & confirmation):
        raise ValueError("both frozen reaction splits must be present")
    y = pd.to_numeric(frame[target], errors="coerce").to_numpy(float)
    if not np.isfinite(y).all():
        raise ValueError("both cue targets must share the same finite reaction population")
    y_dev, y_conf = y[development], y[confirmation]
    baseline_prediction = np.full(y_conf.shape, float(np.mean(y_dev)))
    baseline_metrics = _prediction_metrics(y_conf, baseline_prediction)
    result = [{
        "target": target,
        "model": "development_mean_baseline",
        "predictors": "none",
        "n_development": int(development.sum()),
        "n_confirmation": int(confirmation.sum()),
        **baseline_metrics,
    }]
    fitted = []
    for model_name, columns in PREDICTORS.items():
        estimator = make_pipeline(StandardScaler(), Ridge(alpha=1.0))
        estimator.fit(frame.loc[development, list(columns)].astype(float), y_dev)
        prediction = estimator.predict(frame.loc[confirmation, list(columns)].astype(float))
        fitted.append({
            "target": target,
            "model": model_name,
            "predictors": "+".join(columns),
            "n_development": int(development.sum()),
            "n_confirmation": int(confirmation.sum()),
            "mae": float(mean_absolute_error(y_conf, prediction)),
            "rmse": float(np.sqrt(mean_squared_error(y_conf, prediction))),
            "r2": float(r2_score(y_conf, prediction)) if np.var(y_conf) > 0 else math.nan,
            "mae_improvement_vs_baseline": baseline_metrics["mae"] - float(mean_absolute_error(y_conf, prediction)),
        })
    eta_baseline = next(row for row in fitted if row["model"] == "eta2")
    for row in fitted:
        row["mae_improvement_vs_eta2"] = eta_baseline["mae"] - row["mae"]
        row["rmse_improvement_vs_eta2"] = eta_baseline["rmse"] - row["rmse"]
    result.extend(fitted)
    return result


def explanatory_partial_rank_association(
    uncertainty_reduction: Sequence[float],
    usefulness: Sequence[float],
    eta2: Sequence[float],
) -> dict[str, float | int]:
    """Descriptive association of eta2 with reduction after usefulness adjustment.

    This intentionally has no train/test interpretation and must not feed the
    prospective model fitting above.
    """
    from scipy.stats import rankdata, spearmanr

    y = np.asarray(uncertainty_reduction, dtype=float)
    u = np.asarray(usefulness, dtype=float)
    e = np.asarray(eta2, dtype=float)
    if y.shape != u.shape or y.shape != e.shape or y.ndim != 1:
        raise ValueError("explanatory inputs must be equal-length vectors")
    mask = np.isfinite(y) & np.isfinite(u) & np.isfinite(e)
    if int(mask.sum()) < 3:
        return {"n": int(mask.sum()), "rho_eta2_raw": math.nan, "rho_eta2_usefulness_adjusted": math.nan}
    raw = float(spearmanr(e[mask], y[mask]).statistic)
    yr, ur, er = (rankdata(values[mask], method="average") for values in (y, u, e))
    design = np.column_stack((np.ones(len(ur)), ur))
    y_residual = yr - design @ np.linalg.lstsq(design, yr, rcond=None)[0]
    e_residual = er - design @ np.linalg.lstsq(design, er, rcond=None)[0]
    denominator = float(np.linalg.norm(y_residual) * np.linalg.norm(e_residual))
    rho = float(np.dot(y_residual, e_residual) / denominator) if denominator else math.nan
    return {"n": int(mask.sum()), "rho_eta2_raw": raw, "rho_eta2_usefulness_adjusted": rho}

