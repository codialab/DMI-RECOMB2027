import math

import numpy as np
import pandas as pd

from scripts import dmi_bridge_pl2d_confirmation_core_v1 as core
from scripts import dmi_bridge_pl2d_hypothesis_freeze_core_v1 as hcore


def _contexts(count=30, rho=0.2):
    return [count] * 16, [rho] * 16


def test_frozen_cardinalities():
    assert core.EXPECTED_EVALUATION_REACTION_ROWS == 386_650
    assert core.EXPECTED_CONFIRMATION_CASE_ROWS == 873_620
    assert core.EXPECTED_DEVELOPMENT_CASE_ROWS_TO_SKIP == 2_620_860


def test_primary_gate_confirmed_at_boundaries():
    counts, rhos = _contexts(30, 0.1)
    result = core.classify_confirmation(
        finite_reactions=700, overall_rho=0.50,
        context_finite_counts=counts, context_rhos=rhos,
    )
    assert result.status == "PL2D_CONFIRMED"
    assert result.positive_contexts == 16


def test_primary_gate_29_is_insufficient():
    counts, rhos = _contexts(30, 0.1)
    counts[7] = 29
    result = core.classify_confirmation(
        finite_reactions=1000, overall_rho=0.9,
        context_finite_counts=counts, context_rhos=rhos,
    )
    assert result.status == "PL2D_INSUFFICIENT_SUPPORT"


def test_primary_gate_699_is_insufficient():
    counts, rhos = _contexts()
    assert core.classify_confirmation(
        finite_reactions=699, overall_rho=0.9,
        context_finite_counts=counts, context_rhos=rhos,
    ).status == "PL2D_INSUFFICIENT_SUPPORT"


def test_primary_gate_nonfinite_rho_is_insufficient():
    counts, rhos = _contexts()
    assert core.classify_confirmation(
        finite_reactions=900, overall_rho=math.nan,
        context_finite_counts=counts, context_rhos=rhos,
    ).status == "PL2D_INSUFFICIENT_SUPPORT"


def test_primary_gate_nonfinite_context_rho_is_insufficient():
    counts, rhos = _contexts()
    rhos[4] = math.nan
    assert core.classify_confirmation(
        finite_reactions=900, overall_rho=0.8,
        context_finite_counts=counts, context_rhos=rhos,
    ).status == "PL2D_INSUFFICIENT_SUPPORT"


def test_primary_gate_effect_failure_is_not_confirmed():
    counts, rhos = _contexts()
    assert core.classify_confirmation(
        finite_reactions=900, overall_rho=0.499999,
        context_finite_counts=counts, context_rhos=rhos,
    ).status == "PL2D_NOT_CONFIRMED"


def test_primary_gate_requires_12_strictly_positive_contexts():
    counts = [40] * 16
    rhos = [0.1] * 11 + [0.0] * 5
    assert core.classify_confirmation(
        finite_reactions=900, overall_rho=0.8,
        context_finite_counts=counts, context_rhos=rhos,
    ).status == "PL2D_NOT_CONFIRMED"
    rhos[11] = 1e-300
    assert core.classify_confirmation(
        finite_reactions=900, overall_rho=0.8,
        context_finite_counts=counts, context_rhos=rhos,
    ).status == "PL2D_CONFIRMED"


def test_primary_association_uses_finite_reactions_only():
    rows = [
        {"mean_sign_magnitude_eta2": 1.0, "mean_directionally_useful_fraction": 1.0},
        {"mean_sign_magnitude_eta2": 2.0, "mean_directionally_useful_fraction": 2.0},
        {"mean_sign_magnitude_eta2": None, "mean_directionally_useful_fraction": 3.0},
    ]
    n, rho = core.primary_association(rows)
    assert n == 2
    assert math.isclose(rho, 1.0, rel_tol=0.0, abs_tol=1e-15)


def test_reaction_aggregation_authority_is_pandas_group_mean():
    frame = pd.DataFrame({
        "reaction_id": ["R1", "R1", "R2", "R2"],
        "sign_magnitude_eta2": [0.1, 0.2, 0.2, 0.4],
        "directionally_useful_fraction": [0.0, 1.0, 0.0, 1.0],
        "non_tie_pair_fraction": [0.5, 0.5, 1.0, 1.0],
        "mean_information_advantage": [1.0, 3.0, 2.0, 4.0],
        "directional_entropy3": [0.0, 0.1, 0.2, 0.3],
        "dominant_sign_mass": [1.0, 0.9, 0.8, 0.7],
        "n_supported_sign_states": [1.0, 2.0, 2.0, 3.0],
    })
    rows = hcore.aggregate_reaction_frame(frame)
    by_id = {r["reaction_id"]: r for r in rows}
    assert by_id["R1"]["mean_sign_magnitude_eta2"] == np.mean([0.1, 0.2])
    assert by_id["R2"]["mean_directionally_useful_fraction"] == 0.5


def _bootstrap_rows():
    rows = []
    for i in range(core.EXPECTED_CONFIRMATION_REACTIONS):
        rows.append({
            "reaction_id": f"R{i:04d}",
            "mean_sign_magnitude_eta2": float(i),
            "mean_directionally_useful_fraction": float(i % 17),
        })
    return rows


def test_bootstrap_is_deterministic():
    rows = _bootstrap_rows()
    a = core.bootstrap_primary(rows, replicates=25)
    b = core.bootstrap_primary(rows, replicates=25)
    np.testing.assert_array_equal(a["rho"], b["rho"])
    np.testing.assert_array_equal(a["finite_pair_counts"], b["finite_pair_counts"])
    assert a["rng"] == "numpy.random.Generator(PCG64)"


def test_bootstrap_rejects_wrong_population_size():
    try:
        core.bootstrap_primary(_bootstrap_rows()[:-1], replicates=2)
    except ValueError as exc:
        assert "1045" in str(exc)
    else:
        raise AssertionError("wrong confirmation population size accepted")


def test_contract_constants_freeze_amended_gate():
    frozen = core.contract_constants()
    assert frozen["primary_minimum_finite_reactions"] == 700
    assert frozen["primary_minimum_finite_reactions_per_context"] == 30
    assert frozen["primary_minimum_rho"] == 0.50
    assert frozen["primary_minimum_positive_contexts"] == 12
    assert frozen["bootstrap_replicates"] == 5000
    assert frozen["bootstrap_seed"] == 20260929
