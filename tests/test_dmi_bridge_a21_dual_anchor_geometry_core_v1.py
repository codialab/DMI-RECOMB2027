import math

import numpy as np

from scripts.dmi_bridge_a21_dual_anchor_geometry_core_v1 import (
    A1_COMPARABLE_FEATURES,
    analysis_contract,
    lactate_production,
    target_landscape_metrics_a21,
)
from scripts.dmi_bridge_pl1_predictability_core_v1 import target_landscape_metrics


def _toy():
    hex1_c = np.array([1.0, 2.0, 3.0])
    hex1_g = np.array([0.5, 1.5, 2.5])
    lac_c = np.array([0.0, 2.0, 4.0])
    lac_g = np.array([4.0, 2.0, 0.0])
    target_c = np.array([1.0, 4.0, 9.0])
    target_g = np.array([0.5, 2.0, 8.0])
    wc = np.array([0.2, 0.3, 0.5])
    wg = np.array([0.5, 0.3, 0.2])
    return hex1_c, hex1_g, lac_c, lac_g, target_c, target_g, wc, wg


def test_lactate_production_semantics():
    got = lactate_production([3.0, 0.0, -2.0, -5.0])
    np.testing.assert_allclose(got, [0.0, 0.0, 2.0, 5.0])


def test_a21_preserves_every_pl1_field_value():
    h_c, h_g, l_c, l_g, b_c, b_g, wc, wg = _toy()
    base = target_landscape_metrics(
        anchor_ct2a=h_c,
        anchor_gl261=h_g,
        target_ct2a=b_c,
        target_gl261=b_g,
        strong_weights_ct2a=wc,
        strong_weights_gl261=wg,
    )
    got = target_landscape_metrics_a21(
        hex1_ct2a=h_c,
        hex1_gl261=h_g,
        lactate_ct2a=l_c,
        lactate_gl261=l_g,
        target_ct2a=b_c,
        target_gl261=b_g,
        strong_weights_ct2a=wc,
        strong_weights_gl261=wg,
    )
    for key, value in base.items():
        if isinstance(value, float) and math.isnan(value):
            assert math.isnan(got[key])
        else:
            assert got[key] == value


def test_hex1_aliases_are_exact():
    h_c, h_g, l_c, l_g, b_c, b_g, wc, wg = _toy()
    got = target_landscape_metrics_a21(
        hex1_ct2a=h_c,
        hex1_gl261=h_g,
        lactate_ct2a=l_c,
        lactate_gl261=l_g,
        target_ct2a=b_c,
        target_gl261=b_g,
        strong_weights_ct2a=wc,
        strong_weights_gl261=wg,
    )
    assert got["delta_hex1_mean"] == got["delta_a_mean"]
    assert got["delta_hex1_sd"] == got["delta_a_sd"]
    assert got["delta_hex1_target_correlation"] == got["delta_anchor_target_correlation"]
    assert got["abs_delta_hex1_target_correlation"] == got["abs_delta_anchor_target_correlation"]


def test_lactate_coupling_fields_are_finite_for_nondegenerate_toy():
    h_c, h_g, l_c, l_g, b_c, b_g, wc, wg = _toy()
    got = target_landscape_metrics_a21(
        hex1_ct2a=h_c,
        hex1_gl261=h_g,
        lactate_ct2a=l_c,
        lactate_gl261=l_g,
        target_ct2a=b_c,
        target_gl261=b_g,
        strong_weights_ct2a=wc,
        strong_weights_gl261=wg,
    )
    assert math.isfinite(got["delta_lactate_target_correlation"])
    assert math.isfinite(got["abs_delta_lactate_target_correlation"])
    assert got["lactate_target_coupling_status"] == "OK"


def test_degenerate_lactate_is_retained_explicitly():
    h_c, h_g, _, _, b_c, b_g, wc, wg = _toy()
    got = target_landscape_metrics_a21(
        hex1_ct2a=h_c,
        hex1_gl261=h_g,
        lactate_ct2a=[1.0, 1.0, 1.0],
        lactate_gl261=[1.0, 1.0, 1.0],
        target_ct2a=b_c,
        target_gl261=b_g,
        strong_weights_ct2a=wc,
        strong_weights_gl261=wg,
    )
    assert got["lactate_delta_degenerate"] is True
    assert math.isnan(got["delta_lactate_target_correlation"])
    assert got["lactate_target_coupling_status"] == "DEGENERATE"


def test_contract_freezes_a2_l_only_and_no_composite_score():
    contract = analysis_contract()
    assert contract["qualified_input_arm"] == "A2-L"
    assert contract["a2_g_used"] is False
    assert contract["strong_anchor_set"] == ["HEX1", "LDH_L"]
    assert "composite_dual_anchor_bridgeability_score" in contract["forbidden"]
    assert tuple(contract["a1_comparable_features"]) == A1_COMPARABLE_FEATURES
