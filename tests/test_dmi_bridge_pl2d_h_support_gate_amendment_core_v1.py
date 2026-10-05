from __future__ import annotations

import copy

import pytest

from scripts import dmi_bridge_pl2d_h_support_gate_amendment_core_v1 as m


def _original_primary():
    return {
        "unit": "reaction",
        "predictor": "mean_sign_magnitude_eta2_across_evaluations",
        "response": "mean_directionally_useful_fraction_across_defined_evaluations",
        "expected_direction": "POSITIVE",
        "minimum_finite_reactions": 700,
        "minimum_spearman_rho": 0.50,
        "context_count": 16,
        "minimum_finite_reactions_per_context": 100,
        "minimum_positive_contexts": 12,
    }


def test_exact_original_identity_and_append_only_status():
    assert m.ORIGINAL_MANIFEST_SHA256 == "a1dd76d6b40e741e07efe64c34bc39eaad289a46bf6c870917a4f59cab688427"
    assert m.ORIGINAL_STATUS == "PL2D_HYPOTHESIS_FROZEN_READY_FOR_CONFIRMATION"
    assert m.STATUS == "PL2D_H_SUPPORT_GATE_AMENDED_READY_FOR_CONFIRMATION"
    assert m.ORIGINAL_ARTIFACT_SHA256["BRIDGEPL2DH_CONFIRMATION_REGISTRY.tsv"] == "02a38afe3cb25b8c85a7774c1f20118ad9601ba74f0df0fbff021ba2494bd2e5"


def test_frozen_context_support_inventory_and_sparsest_context():
    counts = [v for _, _, v in m.FROZEN_DEVELOPMENT_CONTEXT_SUPPORT]
    assert len(counts) == 16
    assert sorted(counts)[:4] == [129, 134, 159, 162]
    assert min(counts) == 129


def test_exact_proportional_support_and_amended_floor():
    counts = [v for _, _, v in m.FROZEN_DEVELOPMENT_CONTEXT_SUPPORT]
    result = m.derive_amended_context_floor(counts)
    assert result["minimum_development_context_support"] == 129
    assert result["proportional_confirmation_support_numerator"] == 43
    assert result["proportional_confirmation_support_denominator"] == 1
    assert result["proportional_confirmation_support"] == 43.0
    assert result["derived_context_floor"] == 30
    assert m.AMENDED_MINIMUM_FINITE_REACTIONS_PER_CONTEXT == 30


def test_original_100_floor_is_incompatible_with_sparsest_development_scaled_support():
    expected = m.expected_proportional_support(129)
    assert float(expected) == 43.0
    assert m.ORIGINAL_MINIMUM_FINITE_REACTIONS_PER_CONTEXT > expected


def test_verify_context_rows_is_order_invariant_and_exact():
    rows = [
        {"algorithm": a, "rna_context_key": r, "finite_reactions": n, "rho": "poison"}
        for a, r, n in reversed(m.FROZEN_DEVELOPMENT_CONTEXT_SUPPORT)
    ]
    counts = m.verify_frozen_development_contexts(rows)
    assert counts == [n for _, _, n in m.FROZEN_DEVELOPMENT_CONTEXT_SUPPORT]
    bad = copy.deepcopy(rows)
    bad[0]["finite_reactions"] = int(bad[0]["finite_reactions"]) + 1
    with pytest.raises(ValueError):
        m.verify_frozen_development_contexts(bad)


def test_amendment_changes_exactly_one_primary_field():
    original = _original_primary()
    amended = m.amended_primary(original)
    differing = {key for key in original if original[key] != amended[key]}
    assert differing == {"minimum_finite_reactions_per_context"}
    assert original["minimum_finite_reactions_per_context"] == 100
    assert amended["minimum_finite_reactions_per_context"] == 30
    for key, value in m.UNCHANGED_PRIMARY.items():
        assert amended[key] == value


def test_amendment_rejects_changes_to_other_primary_fields():
    for key in m.UNCHANGED_PRIMARY:
        original = _original_primary()
        original[key] = "changed"
        with pytest.raises(ValueError):
            m.amended_primary(original)
    original = _original_primary()
    original["minimum_finite_reactions_per_context"] = 99
    with pytest.raises(ValueError):
        m.amended_primary(original)


def test_confirmation_classification_uses_30_context_floor():
    good_counts = [30] * 16
    assert m.classify_confirmation(finite_reactions=700, overall_rho=0.50, context_finite_counts=good_counts, positive_contexts=12) == "PL2D_CONFIRMED"
    low = good_counts.copy()
    low[3] = 29
    assert m.classify_confirmation(finite_reactions=900, overall_rho=0.9, context_finite_counts=low, positive_contexts=16) == "PL2D_INSUFFICIENT_SUPPORT"
    assert m.classify_confirmation(finite_reactions=699, overall_rho=0.9, context_finite_counts=good_counts, positive_contexts=16) == "PL2D_INSUFFICIENT_SUPPORT"
    assert m.classify_confirmation(finite_reactions=900, overall_rho=0.49, context_finite_counts=good_counts, positive_contexts=16) == "PL2D_NOT_CONFIRMED"
    assert m.classify_confirmation(finite_reactions=900, overall_rho=0.8, context_finite_counts=good_counts, positive_contexts=11) == "PL2D_NOT_CONFIRMED"


def test_wrong_context_count_or_nonfinite_rho_is_insufficient_support():
    assert m.classify_confirmation(finite_reactions=900, overall_rho=0.8, context_finite_counts=[100] * 15, positive_contexts=15) == "PL2D_INSUFFICIENT_SUPPORT"
    assert m.classify_confirmation(finite_reactions=900, overall_rho=float("nan"), context_finite_counts=[100] * 16, positive_contexts=16) == "PL2D_INSUFFICIENT_SUPPORT"


def test_amendment_contains_no_confirmation_scientific_values():
    # Foundation freezes identities and development support only; there are no
    # constants for holdout finite counts, holdout rho, or holdout gains.
    names = set(vars(m))
    assert "CONFIRMATION_CONTEXT_SUPPORT" not in names
    assert "CONFIRMATION_RHO" not in names
    assert "CONFIRMATION_GAIN" not in names
