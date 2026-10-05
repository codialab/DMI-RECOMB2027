from __future__ import annotations

import copy

import pytest

from scripts import dmi_bridge_m1_manuscript_core_v1 as core


def _primary():
    return {
        "status": "PL2D_CONFIRMED",
        "finite_reactions": 838,
        "overall_rho": 0.8351435823183292,
        "context_count": 16,
        "positive_contexts": 16,
        "context_min_finite_reactions": 34,
    }


def _supportive():
    return {
        "bootstrap": {"replicates": 5000, "finite_replicates": 5000, "rho_q025": 0.8015629118742916, "rho_median": 0.8354193419051692, "rho_q975": 0.8631275423376382},
        "eta2_vs_mean_information_advantage": {"rho": 0.34542078458252307},
        "partial_eta2_vs_usefulness_controlling_non_tie_coverage": {"rho": 0.8056871576315723},
        "directional_entropy3_vs_usefulness": {"rho": 0.8383836292673364},
        "dominant_sign_mass_vs_usefulness": {"rho": -0.8257657485170178},
        "n_supported_sign_states_vs_usefulness": {"rho": 0.8207118788803361},
        "pathway_sensitivity": {
            "exclude_proximal": {"rho": 0.8386984413243371},
            "exclude_same_subsystem_as_HEX1": {"rho": 0.8378450478802439},
            "exclude_either": {"rho": 0.8386984413243371},
        },
    }


def _qc():
    return {
        "counts": {"confirmation_cases": 873620, "non_tie_cases": 359099, "truth_tie_cases": 514521, "confirmation_evaluation_reaction_rows": 386650},
        "primary_result_frozen_before_supportive": True,
        "no_quantile_alias_weighting": True,
        "maximum_random_sign_identity_deviation": 0.0,
        "maximum_information_advantage_identity_deviation": 0.0,
    }


def test_source_hashes_are_frozen():
    assert core.PL2B_MANIFEST_SHA256.startswith("7d9ca97a")
    assert core.PL2C_MANIFEST_SHA256.startswith("6abcf149")
    assert core.PL2D_MANIFEST_SHA256.startswith("9f9eb730")
    assert len(core.PL2D_PRIMARY_SHA256) == 64


def test_development_snapshot_frozen_values():
    core.validate_development_snapshot({
        "finite_primary_reactions": 2455,
        "eta2_usefulness_rho": 0.8570601823118907,
        "positive_contexts": 16,
    })
    with pytest.raises(ValueError):
        core.validate_development_snapshot({"finite_primary_reactions": 2455, "eta2_usefulness_rho": 0.90, "positive_contexts": 16})


def test_confirmation_primary_exact():
    facts = core.validate_confirmation_primary(_primary())
    assert facts.finite_reactions == 838
    assert facts.positive_contexts == 16
    assert facts.minimum_context_n == 34


def test_confirmation_not_confirmed_is_rejected_for_m1():
    x = _primary(); x["status"] = "PL2D_NOT_CONFIRMED"
    with pytest.raises(ValueError):
        core.validate_confirmation_primary(x)


def test_supportive_values_are_frozen():
    core.validate_supportive(_supportive())
    x = copy.deepcopy(_supportive())
    x["eta2_vs_mean_information_advantage"]["rho"] = 0.7
    with pytest.raises(ValueError):
        core.validate_supportive(x)


def test_qc_counts_and_no_alias_weighting():
    core.validate_confirmation_qc(_qc())
    x = copy.deepcopy(_qc()); x["no_quantile_alias_weighting"] = False
    with pytest.raises(ValueError):
        core.validate_confirmation_qc(x)


def test_pl2b_manuscript_rounding_and_random_control():
    rows = [
        {"reliability_q": 0.0, "improved_fraction_of_non_ties": 0.11234, "tied_fraction_of_non_ties": 0.71451},
        {"reliability_q": 0.5, "improved_fraction_of_non_ties": 0.15414, "tied_fraction_of_non_ties": 0.72352},
        {"reliability_q": 1.0, "improved_fraction_of_non_ties": 0.16954, "tied_fraction_of_non_ties": 0.71599},
    ]
    # Synthetic full-precision values sum to one and reproduce the frozen manuscript rounding.
    for row in rows:
        row["harmed_fraction_of_non_ties"] = 1.0 - row["improved_fraction_of_non_ties"] - row["tied_fraction_of_non_ties"]
    core.validate_pl2b_display_rows(rows)


def test_pl2b_missing_random_control_fails():
    rows = [
        {"reliability_q": 0.0, "improved_fraction_of_non_ties": 0.1123, "tied_fraction_of_non_ties": 0.7145, "harmed_fraction_of_non_ties": 0.1732},
        {"reliability_q": 1.0, "improved_fraction_of_non_ties": 0.1695, "tied_fraction_of_non_ties": 0.7160, "harmed_fraction_of_non_ties": 0.1145},
    ]
    with pytest.raises(ValueError):
        core.validate_pl2b_display_rows(rows)


def test_claim_ledger_contains_primary_supportive_and_limitation():
    claims = core.claim_specs()
    ids = {x["claim_id"] for x in claims}
    assert "C2_GEOMETRY_PREDICTS_DIRECTIONAL_USEFULNESS" in ids
    assert "C3_OCCURRENCE_STRONGER_THAN_GAIN_MAGNITUDE" in ids
    assert "L1_DIRECTIONAL_TRUTH_COVERAGE" in ids
    assert any(x["claim_level"] == "LIMITATION" for x in claims)


def test_claim_text_does_not_overclaim_gain_magnitude():
    claim = next(x for x in core.claim_specs() if x["claim_id"] == "C3_OCCURRENCE_STRONGER_THAN_GAIN_MAGNITUDE")
    assert "more strongly" in claim["statement"]
    assert "Do not claim accurate prediction" in claim["boundary"]


def test_figure_panel_plan_is_fixed_and_complete():
    panels = core.figure_panels()
    assert [x["panel"] for x in panels] == ["A", "B", "C", "D", "E"]
    assert "random-sign" in panels[1]["message"]
    assert "never pool" in panels[2]["message"]
    assert "all 16" in panels[3]["message"]


def test_terminal_status_is_manuscript_evidence_only():
    assert core.TERMINAL_STATUS == "BRIDGE_M1_MANUSCRIPT_EVIDENCE_FROZEN"
