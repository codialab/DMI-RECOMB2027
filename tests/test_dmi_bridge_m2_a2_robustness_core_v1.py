import pytest

from scripts.dmi_bridge_m2_a2_robustness_core_v1 import (
    EXPECTED_CASE_KEYS,
    EXPECTED_GEOMETRY_KEYS,
    EXPECTED_NON_TIES,
    EXPECTED_REACTIONS,
    EXPECTED_TRUTH_PAIRS,
    EXPECTED_TRUTH_TIES,
    claim_specs,
    validate_a23_manifest,
    validate_control_rows,
    validate_persistence_rows,
    validate_qc,
    validate_utility_all,
)


def test_validate_a23_manifest():
    validate_a23_manifest({
        "status": "BRIDGE_A23_A1_A2_SYNTHESIS_FROZEN",
        "a2_g_used": False,
        "new_production_performed": False,
        "case_population": {
            "truth_pairs": EXPECTED_TRUTH_PAIRS,
            "reactions": EXPECTED_REACTIONS,
            "case_keys": EXPECTED_CASE_KEYS,
        },
    })


def test_a23_manifest_rejects_a2g():
    with pytest.raises(ValueError):
        validate_a23_manifest({
            "status": "BRIDGE_A23_A1_A2_SYNTHESIS_FROZEN",
            "a2_g_used": True,
            "new_production_performed": False,
            "case_population": {
                "truth_pairs": EXPECTED_TRUTH_PAIRS,
                "reactions": EXPECTED_REACTIONS,
                "case_keys": EXPECTED_CASE_KEYS,
            },
        })


def test_validate_utility_all():
    validate_utility_all({
        "scope": "ALL",
        "all_case_keys": EXPECTED_CASE_KEYS,
        "non_tie_pairs": EXPECTED_NON_TIES,
        "truth_ties": EXPECTED_TRUTH_TIES,
        "a1_correct_mean_gain": 0.014234484268356776,
        "a2_correct_mean_gain": 0.022801865806147577,
        "paired_correct_delta": 0.0085673815377913642,
        "a1_specificity_mean": 0.03309389154965113,
        "a2_specificity_mean": 0.054535148930311576,
        "paired_specificity_delta": 0.021441257380664738,
    })


def test_validate_control_rows():
    validate_control_rows([
        {"scope": "ALL", "arm": "correct_q1", "denominator": EXPECTED_NON_TIES,
         "a1_mean_gain": 0.014234484268356776, "a2_mean_gain": 0.022801865806147577},
        {"scope": "ALL", "arm": "wrong_q0", "denominator": EXPECTED_NON_TIES,
         "a1_mean_gain": -0.018859407281294222, "a2_mean_gain": -0.031733283124166602},
        {"scope": "ALL", "arm": "analytic_random_q0_5", "denominator": EXPECTED_NON_TIES,
         "a1_mean_gain": -0.0023124615064686166, "a2_mean_gain": -0.0044657086590093476},
    ])


def test_validate_persistence_rows():
    validate_persistence_rows([
        {"level": "PL2D_REACTION", "split": "DEVELOPMENT", "context": "ALL",
         "feature": "sign_magnitude_eta2", "a1_rho": 0.8571160663367106,
         "a2_rho": 0.8339440056039364, "a1_finite": 2454, "a2_finite": 2485},
        {"level": "PL2D_REACTION", "split": "CONFIRMATION_HOLDOUT", "context": "ALL",
         "feature": "sign_magnitude_eta2", "a1_rho": 0.8351435823183291,
         "a2_rho": 0.8493657186460156, "a1_finite": 838, "a2_finite": 843},
    ])


def test_validate_qc():
    validate_qc({
        "status": "PASS",
        "matched_truth_pairs": EXPECTED_TRUTH_PAIRS,
        "reactions": EXPECTED_REACTIONS,
        "matched_case_keys": EXPECTED_CASE_KEYS,
        "geometry_keys": EXPECTED_GEOMETRY_KEYS,
        "bootstrap": {
            "replicates": 5000,
            "finite_replicates": 5000,
            "rho_median": 0.8488457925736193,
            "rho_q025": 0.823810132333918,
            "rho_q975": 0.8714723352450245,
        },
    })


def test_claim_ledger_has_primary_and_limitations():
    specs = claim_specs()
    ids = {r["claim_id"] for r in specs}
    assert "R1_MATCHED_CORRECT_DIRECTION_UTILITY" in ids
    assert "R3_GEOMETRY_UTILITY_PERSISTENCE" in ids
    assert "L1_POST_FREEZE_STATUS" in ids
    assert "L2_A2L_ONLY" in ids
