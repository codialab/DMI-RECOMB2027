#!/usr/bin/env python3
"""Frozen manuscript-facing primitives for DMI-BRIDGE-M2."""

from __future__ import annotations

import math
from typing import Mapping, Sequence

SCHEMA = "bridge.m2.a2_robustness_evidence.v1"
TERMINAL_STATUS = "BRIDGE_M2_A2_ROBUSTNESS_EVIDENCE_FROZEN"

M1_MANIFEST_SHA256 = "201ba9f4fc2f74290d991cd294aec558c5e5dd23282877e871a5c0669a44cdde"
A20_MANIFEST_SHA256 = "01f8004a0bf9d5358d7ce7281e827cc1a1fbe00a6750a5c063a1edad2eedf380"
A21_MANIFEST_SHA256 = "db9948d7b99447a695cf8aa21429dde6d642e1fd0632bcce52643deac59e6feb"
A22_MANIFEST_SHA256 = "15b237151005538e5077a352868f873ef00466f4f5ea49307b781a8a8cd0e7cb"
A23_MANIFEST_SHA256 = "b8e57074f3df23fe7708145bc01f51a19f946545b1c16e61628c2678dbcabd18"

A23_ARTIFACT_SHA256 = {
    "BRIDGEA23_QC.json": "ed502de0e91e793d4c0ed3c7a0237355ae30a0bf69d81ccde2e1a818e358ee20",
    "BRIDGEA23_UTILITY_COMPARISON.tsv": "a59b1870d50efc87c75d3be32c5e890f0b63413e38c87a47adb65dc99f00f354",
    "BRIDGEA23_GEOMETRY_COMPARISON.tsv": "e9e0561248a231d5f9480365888a76cd8c99f6d06dd04ab4d92141601fa21f2f",
    "BRIDGEA23_GEOMETRY_UTILITY_PERSISTENCE.tsv": "33af7c3b59a32872aee0070af10f35b72f974537bb6e1745dcc2c5fbd42afc03",
    "BRIDGEA23_CONTROL_SUMMARY.tsv": "e70e9361e68cfd3e8d72eb27095a0781eb564626913951b3bbe5631d83f45cb9",
    "BRIDGEA23_SENSITIVITY_SUMMARY.tsv": "2d33a5358ad38741660585ffbb52f68636907e825467cedaab46ad7df5a6497f",
}

EXPECTED_TRUTH_PAIRS = 836
EXPECTED_REACTIONS = 4179
EXPECTED_CASE_KEYS = 3493644
EXPECTED_NON_TIES = 1440313
EXPECTED_TRUTH_TIES = 2053331
EXPECTED_GEOMETRY_KEYS = 1671600
EXPECTED_ETA2_FINITE_PAIRED = 573882

EXPECTED_ALL_UTILITY = {
    "a1_correct_mean_gain": 0.014234484268356776,
    "a2_correct_mean_gain": 0.022801865806147577,
    "paired_correct_delta": 0.0085673815377913642,
    "a1_specificity_mean": 0.03309389154965113,
    "a2_specificity_mean": 0.054535148930311576,
    "paired_specificity_delta": 0.021441257380664738,
}

EXPECTED_ALL_CONTROLS = {
    "correct_q1": (0.014234484268356776, 0.022801865806147577),
    "wrong_q0": (-0.018859407281294222, -0.031733283124166602),
    "analytic_random_q0_5": (-0.0023124615064686166, -0.0044657086590093476),
}

EXPECTED_PERSISTENCE = {
    "DEVELOPMENT": {
        "a1_rho": 0.8571160663367106,
        "a2_rho": 0.8339440056039364,
        "a1_finite": 2454,
        "a2_finite": 2485,
    },
    "CONFIRMATION_HOLDOUT": {
        "a1_rho": 0.8351435823183291,
        "a2_rho": 0.8493657186460156,
        "a1_finite": 838,
        "a2_finite": 843,
    },
}

EXPECTED_A2_BOOTSTRAP = {
    "replicates": 5000,
    "finite_replicates": 5000,
    "rho_median": 0.8488457925736193,
    "rho_q025": 0.823810132333918,
    "rho_q975": 0.8714723352450245,
}


def _close(actual: object, expected: float, tol: float = 1e-12) -> bool:
    try:
        value = float(actual)
    except (TypeError, ValueError):
        return False
    return math.isfinite(value) and abs(value - expected) <= tol


def validate_a23_manifest(manifest: Mapping[str, object]) -> None:
    if manifest.get("status") != "BRIDGE_A23_A1_A2_SYNTHESIS_FROZEN":
        raise ValueError("A2.3 terminal status mismatch")
    if bool(manifest.get("a2_g_used")):
        raise ValueError("A2-G must not be used")
    if bool(manifest.get("new_production_performed")):
        raise ValueError("A2.3 unexpectedly performed new production")
    pop = manifest["case_population"]
    if int(pop["truth_pairs"]) != EXPECTED_TRUTH_PAIRS:
        raise ValueError("truth-pair count mismatch")
    if int(pop["reactions"]) != EXPECTED_REACTIONS:
        raise ValueError("reaction count mismatch")
    if int(pop["case_keys"]) != EXPECTED_CASE_KEYS:
        raise ValueError("case-key count mismatch")


def validate_utility_all(row: Mapping[str, object]) -> None:
    if row.get("scope") != "ALL":
        raise ValueError("expected ALL utility row")
    if int(row["all_case_keys"]) != EXPECTED_CASE_KEYS:
        raise ValueError("ALL case-key count mismatch")
    if int(row["non_tie_pairs"]) != EXPECTED_NON_TIES:
        raise ValueError("non-tie count mismatch")
    if int(row["truth_ties"]) != EXPECTED_TRUTH_TIES:
        raise ValueError("truth-tie count mismatch")
    for field, expected in EXPECTED_ALL_UTILITY.items():
        if not _close(row[field], expected):
            raise ValueError(f"unexpected ALL utility field {field}")


def validate_control_rows(rows: Sequence[Mapping[str, object]]) -> None:
    indexed = {str(r["arm"]): r for r in rows if r.get("scope") == "ALL"}
    if set(indexed) != set(EXPECTED_ALL_CONTROLS):
        raise ValueError("ALL control-arm inventory mismatch")
    for arm, (a1, a2) in EXPECTED_ALL_CONTROLS.items():
        row = indexed[arm]
        if int(row["denominator"]) != EXPECTED_NON_TIES:
            raise ValueError(f"control denominator mismatch for {arm}")
        if not _close(row["a1_mean_gain"], a1):
            raise ValueError(f"A1 gain mismatch for {arm}")
        if not _close(row["a2_mean_gain"], a2):
            raise ValueError(f"A2 gain mismatch for {arm}")


def validate_persistence_rows(rows: Sequence[Mapping[str, object]]) -> None:
    indexed = {
        str(r["split"]): r
        for r in rows
        if r.get("level") == "PL2D_REACTION"
        and r.get("context") == "ALL"
        and r.get("feature") == "sign_magnitude_eta2"
    }
    if set(indexed) != set(EXPECTED_PERSISTENCE):
        raise ValueError("primary persistence row inventory mismatch")
    for split, expected in EXPECTED_PERSISTENCE.items():
        row = indexed[split]
        for field in ("a1_rho", "a2_rho"):
            if not _close(row[field], expected[field]):
                raise ValueError(f"{split} {field} mismatch")
        for field in ("a1_finite", "a2_finite"):
            if int(row[field]) != int(expected[field]):
                raise ValueError(f"{split} {field} mismatch")


def validate_qc(qc: Mapping[str, object]) -> None:
    if qc.get("status") != "PASS":
        raise ValueError("A2.3 QC is not PASS")
    if int(qc["matched_truth_pairs"]) != EXPECTED_TRUTH_PAIRS:
        raise ValueError("QC truth-pair count mismatch")
    if int(qc["reactions"]) != EXPECTED_REACTIONS:
        raise ValueError("QC reaction count mismatch")
    if int(qc["matched_case_keys"]) != EXPECTED_CASE_KEYS:
        raise ValueError("QC case-key count mismatch")
    if int(qc["geometry_keys"]) != EXPECTED_GEOMETRY_KEYS:
        raise ValueError("QC geometry-key count mismatch")
    b = qc["bootstrap"]
    for field, expected in EXPECTED_A2_BOOTSTRAP.items():
        if isinstance(expected, int):
            if int(b[field]) != expected:
                raise ValueError(f"bootstrap {field} mismatch")
        elif not _close(b[field], expected):
            raise ValueError(f"bootstrap {field} mismatch")


def claim_specs() -> list[dict[str, str]]:
    return [
        {
            "claim_id": "R1_MATCHED_CORRECT_DIRECTION_UTILITY",
            "claim_level": "ROBUSTNESS_PRIMARY",
            "statement": "In the exact matched post-freeze population, mean correct-direction gain is larger under the A2-L dual-anchor baseline than under A1.",
            "boundary": "Descriptive matched robustness only; do not call this independent confirmation or a newly prespecified hypothesis test.",
        },
        {
            "claim_id": "R2_DIRECTIONAL_SPECIFICITY",
            "claim_level": "ROBUSTNESS_PRIMARY",
            "statement": "The A2-L increase in correct-direction gain is accompanied by a larger correct-minus-wrong contrast and more harmful wrong-direction control.",
            "boundary": "Do not describe A2 as a generic regularization benefit; q=0.5 is a random-sign control with no directional information on average.",
        },
        {
            "claim_id": "R3_GEOMETRY_UTILITY_PERSISTENCE",
            "claim_level": "ROBUSTNESS_PRIMARY",
            "statement": "The frozen geometry-to-directional-usefulness association remains strong under A2-L in both the original development and confirmation reaction splits.",
            "boundary": "Do not claim the A1/A2 rho difference is significant or that A2 independently confirms the hypothesis.",
        },
        {
            "claim_id": "R4_RESIDUAL_GEOMETRY_SHIFT",
            "claim_level": "SUPPORTIVE",
            "statement": "Adding the qualified lactate anchor changes several residual-geometry descriptors while leaving substantial directional structure and finite-support limitations.",
            "boundary": "Do not collapse mixed descriptor changes into a universal claim that A2 is simply more constrained or more ambiguous.",
        },
        {
            "claim_id": "L1_POST_FREEZE_STATUS",
            "claim_level": "LIMITATION",
            "statement": "A2.0-A2.3 were performed after the original PL2D confirmation was opened and therefore constitute robustness/generalization evidence rather than untouched confirmation.",
            "boundary": "Never label A2 as an independent confirmation cohort.",
        },
        {
            "claim_id": "L2_A2L_ONLY",
            "claim_level": "LIMITATION",
            "statement": "Only the lactate second-anchor arm A2-L qualified; A2-G failed the quantitative-operator provenance gate.",
            "boundary": "Do not generalize the A2-L result to arbitrary second strong anchors.",
        },
    ]
