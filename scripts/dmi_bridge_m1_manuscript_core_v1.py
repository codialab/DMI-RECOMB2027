#!/usr/bin/env python3
"""Frozen manuscript-facing primitives for DMI-BRIDGE-M1."""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Mapping, Sequence

SCHEMA = "bridge.m1.manuscript_evidence.v1"
TERMINAL_STATUS = "BRIDGE_M1_MANUSCRIPT_EVIDENCE_FROZEN"

PL2B_MANIFEST_SHA256 = "7d9ca97a3f4ea5aa8b0b983e09ee83fc3c82f326d1b1cbf9c6137cfadf25539d"
PL2B_CONTEXT_SUMMARY_SHA256 = "d59b0eec6956081a760a52359ca489e461d4b8dbeac784a223c78dc3d20a9a01"
PL2C_MANIFEST_SHA256 = "6abcf149cd6adb294c0931585cec28f992b453451a5b5f03e10b3ffb2f9b9c39"
PL2DH_MANIFEST_SHA256 = "a1dd76d6b40e741e07efe64c34bc39eaad289a46bf6c870917a4f59cab688427"
PL2DH_DEVELOPMENT_SNAPSHOT_SHA256 = "ebd60bbdd480ec83fe233e068a398776b4bd4ff964cc4bc6bf8c615c66cc28b8"
PL2DHA_MANIFEST_SHA256 = "fb7df9e208da2e6c9841d43a3b572038e59943a4221535705eef9f5dc0ed6a3a"
PL2D_MANIFEST_SHA256 = "9f9eb730473e7ad209c1d5ec3340808ab614087b8ae1121d929eff61be12b699"
PL2D_PRIMARY_SHA256 = "a01ca592b819ab2e6b4463bcbb91ba3bc113bd80dbde2c7aa96b505521bfa1b8"
PL2D_SUPPORTIVE_SHA256 = "fbd1e8398a248eeb043e9e985c09dbdcd5d75f067168d8d35b8d18cbec2735b0"
PL2D_CONTEXT_SHA256 = "c063d9809ee0999ce030f7b557fc8064b53348ccda551f3175dd71ef477c744b"
PL2D_QC_SHA256 = "108befc356110ced2ec19ff0fcec6fc01ce38d3d6dd5dc584f2b9f11915936ea"

EXPECTED_DEVELOPMENT_FINITE = 2455
EXPECTED_DEVELOPMENT_RHO = 0.8570601823118907
EXPECTED_CONFIRMATION_FINITE = 838
EXPECTED_CONFIRMATION_RHO = 0.8351435823183292
EXPECTED_CONFIRMATION_CONTEXTS = 16
EXPECTED_CONFIRMATION_POSITIVE_CONTEXTS = 16
EXPECTED_CONFIRMATION_MIN_CONTEXT_N = 34
EXPECTED_BOOTSTRAP_Q025 = 0.8015629118742916
EXPECTED_BOOTSTRAP_MEDIAN = 0.8354193419051692
EXPECTED_BOOTSTRAP_Q975 = 0.8631275423376382
EXPECTED_INFO_ADVANTAGE_RHO = 0.34542078458252307
EXPECTED_PARTIAL_RHO = 0.8056871576315723
EXPECTED_ENTROPY_RHO = 0.8383836292673364
EXPECTED_DOMINANT_SIGN_RHO = -0.8257657485170178
EXPECTED_SIGN_STATES_RHO = 0.8207118788803361
EXPECTED_PATHWAY_EXCLUDE_PROXIMAL_RHO = 0.8386984413243371
EXPECTED_PATHWAY_EXCLUDE_SAME_SUBSYSTEM_RHO = 0.8378450478802439
EXPECTED_PATHWAY_EXCLUDE_EITHER_RHO = 0.8386984413243371
EXPECTED_CONFIRMATION_CASES = 873620
EXPECTED_NON_TIE_CASES = 359099
EXPECTED_TRUTH_TIE_CASES = 514521
EXPECTED_EVALUATION_REACTION_ROWS = 386650
EXPECTED_DEFINED_ROWS = 164690
EXPECTED_NO_DIRECTIONAL_ROWS = 221960

# Manuscript-rounded PL2B primary-lambda overall reliability rows.
EXPECTED_PL2B_DISPLAY = {
    0.00: (0.1123, 0.7145, 0.1732),
    0.50: (0.1541, 0.7235, 0.1223),
    1.00: (0.1695, 0.7160, 0.1145),
}


@dataclass(frozen=True)
class ConfirmationFacts:
    finite_reactions: int
    rho: float
    positive_contexts: int
    context_count: int
    minimum_context_n: int


def _close(a: float, b: float, tol: float = 1e-12) -> bool:
    return math.isfinite(float(a)) and abs(float(a) - float(b)) <= tol


def validate_development_snapshot(snapshot: Mapping[str, object]) -> None:
    """Validate only already-frozen development statistics; do not reselect a metric."""
    # The adapter may normalize field names from the qualified snapshot before calling this.
    if int(snapshot["finite_primary_reactions"]) != EXPECTED_DEVELOPMENT_FINITE:
        raise ValueError("unexpected development finite-reaction count")
    if not _close(float(snapshot["eta2_usefulness_rho"]), EXPECTED_DEVELOPMENT_RHO):
        raise ValueError("unexpected development eta2/usefulness rho")
    if int(snapshot["positive_contexts"]) != EXPECTED_CONFIRMATION_CONTEXTS:
        raise ValueError("development context-direction count is not 16/16")


def validate_confirmation_primary(primary: Mapping[str, object]) -> ConfirmationFacts:
    if primary.get("status") != "PL2D_CONFIRMED":
        raise ValueError("PL2D primary status is not PL2D_CONFIRMED")
    if int(primary["finite_reactions"]) != EXPECTED_CONFIRMATION_FINITE:
        raise ValueError("unexpected confirmation finite-reaction count")
    if not _close(float(primary["overall_rho"]), EXPECTED_CONFIRMATION_RHO):
        raise ValueError("unexpected confirmation rho")
    if int(primary["context_count"]) != EXPECTED_CONFIRMATION_CONTEXTS:
        raise ValueError("unexpected confirmation context count")
    if int(primary["positive_contexts"]) != EXPECTED_CONFIRMATION_POSITIVE_CONTEXTS:
        raise ValueError("confirmation is not positive in 16/16 contexts")
    if int(primary["context_min_finite_reactions"]) != EXPECTED_CONFIRMATION_MIN_CONTEXT_N:
        raise ValueError("unexpected minimum confirmation context support")
    return ConfirmationFacts(
        finite_reactions=int(primary["finite_reactions"]),
        rho=float(primary["overall_rho"]),
        positive_contexts=int(primary["positive_contexts"]),
        context_count=int(primary["context_count"]),
        minimum_context_n=int(primary["context_min_finite_reactions"]),
    )


def validate_supportive(supportive: Mapping[str, object]) -> None:
    checks = (
        (supportive["bootstrap"]["rho_q025"], EXPECTED_BOOTSTRAP_Q025, "bootstrap q025"),
        (supportive["bootstrap"]["rho_median"], EXPECTED_BOOTSTRAP_MEDIAN, "bootstrap median"),
        (supportive["bootstrap"]["rho_q975"], EXPECTED_BOOTSTRAP_Q975, "bootstrap q975"),
        (supportive["eta2_vs_mean_information_advantage"]["rho"], EXPECTED_INFO_ADVANTAGE_RHO, "information advantage rho"),
        (supportive["partial_eta2_vs_usefulness_controlling_non_tie_coverage"]["rho"], EXPECTED_PARTIAL_RHO, "partial rho"),
        (supportive["directional_entropy3_vs_usefulness"]["rho"], EXPECTED_ENTROPY_RHO, "entropy rho"),
        (supportive["dominant_sign_mass_vs_usefulness"]["rho"], EXPECTED_DOMINANT_SIGN_RHO, "dominant-sign rho"),
        (supportive["n_supported_sign_states_vs_usefulness"]["rho"], EXPECTED_SIGN_STATES_RHO, "sign-state rho"),
        (supportive["pathway_sensitivity"]["exclude_proximal"]["rho"], EXPECTED_PATHWAY_EXCLUDE_PROXIMAL_RHO, "proximal sensitivity rho"),
        (supportive["pathway_sensitivity"]["exclude_same_subsystem_as_HEX1"]["rho"], EXPECTED_PATHWAY_EXCLUDE_SAME_SUBSYSTEM_RHO, "same-subsystem sensitivity rho"),
        (supportive["pathway_sensitivity"]["exclude_either"]["rho"], EXPECTED_PATHWAY_EXCLUDE_EITHER_RHO, "union sensitivity rho"),
    )
    for actual, expected, label in checks:
        if not _close(float(actual), expected):
            raise ValueError(f"unexpected {label}")
    if int(supportive["bootstrap"]["replicates"]) != 5000:
        raise ValueError("unexpected bootstrap replicate count")
    if int(supportive["bootstrap"]["finite_replicates"]) != 5000:
        raise ValueError("bootstrap contains non-finite replicates")


def validate_confirmation_qc(qc: Mapping[str, object]) -> None:
    counts = qc["counts"]
    expected = {
        "confirmation_cases": EXPECTED_CONFIRMATION_CASES,
        "non_tie_cases": EXPECTED_NON_TIE_CASES,
        "truth_tie_cases": EXPECTED_TRUTH_TIE_CASES,
        "confirmation_evaluation_reaction_rows": EXPECTED_EVALUATION_REACTION_ROWS,
    }
    for key, value in expected.items():
        if int(counts[key]) != value:
            raise ValueError(f"unexpected confirmation QC count: {key}")
    if not bool(qc["primary_result_frozen_before_supportive"]):
        raise ValueError("primary result was not frozen before supportive analysis")
    if not bool(qc["no_quantile_alias_weighting"]):
        raise ValueError("quantile aliases were used as weights")
    if float(qc["maximum_random_sign_identity_deviation"]) != 0.0:
        raise ValueError("random-sign identity deviation is nonzero")
    if float(qc["maximum_information_advantage_identity_deviation"]) != 0.0:
        raise ValueError("information-advantage identity deviation is nonzero")


def validate_pl2b_display_rows(rows: Sequence[Mapping[str, object]]) -> None:
    """Verify the three manuscript-facing PL2B rows without treating q=.5 as baseline."""
    indexed = {round(float(r["reliability_q"]), 12): r for r in rows}
    for q, expected in EXPECTED_PL2B_DISPLAY.items():
        if q not in indexed:
            raise ValueError(f"missing PL2B reliability row q={q}")
        row = indexed[q]
        got = tuple(round(float(row[k]), 4) for k in ("improved_fraction_of_non_ties", "tied_fraction_of_non_ties", "harmed_fraction_of_non_ties"))
        if got != expected:
            raise ValueError(f"PL2B manuscript rounding mismatch for q={q}: {got} != {expected}")
        if abs(sum(float(row[k]) for k in ("improved_fraction_of_non_ties", "tied_fraction_of_non_ties", "harmed_fraction_of_non_ties")) - 1.0) > 5e-12:
            raise ValueError(f"PL2B fractions do not sum to one for q={q}")


def claim_specs() -> list[dict[str, str]]:
    """Frozen manuscript-safe claim ledger skeleton."""
    return [
        {
            "claim_id": "C1_SIGN_ONLY_UTILITY_NOT_UNIVERSAL",
            "claim_level": "PRIMARY_CONTEXT",
            "statement": "At primary lambda=0.25, correct directional information improved only a subset of non-tie reaction/truth-pair cases; tied and harmed cases remain explicit, and q=0.5 is a random-sign control rather than the strong-anchor baseline.",
            "boundary": "Do not claim general recovery of reaction-B flux or that most reactions improve.",
        },
        {
            "claim_id": "C2_GEOMETRY_PREDICTS_DIRECTIONAL_USEFULNESS",
            "claim_level": "PRIMARY",
            "statement": "Strong-anchor sign-magnitude geometry predicts where correct weak directional information is useful, with a strong development association that was confirmed in an untouched reaction holdout.",
            "boundary": "Do not pool development and confirmation to create a new primary correlation and do not make a causal claim.",
        },
        {
            "claim_id": "C3_OCCURRENCE_STRONGER_THAN_GAIN_MAGNITUDE",
            "claim_level": "PRIMARY_INTERPRETATION",
            "statement": "Geometry predicts whether directional information is useful more strongly than the magnitude of its numerical information advantage.",
            "boundary": "Do not claim accurate prediction of gain magnitude.",
        },
        {
            "claim_id": "C4_CONTEXT_AND_PATHWAY_STABILITY",
            "claim_level": "SUPPORTIVE",
            "statement": "The confirmed eta2/usefulness association is positive across all 16 reconstruction-by-RNA contexts and remains similar under prespecified pathway exclusions.",
            "boundary": "Supportive analyses cannot be presented as additional confirmation gates.",
        },
        {
            "claim_id": "L1_DIRECTIONAL_TRUTH_COVERAGE",
            "claim_level": "LIMITATION",
            "statement": "A large fraction of benchmark cases are truth ties, so directional usefulness is defined only where a binary directional truth exists.",
            "boundary": "Do not imply that every reaction/context admits a meaningful binary direction.",
        },
    ]


def figure_panels() -> list[dict[str, str]]:
    return [
        {"panel": "A", "title": "BRIDGE concept", "source": "definitions", "message": "Strong-anchor geometry can make weak directional information useful only in selected geometries."},
        {"panel": "B", "title": "Sign-only utility benchmark", "source": "PL2B overall lambda=0.25", "message": "Show improved/tied/harmed fractions for q=0,0.5,1; label q=0.5 as random-sign control."},
        {"panel": "C", "title": "Development and confirmation", "source": "PL2D-H development snapshot + PL2D primary", "message": "Show development and untouched-holdout eta2/usefulness associations separately; never pool them."},
        {"panel": "D", "title": "Context confirmation", "source": "PL2D context results", "message": "Show all 16 context rhos and finite counts, including the weakest context."},
        {"panel": "E", "title": "Specificity and robustness", "source": "PL2D supportive results", "message": "Contrast usefulness with information-advantage association and show prespecified pathway sensitivities."},
    ]
