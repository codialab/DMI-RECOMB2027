#!/usr/bin/env python3
"""Export the frozen, matched Figure 2 D–G geometry/information points."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[4]
OUT = ROOT / "reproduced/derived/fig2/figure2_DG_points.csv.xz"
GAIN_TOL = 1e-12
IDENTITY_TOL = 5e-12
EXPECTED_UPSTREAM_MANIFESTS = {
    "outputs/dmi_bridge_pl1_predictability_landscape_v1_resumable_v1/BRIDGEPL1_MANIFEST.json": "b0f81620ff24f5791dd964a8063416784db86a358b0464b8c68a504e84e5f529",
    "outputs/dmi_bridge_pl2a_sign_only_prepare_v1/BRIDGEPL2A_MANIFEST.json": "8ae3eb60eb23614e379e38e3aa3cb72ec16c6bdb69aa9c74cc94499e31f1c312",
    "outputs/dmi_bridge_pl2b_sign_only_utility_v1/BRIDGEPL2B_MANIFEST.json": "7d9ca97a3f4ea5aa8b0b983e09ee83fc3c82f326d1b1cbf9c6137cfadf25539d",
    "outputs/dmi_bridge_pl2c_geometry_utility_development_v1/BRIDGEPL2C_MANIFEST.json": "6abcf149cd6adb294c0931585cec28f992b453451a5b5f03e10b3ffb2f9b9c39",
    "outputs/dmi_bridge_pl2d_confirmation_v1/BRIDGEPL2D_MANIFEST.json": "9f9eb730473e7ad209c1d5ec3340808ab614087b8ae1121d929eff61be12b699",
    "outputs/dmi_bridge_pl2d_hypothesis_freeze_v1/BRIDGEPL2DH_MANIFEST.json": "a1dd76d6b40e741e07efe64c34bc39eaad289a46bf6c870917a4f59cab688427",
    "outputs/dmi_bridge_pl2d_h_support_gate_amendment_v1/BRIDGEPL2DHA_MANIFEST.json": "fb7df9e208da2e6c9841d43a3b572038e59943a4221535705eef9f5dc0ed6a3a",
    "outputs/dmi_bridge_a21_dual_anchor_geometry_v1/BRIDGEA21_MANIFEST.json": "db9948d7b99447a695cf8aa21429dde6d642e1fd0632bcce52643deac59e6feb",
    "outputs/dmi_bridge_a22_dual_anchor_sign_only_utility_v1/BRIDGEA22_MANIFEST.json": "15b237151005538e5077a352868f873ef00466f4f5ea49307b781a8a8cd0e7cb",
    "outputs/dmi_bridge_a23_a1_a2_synthesis_v1/BRIDGEA23_MANIFEST.json": "b8e57074f3df23fe7708145bc01f51a19f946545b1c16e61628c2678dbcabd18",
}

CASE_USECOLS = [
    "truth_pair_id", "evaluation_id", "reaction_id", "truth_status",
    "lambda_0_25_correct_absolute_error_gain",
    "lambda_0_25_wrong_absolute_error_gain",
    "lambda_0_25_random_sign_expected_gain", "subsystem",
]
POINT_COLUMNS = [
    "schema_version", "panel", "panel_order", "anchor_setting", "reaction_split",
    "reaction_id", "reaction_name", "pathway",
    "x_direction_explained_magnitude_variance", "y_mean_g_info",
    "y_directional_usefulness_fraction", "n_eta2_evaluations",
    "n_g_info_evaluations", "n_directional_usefulness_evaluations",
    "n_non_tie_cases",
]
EXPECTED_USEFULNESS = {
    "D": (2454, 0.8571160663367106),
    "E": (2485, 0.8339440056039364),
    "F": (838, 0.8351435823183291),
    "G": (843, 0.8493657186460156),
}
# The frozen A23 correlation was computed from its persisted float arrays and
# registry-ordered pandas means; allow only last-bit/rank-tie summation drift.
USEFULNESS_RHO_ATOL = 1e-6


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def verify_manifest_artifacts(manifest_path: Path) -> dict[str, Any]:
    manifest = read_json(manifest_path)
    base = manifest_path.parent
    for relative, expected in manifest.get("artifact_sha256", {}).items():
        artifact = base / relative
        if not artifact.is_file() or sha256(artifact) != expected:
            raise RuntimeError(f"artifact hash mismatch or missing file: {artifact}")
    return manifest


def verify_sources() -> dict[str, Any]:
    manifests: dict[str, dict[str, Any]] = {}
    for relative, expected in EXPECTED_UPSTREAM_MANIFESTS.items():
        path = ROOT / relative
        actual = sha256(path)
        if expected and actual != expected:
            raise RuntimeError(f"frozen manifest hash mismatch: {relative}")
        manifests[relative] = verify_manifest_artifacts(path)

    wanted_status = {
        "outputs/dmi_bridge_pl1_predictability_landscape_v1_resumable_v1/BRIDGEPL1_MANIFEST.json": "PL1_PREDICTABILITY_LANDSCAPE_COMPLETE",
        "outputs/dmi_bridge_pl2a_sign_only_prepare_v1/BRIDGEPL2A_MANIFEST.json": "PL2A_READY_FOR_PL2B_SIGN_ONLY_UTILITY",
        "outputs/dmi_bridge_pl2b_sign_only_utility_v1/BRIDGEPL2B_MANIFEST.json": "PL2B_SIGN_ONLY_UTILITY_COMPLETE",
        "outputs/dmi_bridge_pl2c_geometry_utility_development_v1/BRIDGEPL2C_MANIFEST.json": "PL2C_DEVELOPMENT_GEOMETRY_UTILITY_COMPLETE",
        "outputs/dmi_bridge_pl2d_confirmation_v1/BRIDGEPL2D_MANIFEST.json": "PL2D_CONFIRMED",
        "outputs/dmi_bridge_a21_dual_anchor_geometry_v1/BRIDGEA21_MANIFEST.json": "BRIDGE_A21_DUAL_ANCHOR_GEOMETRY_FROZEN",
        "outputs/dmi_bridge_a22_dual_anchor_sign_only_utility_v1/BRIDGEA22_MANIFEST.json": "BRIDGE_A22_SIGN_ONLY_UTILITY_FROZEN",
        "outputs/dmi_bridge_a23_a1_a2_synthesis_v1/BRIDGEA23_MANIFEST.json": "BRIDGE_A23_A1_A2_SYNTHESIS_FROZEN",
    }
    for name, status in wanted_status.items():
        if manifests[name].get("status") != status:
            raise RuntimeError(f"unexpected terminal status in {name}")
    if manifests["outputs/dmi_bridge_a21_dual_anchor_geometry_v1/BRIDGEA21_MANIFEST.json"].get("a2_g_used") is not False:
        raise RuntimeError("A2.1 does not certify A2-G exclusion")
    if manifests["outputs/dmi_bridge_a22_dual_anchor_sign_only_utility_v1/BRIDGEA22_MANIFEST.json"].get("a2_g_used") is not False:
        raise RuntimeError("A2.2 does not certify A2-G exclusion")

    a23 = manifests["outputs/dmi_bridge_a23_a1_a2_synthesis_v1/BRIDGEA23_MANIFEST.json"]
    for relative, expected in a23["source_sha256"].items():
        if sha256(ROOT / relative) != expected:
            raise RuntimeError(f"A2.3 source-manifest binding mismatch: {relative}")
    return manifests


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def load_cases(
    case_dir: Path,
    prefix: str,
    pair_to_eval: dict[str, str],
    target_reactions: set[str],
) -> tuple[pd.DataFrame, dict[str, str]]:
    parts = sorted(case_dir.glob(f"{prefix}_CASE_OUTCOMES.part-*.tsv.xz"))
    if not parts:
        raise RuntimeError(f"no case parts found under {case_dir}")
    frames: list[pd.DataFrame] = []
    annotation: dict[str, str] = {}
    n_rows = 0
    for part in parts:
        for chunk in pd.read_csv(part, sep="\t", compression="xz", usecols=CASE_USECOLS, dtype=str, chunksize=100_000):
            chunk = chunk.loc[
                chunk["truth_pair_id"].isin(pair_to_eval)
                & chunk["reaction_id"].isin(target_reactions)
            ].copy()
            if chunk.empty:
                continue
            if not (chunk["evaluation_id"].to_numpy() == chunk["truth_pair_id"].map(pair_to_eval).to_numpy()).all():
                raise RuntimeError("case evaluation_id disagrees with the A2.3 matched truth-pair registry")
            n_rows += len(chunk)
            for rid, sub in chunk.groupby("reaction_id", sort=False):
                value = next((v for v in sub["subsystem"].tolist() if v), "")
                if value:
                    prior = annotation.setdefault(rid, value)
                    if prior != value:
                        raise RuntimeError(f"conflicting subsystem annotation for {rid}")
            frames.append(chunk)
    expected = len(pair_to_eval) * len(target_reactions)
    if n_rows != expected:
        raise RuntimeError(f"matched case population is {n_rows}, expected {expected} ({prefix})")
    frame = pd.concat(frames, ignore_index=True)
    if frame.duplicated(["truth_pair_id", "reaction_id"]).any():
        raise RuntimeError(f"duplicate matched truth-pair/reaction case keys ({prefix})")
    per_pair = frame.groupby("truth_pair_id", sort=False)["reaction_id"].nunique()
    if len(per_pair) != len(pair_to_eval) or not (per_pair == len(target_reactions)).all():
        raise RuntimeError(f"matched cases do not cover every reaction for every truth pair ({prefix})")

    ties = frame["truth_status"].eq("TRUTH_TIE")
    valid = frame["truth_status"].eq("NON_TIE")
    if not (ties | valid).all():
        raise RuntimeError("unexpected truth_status value")
    response = frame.loc[valid].copy()
    for col in CASE_USECOLS[-4:-1]:
        response[col] = pd.to_numeric(response[col], errors="coerce")
    gc = response["lambda_0_25_correct_absolute_error_gain"].to_numpy(dtype=float)
    gw = response["lambda_0_25_wrong_absolute_error_gain"].to_numpy(dtype=float)
    gr = response["lambda_0_25_random_sign_expected_gain"].to_numpy(dtype=float)
    if not (np.isfinite(gc) & np.isfinite(gw) & np.isfinite(gr)).all():
        raise RuntimeError("non-finite gain on a matched non-tie case")
    if (np.abs(gr - 0.5 * (gc + gw)) > IDENTITY_TOL).any():
        raise RuntimeError("random-sign gain identity failed")
    response["g_info"] = gc - gr
    if (np.abs(response["g_info"].to_numpy() - 0.5 * (gc - gw)) > IDENTITY_TOL).any():
        raise RuntimeError("information-advantage identity failed")
    response["directionally_useful"] = (gc > GAIN_TOL) & (response["g_info"].to_numpy() > GAIN_TOL)
    return response[["evaluation_id", "reaction_id", "truth_pair_id", "g_info", "directionally_useful"]], annotation


def main() -> None:
    manifests = verify_sources()

    a23_dir = ROOT / "outputs/dmi_bridge_a23_a1_a2_synthesis_v1"
    matched = read_tsv(a23_dir / "BRIDGEA23_MATCHED_POPULATION_AUDIT.tsv")
    matched = matched.loc[matched["matched"].str.lower().eq("true")]
    if len(matched) != 836 or matched["truth_pair_id"].nunique() != 836:
        raise RuntimeError("A2.3 matched truth-pair population is not exactly 836")
    pair_to_eval = dict(zip(matched["truth_pair_id"], matched["evaluation_id"]))
    matched_evaluations = set(pair_to_eval.values())

    a21_dir = ROOT / "outputs/dmi_bridge_a21_dual_anchor_geometry_v1"
    a21_features = sorted(a21_dir.glob("BRIDGEA21_REACTION_EVALUATION_FEATURES.part-*.tsv.xz"))
    if not a21_features:
        raise RuntimeError("A2.1 feature parts are missing")
    target_reactions: set[str] = set()
    for path in a21_features:
        target_reactions.update(pd.read_csv(path, sep="\t", usecols=["reaction_id"], dtype=str)["reaction_id"])
    if len(target_reactions) != 4179:
        raise RuntimeError(f"A2-L reaction universe is {len(target_reactions)}, expected 4179")
    a2_dir = ROOT / "outputs/dmi_bridge_a22_dual_anchor_sign_only_utility_v1"
    reaction_registry = read_tsv(a2_dir / "BRIDGEA22_REACTION_REGISTRY.tsv")
    reaction_order = {reaction_id: index for index, reaction_id in enumerate(reaction_registry["reaction_id"])}
    if set(reaction_order) != target_reactions:
        raise RuntimeError("matched feature reactions disagree with the frozen A2.2 reaction registry")

    split_path = ROOT / "outputs/dmi_bridge_pl2c_geometry_utility_development_v1/BRIDGEPL2C_REACTION_SPLIT.tsv"
    split_frame = read_tsv(split_path)
    split_map = dict(zip(split_frame["reaction_id"], split_frame["analysis_split"]))
    if not target_reactions.issubset(split_map):
        raise RuntimeError("matched target reactions do not all occur in the frozen PL2C split")
    split_by_reaction = {rid: split_map[rid] for rid in target_reactions}
    split_counts = pd.Series(split_by_reaction).value_counts().to_dict()
    if split_counts != {"DEVELOPMENT": 3134, "CONFIRMATION_HOLDOUT": 1045}:
        raise RuntimeError(f"unexpected matched split counts: {split_counts}")

    pl1_dir = ROOT / "outputs/dmi_bridge_pl1_predictability_landscape_v1_resumable_v1"
    pl1_registry = read_tsv(pl1_dir / "BRIDGEPL1_EVALUATION_REGISTRY.tsv")
    a21_registry = read_tsv(a21_dir / "BRIDGEA21_EVALUATION_REGISTRY.tsv")
    if set(pl1_registry["evaluation_id"]) != set(a21_registry["evaluation_id"]) or len(pl1_registry) != 400 or len(a21_registry) != 400:
        raise RuntimeError("A1 and A2-L evaluation identities differ")
    evaluation_order = {evaluation_id: index for index, evaluation_id in enumerate(pl1_registry["evaluation_id"])}

    # A2.1 contains all 16 evaluation strata in each of its feature parts.
    a2_frames = []
    for path in a21_features:
        frame = pd.read_csv(path, sep="\t", compression="xz", dtype={"evaluation_id": str, "reaction_id": str})
        frame = frame.loc[frame["reaction_id"].isin(target_reactions) & frame["evaluation_id"].isin(matched_evaluations), ["evaluation_id", "reaction_id", "sign_magnitude_eta2", "subsystem"]].copy()
        frame["sign_magnitude_eta2"] = pd.to_numeric(frame["sign_magnitude_eta2"], errors="coerce")
        frame.loc[~np.isfinite(frame["sign_magnitude_eta2"]), "sign_magnitude_eta2"] = np.nan
        a2_frames.append(frame)
    a2_feature_frame = pd.concat(a2_frames, ignore_index=True)
    a2_feature_frame["_evaluation_order"] = a2_feature_frame["evaluation_id"].map(evaluation_order)
    a2_feature_frame["_reaction_order"] = a2_feature_frame["reaction_id"].map(reaction_order)
    a2_feature_frame = a2_feature_frame.sort_values(["_evaluation_order", "_reaction_order"], kind="stable")
    if a2_feature_frame.duplicated(["evaluation_id", "reaction_id"]).any():
        raise RuntimeError("duplicate A2.1 evaluation/reaction feature keys")
    a2_x_series = a2_feature_frame.groupby("reaction_id", sort=True)["sign_magnitude_eta2"].mean()
    a2_nx = a2_feature_frame.groupby("reaction_id", sort=True)["sign_magnitude_eta2"].count().astype(int).to_dict()
    a2_x = a2_x_series.to_dict()
    a2_pathway = a2_feature_frame.groupby("reaction_id", sort=True)["subsystem"].agg(lambda values: next((v for v in values if v), "")).to_dict()

    # PL1 part files are partitioned by context, not by evaluation. Recompute the
    # A1 mean from concatenated rows so pandas grouped-mean order is preserved.
    a1_frames = []
    for part in range(4):
        path = pl1_dir / f"BRIDGEPL1_REACTION_EVALUATION_FEATURES.part-{part:03d}.tsv.xz"
        frame = pd.read_csv(path, sep="\t", compression="xz", dtype={"evaluation_id": str, "reaction_id": str})
        frame = frame.loc[frame["reaction_id"].isin(target_reactions) & frame["evaluation_id"].isin(matched_evaluations), ["evaluation_id", "reaction_id", "sign_magnitude_eta2", "subsystem"]].copy()
        frame["sign_magnitude_eta2"] = pd.to_numeric(frame["sign_magnitude_eta2"], errors="coerce")
        frame.loc[~np.isfinite(frame["sign_magnitude_eta2"]), "sign_magnitude_eta2"] = np.nan
        a1_frames.append(frame)
    a1_feature_frame = pd.concat(a1_frames, ignore_index=True)
    a1_feature_frame["_evaluation_order"] = a1_feature_frame["evaluation_id"].map(evaluation_order)
    a1_feature_frame["_reaction_order"] = a1_feature_frame["reaction_id"].map(reaction_order)
    a1_feature_frame = a1_feature_frame.sort_values(["_evaluation_order", "_reaction_order"], kind="stable")
    if a1_feature_frame.duplicated(["evaluation_id", "reaction_id"]).any():
        raise RuntimeError("duplicate A1 evaluation/reaction feature keys")
    a1_x = a1_feature_frame.groupby("reaction_id", sort=True)["sign_magnitude_eta2"].mean().to_dict()
    a1_nx = a1_feature_frame.groupby("reaction_id", sort=True)["sign_magnitude_eta2"].count().astype(int).to_dict()
    a1_pathway = a1_feature_frame.groupby("reaction_id", sort=True)["subsystem"].agg(lambda values: next((v for v in values if v), "")).to_dict()

    pair_dir = ROOT / "outputs/dmi_bridge_pl2b_sign_only_utility_v1"
    a1_y_rows, a1_case_pathways = load_cases(pair_dir, "BRIDGEPL2B", pair_to_eval, target_reactions)
    a2_y_rows, a2_case_pathways = load_cases(a2_dir, "BRIDGEA22", pair_to_eval, target_reactions)
    if a1_y_rows["truth_pair_id"].nunique() != 836 or a2_y_rows["truth_pair_id"].nunique() != 836:
        raise RuntimeError("matched case input did not retain exactly 836 truth pairs")

    def y_summary(rows: pd.DataFrame) -> tuple[dict[str, float], dict[str, int], dict[str, int], dict[str, float], dict[str, int]]:
        per_eval = rows.groupby(["evaluation_id", "reaction_id"], sort=True)["g_info"].mean().reset_index()
        per_eval["_evaluation_order"] = per_eval["evaluation_id"].map(evaluation_order)
        per_eval["_reaction_order"] = per_eval["reaction_id"].map(reaction_order)
        per_eval = per_eval.sort_values(["_evaluation_order", "_reaction_order"], kind="stable")
        per_reaction = per_eval.groupby("reaction_id", sort=True)["g_info"].mean()
        n_eval = per_eval.groupby("reaction_id", sort=True)["evaluation_id"].nunique().astype(int)
        n_cases = rows.groupby("reaction_id", sort=True)["truth_pair_id"].nunique().astype(int)
        usefulness_per_eval = rows.groupby(["evaluation_id", "reaction_id"], sort=True)["directionally_useful"].mean().reset_index()
        usefulness_per_eval["_evaluation_order"] = usefulness_per_eval["evaluation_id"].map(evaluation_order)
        usefulness_per_eval["_reaction_order"] = usefulness_per_eval["reaction_id"].map(reaction_order)
        usefulness_per_eval = usefulness_per_eval.sort_values(["_evaluation_order", "_reaction_order"], kind="stable")
        usefulness_per_reaction = usefulness_per_eval.groupby("reaction_id", sort=True)["directionally_useful"].mean()
        n_usefulness_eval = usefulness_per_eval.groupby("reaction_id", sort=True)["evaluation_id"].nunique().astype(int)
        return per_reaction.to_dict(), n_eval.to_dict(), n_cases.to_dict(), usefulness_per_reaction.to_dict(), n_usefulness_eval.to_dict()

    a1_y, a1_ny, a1_nc, a1_u, a1_nu = y_summary(a1_y_rows)
    a2_y, a2_ny, a2_nc, a2_u, a2_nu = y_summary(a2_y_rows)
    pathways = {rid: a1_pathway.get(rid) or a2_pathway.get(rid) or a1_case_pathways.get(rid) or a2_case_pathways.get(rid) or "" for rid in target_reactions}

    panel_defs = [
        ("D", 1, "A1", "DEVELOPMENT", a1_x, a1_nx, a1_y, a1_ny, a1_nc, a1_u, a1_nu),
        ("E", 2, "A2-L", "DEVELOPMENT", a2_x, a2_nx, a2_y, a2_ny, a2_nc, a2_u, a2_nu),
        ("F", 3, "A1", "CONFIRMATION_HOLDOUT", a1_x, a1_nx, a1_y, a1_ny, a1_nc, a1_u, a1_nu),
        ("G", 4, "A2-L", "CONFIRMATION_HOLDOUT", a2_x, a2_nx, a2_y, a2_ny, a2_nc, a2_u, a2_nu),
    ]
    output_rows = []
    for panel, order, anchor, split, xs, nxs, ys, nys, ncs, us, nus in panel_defs:
        for rid in sorted(target_reactions):
            if split_by_reaction[rid] != split:
                continue
            x, y = xs.get(rid, math.nan), ys.get(rid, math.nan)
            if not (math.isfinite(x) and math.isfinite(y)):
                continue
            output_rows.append({
                "schema_version": "recomb.fig2.dg.points.v1",
                "panel": panel,
                "panel_order": order,
                "anchor_setting": anchor,
                "reaction_split": "development" if split == "DEVELOPMENT" else "held_out",
                "reaction_id": rid,
                "reaction_name": "",
                "pathway": pathways.get(rid, ""),
                "x_direction_explained_magnitude_variance": x,
                "y_mean_g_info": y,
                "y_directional_usefulness_fraction": us.get(rid, math.nan),
                "n_eta2_evaluations": int(nxs.get(rid, 0)),
                "n_g_info_evaluations": int(nys.get(rid, 0)),
                "n_directional_usefulness_evaluations": int(nus.get(rid, 0)),
                "n_non_tie_cases": int(ncs.get(rid, 0)),
            })
    result = pd.DataFrame(output_rows, columns=POINT_COLUMNS)
    if result.duplicated(["panel", "reaction_id"]).any() or not np.isfinite(result[["x_direction_explained_magnitude_variance", "y_mean_g_info", "y_directional_usefulness_fraction"]].to_numpy(dtype=float)).all():
        raise RuntimeError("final export has duplicate keys or non-finite plotted coordinates")
    if set(result["panel"]) != {"D", "E", "F", "G"}:
        raise RuntimeError("one or more panels have no finite paired points")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(OUT, index=False, encoding="utf-8", lineterminator="\n", compression={"method": "xz", "preset": 9}, float_format="%.17g")
    for panel in "DEFG":
        panel_data = result.loc[result["panel"].eq(panel)]
        finite = np.isfinite(panel_data["x_direction_explained_magnitude_variance"]) & np.isfinite(panel_data["y_directional_usefulness_fraction"])
        plotted = panel_data.loc[finite]
        rho = spearmanr(plotted["x_direction_explained_magnitude_variance"], plotted["y_directional_usefulness_fraction"]).statistic
        expected_n, expected_rho = EXPECTED_USEFULNESS[panel]
        if len(plotted) != expected_n or not np.isclose(rho, expected_rho, rtol=0, atol=USEFULNESS_RHO_ATOL):
            raise RuntimeError(f"directional-usefulness parity mismatch for {panel}: n={len(plotted)} rho={rho:.16g}; expected n={expected_n} rho={expected_rho:.16g}")
        print(f"panel={panel} usefulness_n={len(plotted)} usefulness_spearman_rho={rho:.16g}")
    print(f"wrote={OUT.relative_to(ROOT)} rows={len(result)} sha256={sha256(OUT)}")


if __name__ == "__main__":
    main()
