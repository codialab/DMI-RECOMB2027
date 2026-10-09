#!/usr/bin/env python3
"""Post-hoc reviewer sensitivity analyses using released compact tables only.

This script does not invoke production pipelines or change scientific definitions.
Run from the repository root with:

    python analyses/reviewer_checks/run_reviewer_checks.py

Tables and an input-hash manifest are written under reproduced/reviewer_checks/.
"""
from __future__ import annotations

import hashlib
import io
import json
import lzma
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "reproduced/reviewer_checks"
FIG2 = ROOT / "data/figure_inputs/fig2/figure2_DG_points.csv.xz"
FIG3 = ROOT / "data/figure_inputs/fig3/fig3_gain_reaction_evaluation.parquet.xz"
FIG4_DIR = ROOT / "data/figure_inputs/fig4/fig4_geometry_paired"
SEED = 20271027
BOOTSTRAP_REPLICATES = 2000
PATHWAY_BOOTSTRAP_REPLICATES = 2000
FEATURES = ("eta2", "H_dir", "D_dir", "K_dir")
MODEL_FEATURES = {
    "M1": ("H_dir",),
    "M2": ("eta2",),
    "M3": ("H_dir", "eta2"),
    "M4": ("H_dir", "D_dir"),
    "M5": ("H_dir", "D_dir", "eta2"),
}
PANELS = {
    "D": ("A1", "development"),
    "E": ("A2", "development"),
    "F": ("A1", "held_out"),
    "G": ("A2", "held_out"),
}
GEOMETRY_COLUMNS = {
    "eta2": "direction_explained_magnitude_variance",
    "H_dir": "H_dir",
    "D_dir": "dominant_direction_mass",
    "K_dir": "supported_sign_state_count",
    "non_tie_coverage": "non_tie_coverage",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def finite_pair(frame: pd.DataFrame, x: str, y: str) -> pd.DataFrame:
    mask = np.isfinite(frame[x].to_numpy(dtype=float)) & np.isfinite(frame[y].to_numpy(dtype=float))
    return frame.loc[mask].copy()


def holdout_product_weights(left_weights, right_weights, left_index: int, right_index: int):
    """Exclude one truth vector per condition, then independently renormalize."""
    left=np.asarray(left_weights,dtype=float).copy(); right=np.asarray(right_weights,dtype=float).copy()
    if left.ndim!=1 or right.ndim!=1 or len(left)<2 or len(right)<2:
        raise ValueError("each condition must have at least two candidate weights")
    if not (0<=left_index<len(left) and 0<=right_index<len(right)):
        raise IndexError("held-out truth index outside candidate support")
    if not np.isfinite(left).all() or not np.isfinite(right).all() or (left<0).any() or (right<0).any():
        raise ValueError("candidate weights must be finite and non-negative")
    left[left_index]=0.; right[right_index]=0.
    if left.sum()<=0 or right.sum()<=0:
        raise ValueError("truth holdout leaves no positive support")
    return left/left.sum(),right/right.sum()


def spearman(x: Iterable[float], y: Iterable[float]) -> float:
    x, y = np.asarray(list(x), dtype=float), np.asarray(list(y), dtype=float)
    if len(x) < 2 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return float("nan")
    return float(spearmanr(x, y).statistic)


def partial_rank(x: Iterable[float], y: Iterable[float], controls: Iterable[Iterable[float]]) -> float:
    """Pearson association between rank residuals after linear control adjustment."""
    x, y = np.asarray(list(x), float), np.asarray(list(y), float)
    z = np.column_stack([np.asarray(list(v), float) for v in controls])
    xr, yr = rankdata(x), rankdata(y)
    zr = np.column_stack([rankdata(z[:, i]) for i in range(z.shape[1])])
    design = np.column_stack([np.ones(len(x)), zr])
    rx = xr - design @ np.linalg.lstsq(design, xr, rcond=None)[0]
    ry = yr - design @ np.linalg.lstsq(design, yr, rcond=None)[0]
    if np.ptp(rx) == 0 or np.ptp(ry) == 0:
        return float("nan")
    return float(np.corrcoef(rx, ry)[0, 1])


def bootstrap_corr_delta(
    frame: pd.DataFrame, left: str, right: str, y: str, *, seed: int, replicates: int
) -> tuple[float, float, int]:
    """Paired reaction bootstrap interval for rho(left,y)-rho(right,y)."""
    d = frame[[left, right, y]].dropna().to_numpy(float)
    d = d[np.isfinite(d).all(axis=1)]
    rng = np.random.default_rng(seed)
    vals = np.full(replicates, np.nan)
    for i in range(replicates):
        sample = d[rng.integers(0, len(d), len(d))]
        vals[i] = spearman(sample[:, 0], sample[:, 2]) - spearman(sample[:, 1], sample[:, 2])
    vals = vals[np.isfinite(vals)]
    if not len(vals):
        return float("nan"), float("nan"), 0
    return tuple(map(float, np.quantile(vals, [0.025, 0.975]))) + (len(vals),)


def bootstrap_rho(frame: pd.DataFrame, x: str, y: str, *, seed: int, replicates: int) -> tuple[float,float,int]:
    d=frame[[x,y]].dropna().to_numpy(float)
    d=d[np.isfinite(d).all(axis=1)]
    rng=np.random.default_rng(seed)
    vals=np.full(replicates,np.nan)
    for i in range(replicates):
        sample=d[rng.integers(0,len(d),len(d))]
        vals[i]=spearman(sample[:,0],sample[:,1])
    vals=vals[np.isfinite(vals)]
    if not len(vals): return float("nan"),float("nan"),0
    lo,hi=np.quantile(vals,[.025,.975])
    return float(lo),float(hi),len(vals)


def bootstrap_strength_delta(frame: pd.DataFrame, first: str, second: str, y: str, *, seed: int, replicates: int) -> tuple[float,float,int,float]:
    """Paired bootstrap for |rho(first,y)| - |rho(second,y)|."""
    d=frame[[first,second,y]].dropna().to_numpy(float)
    d=d[np.isfinite(d).all(axis=1)]
    rng=np.random.default_rng(seed)
    values=np.full(replicates,np.nan)
    point=abs(spearman(d[:,0],d[:,2]))-abs(spearman(d[:,1],d[:,2]))
    for i in range(replicates):
        sample=d[rng.integers(0,len(d),len(d))]
        values[i]=abs(spearman(sample[:,0],sample[:,2]))-abs(spearman(sample[:,1],sample[:,2]))
    values=values[np.isfinite(values)]
    if not len(values): return float("nan"),float("nan"),0,point
    lo,hi=np.quantile(values,[.025,.975])
    return float(lo),float(hi),len(values),point


def paired_anchor_rho_bootstrap(frame: pd.DataFrame, descriptor: str, *, seed: int, replicates: int) -> tuple[float, float, int, float]:
    """Paired reaction bootstrap of rho(A1 descriptor, A1 outcome)-rho(A2...)."""
    cols = [f"{descriptor}_A1", "usefulness_A1", f"{descriptor}_A2", "usefulness_A2"]
    d = frame[cols].dropna().to_numpy(float)
    d = d[np.isfinite(d).all(axis=1)]
    rng = np.random.default_rng(seed)
    values = np.full(replicates, np.nan)
    point = spearman(d[:, 0], d[:, 1]) - spearman(d[:, 2], d[:, 3])
    for i in range(replicates):
        sample = d[rng.integers(0, len(d), len(d))]
        values[i] = spearman(sample[:, 0], sample[:, 1]) - spearman(sample[:, 2], sample[:, 3])
    values = values[np.isfinite(values)]
    if not len(values):
        return float("nan"), float("nan"), 0, point
    return tuple(map(float, np.quantile(values, [0.025, 0.975]))) + (len(values), point)


def pathway_block_bootstraps(frame: pd.DataFrame, x: str, y: str, *, seed: int) -> tuple[float, float, int, float, float, int]:
    """Intervals for reaction-weighted rho and equal-pathway rho, resampling blocks."""
    usable = frame.loc[frame["pathway"].notna() & frame["pathway"].ne("")].copy()
    groups = usable["pathway"].drop_duplicates().to_numpy()
    rng = np.random.default_rng(seed)
    raw_blocks = {g: usable.loc[usable["pathway"].eq(g), [x, y]].to_numpy(float) for g in groups}
    summary = usable.groupby("pathway", sort=True).agg(x=(x,"mean"), y=(y,"mean"))
    summary_blocks = summary[["x","y"]].to_numpy(float)
    reaction_values = np.full(PATHWAY_BOOTSTRAP_REPLICATES, np.nan)
    summary_values = np.full(PATHWAY_BOOTSTRAP_REPLICATES, np.nan)
    for i in range(PATHWAY_BOOTSTRAP_REPLICATES):
        selected = rng.integers(0, len(groups), len(groups))
        raw = np.concatenate([raw_blocks[groups[j]] for j in selected], axis=0)
        reaction_values[i] = spearman(raw[:,0], raw[:,1])
        ps = summary_blocks[selected]
        summary_values[i] = spearman(ps[:,0], ps[:,1])
    reaction_values = reaction_values[np.isfinite(reaction_values)]
    summary_values = summary_values[np.isfinite(summary_values)]
    rlo, rhi = np.quantile(reaction_values,[.025,.975]) if len(reaction_values) else (np.nan,np.nan)
    slo, shi = np.quantile(summary_values,[.025,.975]) if len(summary_values) else (np.nan,np.nan)
    return float(rlo),float(rhi),len(reaction_values),float(slo),float(shi),len(summary_values)


def load_geometry(eval_ids: dict[str, set[str]]) -> pd.DataFrame:
    """Read only columns needed from compressed partitions; never persist expanded data."""
    pieces = []
    for path in sorted(FIG4_DIR.glob("algorithm=*/part.parquet.xz")):
        raw = lzma.decompress(path.read_bytes())
        # All rows are used for their respective matched evaluation IDs only.
        cols = ["evaluation_id", "reaction_id", "subsystem"]
        for prefix in ("A1", "A2"):
            cols += [f"{prefix}_{v}" for v in GEOMETRY_COLUMNS.values()]
        df = pd.read_parquet(io.BytesIO(raw), columns=cols)
        for anchor in ("A1", "A2"):
            wanted = eval_ids[anchor]
            sub = df.loc[df["evaluation_id"].astype(str).isin(wanted), ["evaluation_id", "reaction_id", "subsystem", *[f"{anchor}_{v}" for v in GEOMETRY_COLUMNS.values()]]].copy()
            sub = sub.rename(columns={f"{anchor}_{v}": k for k, v in GEOMETRY_COLUMNS.items()})
            sub["anchor"] = anchor
            pieces.append(sub)
    if not pieces:
        raise RuntimeError("no geometry observations found")
    result = pd.concat(pieces, ignore_index=True)
    if result.duplicated(["anchor", "evaluation_id", "reaction_id"]).any():
        raise RuntimeError("duplicate anchor/evaluation/reaction geometry key")
    for col in (*FEATURES, "non_tie_coverage"):
        result[col] = pd.to_numeric(result[col], errors="coerce")
        result.loc[~np.isfinite(result[col]), col] = np.nan
    return result


def aggregate_geometry(geometry: pd.DataFrame) -> pd.DataFrame:
    """Original matched-evaluation aggregation: equal mean over finite evaluations."""
    ordered = geometry.sort_values(["anchor", "evaluation_id", "reaction_id"], kind="stable")
    agg = ordered.groupby(["anchor", "reaction_id"], sort=True).agg(
        **{name: (name, "mean") for name in (*FEATURES, "non_tie_coverage")},
        n_descriptor_evaluations=("eta2", "count"),
        pathway=("subsystem", lambda s: next((v for v in s if isinstance(v, str) and v), "")),
    ).reset_index()
    return agg


def cv_predictions(frame: pd.DataFrame, feature_names: tuple[str, ...], folds: list, *, fit_frame: pd.DataFrame | None = None) -> np.ndarray:
    """OOF predictions for development; optionally fit on one frame and predict another."""
    y = frame["usefulness"].to_numpy(float)
    if not feature_names:
        pred = np.empty(len(frame), dtype=float)
        if fit_frame is None:
            for train, test in folds:
                pred[test] = y[train].mean()
        else:
            pred.fill(float(fit_frame["usefulness"].mean()))
        return pred
    x = frame[list(feature_names)].to_numpy(float)
    pred = np.empty(len(frame), dtype=float)
    if fit_frame is None:
        for train, test in folds:
            model = make_pipeline(StandardScaler(), Ridge(alpha=1.0))
            model.fit(x[train], y[train])
            pred[test] = model.predict(x[test])
    else:
        model = make_pipeline(StandardScaler(), Ridge(alpha=1.0))
        model.fit(fit_frame[list(feature_names)].to_numpy(float), fit_frame["usefulness"].to_numpy(float))
        pred[:] = model.predict(x)
    return pred


def make_group_folds(frame: pd.DataFrame, *, n_splits: int = 5) -> list[tuple[np.ndarray,np.ndarray]]:
    """Create GroupKFold indices and assert complete coverage and zero leakage."""
    groups=frame["pathway"].fillna("").astype(str).to_numpy()
    if len(set(groups)) < n_splits:
        raise RuntimeError(f"need at least {n_splits} groups; found {len(set(groups))}")
    folds=list(GroupKFold(n_splits=n_splits).split(np.zeros(len(frame)),frame["usefulness"],groups))
    test_counts=np.zeros(len(frame),dtype=int)
    for train,test in folds:
        if set(groups[train]) & set(groups[test]):
            raise RuntimeError("pathway leakage detected between train and test fold")
        test_counts[test]+=1
    if not np.all(test_counts==1):
        raise RuntimeError("GroupKFold test folds do not cover every reaction exactly once")
    return folds


def evaluate_models(dev: pd.DataFrame, confirm: pd.DataFrame | None, panel: str, population_scope: str, folds: list[tuple[np.ndarray,np.ndarray]], confirmation_panel: str | None = None) -> tuple[list[dict], list[pd.DataFrame]]:
    specs = {
        "M0": ((),"linear_entropy"),
        "M1": (("H_dir",),"linear_entropy"),
        "M2": (("eta2",),"linear_entropy"),
        "M3": (("H_dir","eta2"),"linear_entropy"),
        "M4": (("H_dir","D_dir"),"linear_entropy"),
        "M5": (("H_dir","D_dir","eta2"),"linear_entropy"),
        "N1_entropy_cubic": (("H_dir","H_dir_sq","H_dir_cu"),"nonlinear_entropy"),
        "N2_entropy_cubic_plus_eta2": (("H_dir","H_dir_sq","H_dir_cu","eta2"),"nonlinear_entropy"),
        "C1_entropy_coverage": (("H_dir","non_tie_coverage"),"coverage_adjusted"),
        "C2_entropy_coverage_plus_eta2": (("H_dir","non_tie_coverage","eta2"),"coverage_adjusted"),
    }
    for frame in (dev,confirm):
        if frame is not None:
            frame["H_dir_sq"]=frame.H_dir**2
            frame["H_dir_cu"]=frame.H_dir**3
    out, prediction_frames = [], []
    baseline_names={"linear_entropy":"M1","nonlinear_entropy":"N1_entropy_cubic","coverage_adjusted":"C1_entropy_coverage"}
    eval_frames=[("development_pathway_grouped_oof",panel,dev,None)]
    if confirm is not None and len(confirm):
        eval_frames.append(("confirmation_posthoc_fixed_development_fit",confirmation_panel,confirm,dev))
    for evaluation,evaluation_panel,eval_frame,fit_frame in eval_frames:
        eval_preds=[]
        for name,(predictors,comparison) in specs.items():
            pred=cv_predictions(eval_frame,predictors,folds if fit_frame is None else [],fit_frame=fit_frame)
            eval_preds.append((name,comparison,pred))
        for name,comparison,pred in eval_preds:
            base_name=baseline_names[comparison]
            base_pred=next(p for n,c,p in eval_preds if n==base_name)
            mae=float(np.mean(np.abs(eval_frame.usefulness.to_numpy()-pred)))
            base_mae=float(np.mean(np.abs(eval_frame.usefulness.to_numpy()-base_pred)))
            rho=spearman(pred,eval_frame.usefulness)
            base_rho=spearman(base_pred,eval_frame.usefulness)
            out.append({"panel":panel,"fit_panel":panel,"evaluation_panel":evaluation_panel,"population_scope":population_scope,"evaluation":evaluation,"model":name,
                        "comparison":comparison,"n":len(eval_frame),"mae":mae,"spearman_prediction_usefulness":rho,
                        "baseline_model":base_name,"mae_improvement_vs_baseline":base_mae-mae,
                        "spearman_change_vs_baseline":rho-base_rho})
            if evaluation=="development_pathway_grouped_oof":
                prediction_frames.append(pd.DataFrame({"panel":panel,"population_scope":population_scope,"reaction_id":dev.reaction_id,"model":name,"observed":dev.usefulness,"prediction":pred}))
    return out, prediction_frames


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    # Match the original plotter and `reproduce.py validate` numeric parser;
    # default float parsing perturbs exact ties in the reaction-level response.
    points = pd.read_csv(FIG2, compression="xz", float_precision="round_trip")
    gain_raw = pd.read_parquet(io.BytesIO(lzma.decompress(FIG3.read_bytes())), columns=["anchor_setting", "evaluation_id"])
    # Only the 370 evaluations represented in the original eligible utility table
    # enter the original matched-pair descriptor aggregation.
    # The frozen schema uses the historical label A2-L; normalize only the
    # internal join key and preserve the source label in the manifest/report.
    eval_ids = {
        "A1": set(gain_raw.loc[gain_raw.anchor_setting.eq("A1"), "evaluation_id"].astype(str)),
        "A2": set(gain_raw.loc[gain_raw.anchor_setting.eq("A2-L"), "evaluation_id"].astype(str)),
    }
    if len(eval_ids["A1"]) != 370 or len(eval_ids["A2"]) != 370:
        raise RuntimeError(f"unexpected eligible evaluation counts: { {k: len(v) for k,v in eval_ids.items()} }")
    geometry = load_geometry(eval_ids)
    aggregated = aggregate_geometry(geometry)

    panel_frames: dict[str, pd.DataFrame] = {}
    parity = []
    published_statistics = pd.read_csv(ROOT / "tables/manuscript/fig2_panel_statistics.tsv",sep="\t").set_index("panel")
    association_rows = []
    for panel, (anchor, split) in PANELS.items():
        original = points.loc[points.panel.eq(panel)].copy()
        expected_split = "development" if panel in ("D","E") else "held_out"
        expected_anchor = "A1" if panel in ("D","F") else "A2-L"
        if set(original.reaction_split.astype(str)) != {expected_split} or set(original.anchor_setting.astype(str)) != {expected_anchor}:
            raise RuntimeError(f"{panel}: frozen anchor/split assignment mismatch")
        if original.reaction_id.duplicated().any():
            raise RuntimeError(f"{panel}: duplicate frozen reaction identity")
        merged = original.merge(aggregated.loc[aggregated.anchor.eq(anchor)], on="reaction_id", how="left", validate="one_to_one")
        # The Figure 2 pathway annotation is the frozen reaction-level label.
        # Retain it as the grouping authority while checking the geometry label.
        if "pathway_x" in merged:
            left = merged.pathway_x.fillna("").astype(str)
            right = merged.pathway_y.fillna("").astype(str)
            conflicts = left.ne("") & right.ne("") & left.ne(right)
            if conflicts.any():
                raise RuntimeError(f"{panel}: pathway labels disagree for {int(conflicts.sum())} reactions")
            merged["pathway"] = left.where(left.ne(""), right)
        if merged["eta2"].notna().sum() != len(original):
            raise RuntimeError(f"{panel}: geometry join lost original eligible Figure 2 reactions")
        if merged["pathway"].isna().any() or merged.pathway.astype(str).eq("").any():
            raise RuntimeError(f"{panel}: incomplete reaction-to-pathway join")
        if not np.array_equal(merged.n_descriptor_evaluations.to_numpy(int),original.n_eta2_evaluations.to_numpy(int)):
            raise RuntimeError(f"{panel}: eta2 evaluation aggregation denominator mismatch")
        delta = np.abs(merged["eta2"].to_numpy(float) - merged["x_direction_explained_magnitude_variance"].to_numpy(float))
        rho_from_points=spearman(original.x_direction_explained_magnitude_variance, original.y_directional_usefulness_fraction)
        published_rho=float(published_statistics.loc[panel,"spearman_rho"])
        parity.append({"panel": panel, "anchor": anchor, "split": split, "frozen_fig2_n": len(original),
                       "eta2_join_n": int(np.isfinite(merged.eta2).sum()), "max_abs_eta2_difference": float(np.max(delta)),
                       "frozen_spearman": rho_from_points,"published_table_spearman":published_rho,"rho_delta_frozen_minus_published":rho_from_points-published_rho,
                       "eta2_recomputed_spearman": spearman(merged.eta2, merged.y_directional_usefulness_fraction),
                       "population_and_eta2_parity_pass": bool(np.max(delta) <= 1e-10 and len(original) == {"D":2454,"E":2485,"F":838,"G":843}[panel]),
                       "published_rho_parity_pass": bool(abs(rho_from_points-published_rho)<=1e-10),
                       "eta2_evaluation_counts_match":True,"reaction_pathway_join_complete":True})
        merged = merged.rename(columns={"y_directional_usefulness_fraction": "usefulness"})
        # Keep each panel's original eligibility; do not intersect anchor or split populations here.
        panel_frames[panel] = merged

    # All frozen population, association, aggregation, and identity checks must
    # pass before any new descriptor, bootstrap, pathway, or predictive analysis.
    expected_n={"D":2454,"E":2485,"F":838,"G":843}
    if set(points.panel.unique()) != set(PANELS):
        raise RuntimeError("Figure 2 panel set differs from frozen D/E/F/G populations")
    if any(len(panel_frames[p]) != expected_n[p] for p in PANELS):
        raise RuntimeError("frozen Figure 2 panel population count mismatch")
    if set(panel_frames["D"].reaction_id)&set(panel_frames["F"].reaction_id):
        raise RuntimeError("A1 development and original confirmation assignments overlap")
    if set(panel_frames["E"].reaction_id)&set(panel_frames["G"].reaction_id):
        raise RuntimeError("A2 development and confirmation-split assignments overlap")
    if not all(row["population_and_eta2_parity_pass"] and row["published_rho_parity_pass"] and row["eta2_evaluation_counts_match"] and row["reaction_pathway_join_complete"] for row in parity):
        raise RuntimeError("frozen Figure 2, descriptor, aggregation, or identity parity failed")

    required=[*FEATURES,"non_tie_coverage","usefulness","pathway"]
    grouped_folds={}; fold_records=[]
    for dev_panel in ("D","E"):
        for scope in ("FULL","MAPPED_ONLY"):
            source=panel_frames[dev_panel]
            if scope=="MAPPED_ONLY": source=source.loc[source.pathway.ne("UNMAPPED_PATHWAY")]
            development=source.dropna(subset=required).copy()
            folds=make_group_folds(development,n_splits=5)
            grouped_folds[(dev_panel,scope)]=folds
            for fold_index,(_,test) in enumerate(folds):
                fold_records.extend({"panel":dev_panel,"population_scope":scope,"reaction_id":development.iloc[i].reaction_id,"pathway":development.iloc[i].pathway,"fold":fold_index} for i in test)

    # Matched finite populations for descriptor comparisons. FULL preserves every
    # original Figure 2 row; MAPPED_ONLY excludes only UNMAPPED_PATHWAY.
    association_rows=[]
    for panel,(anchor,split) in PANELS.items():
        original=panel_frames[panel]
        populations={"FULL":original,"MAPPED_ONLY":original.loc[original.pathway.ne("UNMAPPED_PATHWAY")].copy()}
        for scope,pop in populations.items():
            common=pop.dropna(subset=[*FEATURES,"non_tie_coverage","usefulness"]).copy()
            for descriptor in FEATURES:
                lo,hi,nboot=bootstrap_rho(common,descriptor,"usefulness",seed=SEED+list(PANELS).index(panel)*20+FEATURES.index(descriptor)+(0 if scope=="FULL" else 100),replicates=BOOTSTRAP_REPLICATES)
                assoc=spearman(common[descriptor],common.usefulness)
                sign_aligned=common[descriptor] * (-1 if descriptor=="D_dir" else 1)
                aligned_frame=common.assign(_sign_aligned=sign_aligned)
                alo,ahi,an=bootstrap_rho(aligned_frame,"_sign_aligned","usefulness",seed=SEED+list(PANELS).index(panel)*20+FEATURES.index(descriptor)+(200 if scope=="MAPPED_ONLY" else 0),replicates=BOOTSTRAP_REPLICATES)
                association_rows.append({"panel":panel,"anchor":anchor,"split":split,"population_scope":scope,"descriptor":descriptor,
                    "n_original_panel":len(pop),"n_common_descriptors_response":len(common),"missing_from_common":len(pop)-len(common),
                    "spearman":assoc,"spearman_strength_abs":abs(assoc),"sign_aligned_spearman":(-assoc if descriptor=="D_dir" else assoc),
                    "signed_rho_bootstrap_q025":lo,"signed_rho_bootstrap_q975":hi,"sign_aligned_rho_bootstrap_q025":alo,"sign_aligned_rho_bootstrap_q975":ahi,
                    "finite_bootstrap_replicates":min(nboot,an)})
            association_rows.append({"panel":panel,"anchor":anchor,"split":split,"population_scope":scope,"descriptor":"eta2_partial_H_dir",
                "n_original_panel":len(pop),"n_common_descriptors_response":len(common),"missing_from_common":len(pop)-len(common),
                "spearman":partial_rank(common.eta2,common.usefulness,[common.H_dir]) if len(common) else np.nan})
            association_rows.append({"panel":panel,"anchor":anchor,"split":split,"population_scope":scope,"descriptor":"eta2_partial_H_dir_non_tie_coverage",
                "n_original_panel":len(pop),"n_common_descriptors_response":len(common),"missing_from_common":len(pop)-len(common),
                "spearman":partial_rank(common.eta2,common.usefulness,[common.H_dir,common.non_tie_coverage]) if len(common) else np.nan})

    # Paired descriptor differences are restricted to reaction IDs available in both
    # anchor panels for the same frozen split; these responses derive from the 836
    # matched truth-pair identities in the released Figure 2 construction.
    paired_corr = []
    for p1, p2, label in (("D", "E", "development"), ("F", "G", "confirmation")):
        a, b = panel_frames[p1], panel_frames[p2]
        joined = a.merge(b, on="reaction_id", suffixes=("_A1", "_A2"), validate="one_to_one")
        for variable in (*FEATURES, "usefulness"):
            paired_corr.append({"split": label, "variable": variable, "n_matched_reactions": int(joined[[f"{variable}_A1", f"{variable}_A2"]].dropna().shape[0]),
                                "mean_A2_minus_A1": float((joined[f"{variable}_A2"] - joined[f"{variable}_A1"]).mean())})
        common = joined.dropna(subset=[*[f"{f}_A1" for f in FEATURES], *[f"{f}_A2" for f in FEATURES], "usefulness_A1", "usefulness_A2"])
        for descriptor in FEATURES:
            lo, hi, nboot, rho_delta = paired_anchor_rho_bootstrap(common,descriptor,seed=SEED + (0 if label == "development" else 1) + FEATURES.index(descriptor),replicates=BOOTSTRAP_REPLICATES)
            paired_corr.append({"split": label, "variable": f"rho_A1_minus_A2_{descriptor}", "n_matched_reactions": len(common), "correlation_difference_A1_minus_A2": rho_delta, "paired_bootstrap_q025":lo,"paired_bootstrap_q975":hi,"finite_bootstrap_replicates":nboot,"truth_pairs_per_anchor":836})

    # Development models are specified only on A1 development. A1 confirmation
    # receives no tuning or refitting; A2 development and holdout are robustness
    # prediction exercises fit on their own development panel, marked post-freeze.
    model_rows, prediction_frames = [], []
    for dev_panel, confirm_panel in (("D", "F"), ("E", "G")):
        for scope in ("FULL","MAPPED_ONLY"):
            dev_source=panel_frames[dev_panel]
            conf_source=panel_frames[confirm_panel]
            if scope=="MAPPED_ONLY":
                dev_source=dev_source.loc[dev_source.pathway.ne("UNMAPPED_PATHWAY")]
                conf_source=conf_source.loc[conf_source.pathway.ne("UNMAPPED_PATHWAY")]
            dev=dev_source.dropna(subset=required).copy()
            conf=conf_source.dropna(subset=required).copy()
            rows,preds=evaluate_models(dev,conf,dev_panel,scope,grouped_folds[(dev_panel,scope)],confirmation_panel=confirm_panel)
            model_rows.extend(rows); prediction_frames.extend(preds)

    # Full and mapped-only pathway summaries use each panel's own finite population.
    pathway_rows, pathway_assocs = [], []
    for panel, frame in panel_frames.items():
        for scope,pop in {"FULL":frame,"MAPPED_ONLY":frame.loc[frame.pathway.ne("UNMAPPED_PATHWAY")].copy()}.items():
            common=pop.dropna(subset=[*FEATURES,"non_tie_coverage","usefulness","pathway"]).copy()
            sizes=common.groupby("pathway").size()
            for descriptor in FEATURES:
                for pathway,group in common.groupby("pathway",sort=True):
                    pathway_rows.append({"panel":panel,"population_scope":scope,"pathway":pathway,"pathway_type":"UNMAPPED" if pathway=="UNMAPPED_PATHWAY" else "MAPPED_SUBSYSTEM","n_reactions":len(group),
                        "descriptor":descriptor,"spearman":spearman(group[descriptor],group.usefulness) if len(group)>=5 else np.nan,"sufficient_n_ge_5":len(group)>=5})
                path_summary=common.groupby("pathway",sort=True).agg(**{descriptor:(descriptor,"mean"),"usefulness_mean":("usefulness","mean"),"n_reactions":("reaction_id","size")}).reset_index()
                rho=spearman(path_summary[descriptor],path_summary.usefulness_mean)
                rlo,rhi,rn,slo,shi,sn=pathway_block_bootstraps(common,descriptor,"usefulness",seed=SEED+list(PANELS).index(panel)*10+FEATURES.index(descriptor)+(100 if scope=="MAPPED_ONLY" else 0))
                mapped=common.loc[common.pathway.ne("UNMAPPED_PATHWAY")]
                pathway_assocs.append({"panel":panel,"population_scope":scope,"descriptor":descriptor,"n_reactions":len(common),"n_pathways":len(path_summary),"n_mapped_pathways":int(path_summary.pathway.ne("UNMAPPED_PATHWAY").sum()),"n_unmapped_reactions":int(common.pathway.eq("UNMAPPED_PATHWAY").sum()),
                    "min_pathway_size":int(sizes.min()),"median_pathway_size":float(sizes.median()),"max_pathway_size":int(sizes.max()),
                    "reaction_level_spearman":spearman(common[descriptor],common.usefulness),"reaction_level_pathway_block_bootstrap_q025":rlo,"reaction_level_pathway_block_bootstrap_q975":rhi,"finite_reaction_bootstrap_replicates":rn,
                    "mapped_only_reaction_level_spearman":spearman(mapped[descriptor],mapped.usefulness),
                    "pathway_summary_spearman_equal_pathway":rho,"pathway_summary_bootstrap_q025":slo,"pathway_summary_bootstrap_q975":shi,"finite_pathway_bootstrap_replicates":sn,
                    "interpretation":"post-hoc pathway-block stability; pathway grouping does not guarantee independence"})

    # Paired correlation differences within each panel on the exact common cases.
    strength_rows=[]
    for panel, frame in panel_frames.items():
        for scope,pop in {"FULL":frame,"MAPPED_ONLY":frame.loc[frame.pathway.ne("UNMAPPED_PATHWAY")].copy()}.items():
            common=pop.dropna(subset=[*FEATURES,"non_tie_coverage","usefulness"]).copy()
            for competitor in ("H_dir","D_dir"):
                lo,hi,nboot=bootstrap_corr_delta(common,"eta2",competitor,"usefulness",seed=SEED+100+list(PANELS).index(panel)+(100 if scope=="MAPPED_ONLY" else 0),replicates=BOOTSTRAP_REPLICATES)
                paired_corr.append({"split":PANELS[panel][1],"panel":panel,"population_scope":scope,"variable":f"signed_rho_eta2_minus_{competitor}","n_matched_reactions":len(common),"signed_correlation_difference_eta2_minus_competitor":spearman(common.eta2,common.usefulness)-spearman(common[competitor],common.usefulness),"paired_bootstrap_q025":lo,"paired_bootstrap_q975":hi,"finite_bootstrap_replicates":nboot})
            lo,hi,nboot,diff=bootstrap_strength_delta(common,"eta2","D_dir","usefulness",seed=SEED+200+list(PANELS).index(panel)+(100 if scope=="MAPPED_ONLY" else 0),replicates=BOOTSTRAP_REPLICATES)
            rho_eta=spearman(common.eta2,common.usefulness); rho_d=spearman(common.D_dir,common.usefulness)
            strength_rows.append({"panel":panel,"anchor":PANELS[panel][0],"split":PANELS[panel][1],"population_scope":scope,"n_matched_reactions":len(common),
                "rho_eta2_signed":rho_eta,"rho_D_dir_signed":rho_d,"rho_D_dir_sign_aligned":-rho_d,
                "abs_rho_eta2":abs(rho_eta),"abs_rho_D_dir":abs(rho_d),"strength_difference_abs_rho_eta2_minus_abs_rho_D_dir":diff,
                "paired_bootstrap_q025":lo,"paired_bootstrap_q975":hi,"finite_bootstrap_replicates":nboot})

    # Hashes of exact frozen input bytes; record deferred input needs without
    # opening anything outside this repository.
    stats_path=ROOT/"tables/manuscript/fig2_panel_statistics.tsv"
    hashes = {str(p.relative_to(ROOT)): sha256(p) for p in [FIG2, FIG3, stats_path, *sorted(FIG4_DIR.glob("algorithm=*/part.parquet.xz"))]}
    manifest = {"analysis":"reviewer_checks_v2","seed":SEED,"reaction_bootstrap_replicates":BOOTSTRAP_REPLICATES,"pathway_bootstrap_replicates":PATHWAY_BOOTSTRAP_REPLICATES,
                "model_specification":{"estimator":"Ridge","alpha":1.0,"folds":5,"splitter":"GroupKFold(pathway)","preprocessing":"StandardScaler fit on training folds only","entropy_polynomial":"H_dir, H_dir^2, H_dir^3"},
                "input_sha256":hashes,"populations":{"fig2_panel_rows":{p:int((points.panel==p).sum()) for p in PANELS},"matched_truth_pairs_for_A1_A2":836,"matched_truth_pair_usefulness_evaluations_per_anchor":370,
                    "pathway_counts_by_panel":{p:{"all":int(panel_frames[p].pathway.nunique()),"mapped":int(panel_frames[p].loc[panel_frames[p].pathway.ne("UNMAPPED_PATHWAY"),"pathway"].nunique()),"unmapped_reactions":int(panel_frames[p].pathway.eq("UNMAPPED_PATHWAY").sum())} for p in PANELS}},
                "deferred_experiments":{"3":"Original whole candidate vectors and exact frozen candidate identities are absent; required source artifacts include the frozen 20-vector ensemble flux cache plus candidate ID/weight/truth-pair registries used by PL2A/A22.",
                                         "4":"Original candidate observables, per-mouse ranks, and 320-candidate weighting inputs/registries are absent; do not infer recalibrated weights from compact summaries."},
                "caveats":["post-hoc robustness analysis","400 evaluations are computational, not biological replicates","pathway grouping does not guarantee independence","A2 is post-freeze robustness, not independent confirmation"]}
    outputs = {
        "parity_checks.tsv.gz":pd.DataFrame(parity),
        "reaction_level_analysis.tsv.gz":pd.concat([
            frame.assign(panel=panel,anchor=PANELS[panel][0],split=PANELS[panel][1])
            for panel,frame in panel_frames.items()
        ],ignore_index=True)[["panel","anchor","split","reaction_id","pathway",*FEATURES,"non_tie_coverage","usefulness","n_descriptor_evaluations","n_directional_usefulness_evaluations","n_non_tie_cases"]],
        "eligible_evaluations.tsv.gz":gain_raw.groupby("anchor_setting").evaluation_id.unique().explode().rename("evaluation_id").reset_index().rename(columns={"anchor_setting":"source_anchor_label"}),
        "matched_descriptor_comparisons.tsv.gz":pd.DataFrame(association_rows),
        "paired_anchor_comparisons.tsv.gz":pd.DataFrame(paired_corr),
        "predictive_models.tsv.gz":pd.DataFrame(model_rows),
        "pathway_within_associations.tsv.gz":pd.DataFrame(pathway_rows),
        "pathway_blocked_robustness.tsv.gz":pd.DataFrame(pathway_assocs),
        "fold_assignments.tsv.gz":pd.DataFrame(fold_records),
        "correlation_strength_comparisons.tsv.gz":pd.DataFrame(strength_rows),
    }
    for name, table in outputs.items():
        table.to_csv(OUT/name,sep="\t",index=False,compression="gzip",float_format="%.12g")
    if prediction_frames:
        pd.concat([x for x in prediction_frames if isinstance(x,pd.DataFrame)],ignore_index=True).to_csv(OUT/"development_oof_predictions.tsv.gz",sep="\t",index=False,compression="gzip",float_format="%.12g")
    (OUT/"manifest.json").write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":"complete","output":str(OUT.relative_to(ROOT)),"parity":parity,"tables":list(outputs)},indent=2))


if __name__ == "__main__":
    main()
