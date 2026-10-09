import numpy as np
import pandas as pd
import pytest
from sklearn.model_selection import GroupKFold

from analyses.reviewer_checks import run_reviewer_checks as rc
from scripts import dmi_bridge_pl1_predictability_core_v1 as pl1


def test_descriptor_formulas_match_frozen_core():
    # Weighted sign groups have distinct magnitude means, so eta2 is defined.
    delta=np.array([-3.,-1.,2.,4.])
    weights=np.array([.25,.25,.25,.25])
    eta=pl1.sign_magnitude_eta2(delta,weights)
    assert eta == pytest.approx(.2)
    metrics=pl1.target_landscape_metrics(
        anchor_ct2a=[0.,1.],anchor_gl261=[0.,0.],
        target_ct2a=[-1.,2.],target_gl261=[0.,0.],
        strong_weights_ct2a=[.5,.5],strong_weights_gl261=[.5,.5],
    )
    assert 0. <= metrics["sign_magnitude_eta2"] <= 1.
    assert 0. <= metrics["directional_entropy3"] <= np.log(3.)
    assert metrics["n_supported_sign_states"] in (1,2,3)


def test_truth_holdout_renormalizes_both_condition_weights():
    left,right=rc.holdout_product_weights([1,2,1],[2,1,1],1,0)
    assert left == pytest.approx([.5,0.,.5])
    assert right == pytest.approx([0.,.5,.5])
    assert left.sum()==pytest.approx(1.) and right.sum()==pytest.approx(1.)
    with pytest.raises(ValueError):
        rc.holdout_product_weights([1,0],[1,1],0,0)


def test_reaction_aggregation_keeps_panel_reaction_means_and_pathway():
    frame=pd.DataFrame({
        "anchor":["A1"]*4,"evaluation_id":["E1","E1","E2","E2"],
        "reaction_id":["R1","R2","R1","R2"],"subsystem":["P","P","P","P"],
        "eta2":[.2,.4,.6,.8],"H_dir":[.1,.2,.3,.4],"D_dir":[.9,.8,.7,.6],
        "K_dir":[2,2,3,3],"non_tie_coverage":[1,1,.5,.5],
    })
    out=rc.aggregate_geometry(frame).set_index("reaction_id")
    assert out.loc["R1","eta2"]==pytest.approx(.4)
    assert out.loc["R2","H_dir"]==pytest.approx(.3)
    assert out.loc["R1","n_descriptor_evaluations"]==2
    assert out.loc["R1","pathway"]=="P"


def test_bootstrap_is_deterministic_and_paired():
    frame=pd.DataFrame({"eta2":[1,2,3,4,5],"H_dir":[5,4,3,2,1],"usefulness":[1,2,2,4,5]})
    a=rc.bootstrap_corr_delta(frame,"eta2","H_dir","usefulness",seed=4,replicates=100)
    b=rc.bootstrap_corr_delta(frame,"eta2","H_dir","usefulness",seed=4,replicates=100)
    assert a==b
    assert 0 < a[2] <= 100


def test_grouped_folds_never_split_pathways():
    frame=pd.DataFrame({"pathway":["A","A","B","B","C","C","D","D","E","E"],"usefulness":np.arange(10.)})
    folds=rc.make_group_folds(frame,n_splits=5)
    for train,test in folds:
        assert not set(frame.pathway.iloc[train]) & set(frame.pathway.iloc[test])
    assert sorted(np.concatenate([test for _,test in folds]).tolist())==list(range(len(frame)))


def test_equal_pathway_bootstrap_repeats_with_seed():
    frame=pd.DataFrame({"pathway":["A","A","B","B","C","C"],"x":[1,2,3,4,5,6],"y":[2,1,4,3,6,5]})
    a=rc.pathway_block_bootstraps(frame,"x","y",seed=12)
    b=rc.pathway_block_bootstraps(frame,"x","y",seed=12)
    assert a==b
    assert a[2]==rc.PATHWAY_BOOTSTRAP_REPLICATES
    assert 0 < a[5] <= rc.PATHWAY_BOOTSTRAP_REPLICATES


def test_dominant_mass_strength_uses_absolute_correlation():
    frame=pd.DataFrame({"eta2":np.arange(1,9,dtype=float),"D_dir":np.arange(8,0,-1,dtype=float),"usefulness":np.arange(1,9,dtype=float)})
    lo,hi,n,point=rc.bootstrap_strength_delta(frame,"eta2","D_dir","usefulness",seed=7,replicates=100)
    assert point==pytest.approx(0.)
    assert lo==pytest.approx(0.) and hi==pytest.approx(0.)
    assert n==100


def test_nested_models_share_population_and_grouped_folds():
    rows=[]
    rng=np.random.default_rng(8)
    for group in range(6):
        for i in range(8):
            h=float(i+group)/14
            rows.append({"pathway":f"P{group}","reaction_id":f"R{group}_{i}","H_dir":h,"eta2":h*h,
                         "D_dir":1-h,"K_dir":2,"non_tie_coverage":h,"usefulness":.1+.6*h+rng.normal(0,.01)})
    frame=pd.DataFrame(rows)
    folds=rc.make_group_folds(frame,n_splits=5)
    result,_=rc.evaluate_models(frame.copy(),None,"D","FULL",folds)
    assert len({row["n"] for row in result})==1
    assert {row["model"] for row in result}>={"N1_entropy_cubic","N2_entropy_cubic_plus_eta2","C1_entropy_coverage","C2_entropy_coverage_plus_eta2"}
    assert all(np.isfinite(row["mae"]) for row in result)


def test_unmapped_pathway_filter_preserves_mapped_rows_only():
    frame=pd.DataFrame({"pathway":["A","UNMAPPED_PATHWAY","B"],"reaction_id":["1","2","3"]})
    mapped=frame.loc[frame.pathway.ne("UNMAPPED_PATHWAY")]
    assert mapped.reaction_id.tolist()==["1","3"]
