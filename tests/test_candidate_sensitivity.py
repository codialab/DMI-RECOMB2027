import importlib.util
import os
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
SOURCE_ROOT = Path(os.environ['MRI_BRAIN_TUMOR_ROOT']).resolve() if os.environ.get('MRI_BRAIN_TUMOR_ROOT') else None
spec = importlib.util.spec_from_file_location('candidate_sensitivity', ROOT / 'analyses/reviewer_checks/candidate_sensitivity.py')
cs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cs)


def test_vectorized_descriptors_match_frozen_pl1_core():
    if SOURCE_ROOT is None or not (SOURCE_ROOT / 'scripts/dmi_bridge_pl1_predictability_core_v1.py').is_file():
        pytest.skip('set MRI_BRAIN_TUMOR_ROOT to run upstream parity checks')
    sys.path.insert(0, str(SOURCE_ROOT / 'scripts'))
    import dmi_bridge_pl1_predictability_core_v1 as pl1
    left = np.array([[0.0, 1.0, 3.0], [1.0, 2.0, 2.5], [0.0, -1.0, 2.0]])
    right = np.array([[0.0, 2.0, 1.0], [1.0, 1.0, 0.5]])
    wl = np.array([0.2, 0.3, 0.5]); wr = np.array([0.7, 0.3])
    got = cs.descriptor_block(left, right, wl, wr)
    for j in range(left.shape[1]):
        expected = pl1.target_landscape_metrics(
            anchor_ct2a=left[:, j], anchor_gl261=right[:, j],
            target_ct2a=left[:, j], target_gl261=right[:, j],
            strong_weights_ct2a=wl, strong_weights_gl261=wr)
        assert got['eta2'][j] == pytest.approx(expected['sign_magnitude_eta2'], abs=1e-12)
        assert got['H_dir'][j] == pytest.approx(expected['directional_entropy3'], abs=1e-12)
        assert got['D_dir'][j] == pytest.approx(expected['dominant_sign_mass'], abs=1e-12)


def test_vectorized_weak_cue_gains_match_frozen_scalar_core():
    if SOURCE_ROOT is None or not (SOURCE_ROOT / 'scripts/dmi_bridge_pl2_sign_only_core_v1.py').is_file():
        pytest.skip('set MRI_BRAIN_TUMOR_ROOT to run upstream parity checks')
    sys.path.insert(0, str(SOURCE_ROOT / 'scripts'))
    import dmi_bridge_pl2_sign_only_core_v1 as pl2
    left = np.array([[1.0, 3.0], [2.0, 1.0], [4.0, 2.0]])
    right = np.array([[0.0, 1.0], [2.0, 0.5]])
    wl = np.array([0.2, 0.3, 0.5]); wr = np.array([0.6, 0.4])
    truth = np.array([1.5, -2.0])
    got = cs.vector_outcomes(left, right, wl, wr, truth)
    for j in range(2):
        delta = (left[:, j, None] - right[None, :, j]).reshape(-1)
        joint = (wl[:, None] * wr[None, :]).reshape(-1)
        expected = pl2.evaluate_truth_case(delta_b=delta, baseline_weights=joint,
                                           truth_delta_b=truth[j], lambda_grid=(0.25,))
        row = expected['lambda_results'][0]
        assert got['correct_gain'][j] == pytest.approx(row['correct_abs_error_gain'], abs=1e-12)
        assert got['wrong_gain'][j] == pytest.approx(row['wrong_abs_error_gain'], abs=1e-12)


def test_truth_exclusion_operates_on_whole_vectors_and_renormalizes():
    weights = np.array([0.1, 0.2, 0.3, 0.4])
    sample_ids = np.array([10, 11, 12, 13])
    keep = sample_ids != 12
    remaining = weights[keep] / weights[keep].sum()
    assert sample_ids[keep].tolist() == [10, 11, 13]
    assert remaining.sum() == pytest.approx(1.0)
    assert remaining.tolist() == pytest.approx([1/7, 2/7, 4/7])


def test_global_pool_temperature_operator_reproduces_target_ess():
    if SOURCE_ROOT is None or not (SOURCE_ROOT / 'scripts/dmi_bridge_a20_dual_anchor_core_v1.py').is_file():
        pytest.skip('set MRI_BRAIN_TUMOR_ROOT to run upstream parity checks')
    sys.path.insert(0, str(SOURCE_ROOT / 'scripts'))
    import dmi_bridge_a20_dual_anchor_core_v1 as a20
    distance = np.linspace(0.0, 0.2, 320) ** 2
    for target in (10.0, 20.0, 40.0):
        result = a20.solve_temperature_for_ess(distance, target_ess=target)
        assert result.achieved_ess == pytest.approx(target, abs=1e-9)
        assert result.weights.sum() == pytest.approx(1.0, abs=1e-14)


def test_two_stage_aggregation_equal_evaluation_weighting_and_ties():
    import pandas as pd
    rows = []
    def add(evaluation, pair, value, tie=False):
        rows.append({'scenario':'s','arm':'A1','target_ess':20.,'evaluation_id':evaluation,
            'reaction_id':'R','truth_pair_id':pair,'truth_tie':tie,
            'usefulness':np.nan if tie else value,'correct_gain':np.nan if tie else value,
            'wrong_gain':np.nan if tie else 0.,'information_advantage':np.nan if tie else value/2,
            'correct_minus_wrong':np.nan if tie else value})
    add('e1','p1',1.); add('e1','p2',3.); add('e1','tie',0.,True)
    add('e2','p3',9.)
    # Duplicate pair identities are collapsed before within-evaluation means.
    add('e2','p3',99.)
    evals, reaction = cs.aggregate_evaluation_reaction_records(pd.DataFrame(rows))
    got = reaction[0]
    assert got['usefulness'] == pytest.approx(((1.+3.)/2 + 9.)/2)
    assert got['usefulness'] != pytest.approx((1.+3.+9.)/3)
    assert got['n_eligible_evaluations'] == 2
    assert got['n_defined_response_evaluations'] == 2
    assert got['n_distinct_truth_pairs'] == 4
    assert got['n_non_tie_truth_pairs'] == 3
    assert got['n_truth_ties'] == 1
    assert evals.loc[evals.evaluation_id.eq('e1'),'usefulness'].iloc[0] == pytest.approx(2.)


def test_fixed_common_ess_population_is_exact_intersection():
    common = {10:{'a','b','c'},20:{'b','c','d'},40:{'b','c','e'}}
    assert cs.fixed_common_population(common) == {'b','c'}


def test_frozen_pl2_pair_operator_parity_on_synthetic_records():
    if SOURCE_ROOT is None or not (SOURCE_ROOT / 'scripts/dmi_bridge_pl2c_geometry_utility_core_v1.py').is_file():
        pytest.skip('set MRI_BRAIN_TUMOR_ROOT to run upstream PL2 aggregation parity')
    sys.path.insert(0, str(SOURCE_ROOT / 'scripts'))
    import dmi_bridge_pl2c_geometry_utility_core_v1 as pl2c
    source=[{'truth_pair_id':'p1','truth_status':'NON_TIE','correct_gain':0.4,'wrong_gain':0.1,'random_gain':0.25},
            {'truth_pair_id':'p2','truth_status':'NON_TIE','correct_gain':-0.1,'wrong_gain':-0.3,'random_gain':-0.2},
            {'truth_pair_id':'pt','truth_status':'TRUTH_TIE','correct_gain':np.nan,'wrong_gain':np.nan,'random_gain':np.nan}]
    expected=pl2c.aggregate_distinct_pair_responses(source)
    records=[]
    for x in source:
        tie=x['truth_status']=='TRUTH_TIE'
        info=(x['correct_gain']-x['wrong_gain'])/2 if not tie else np.nan
        records.append({'scenario':'s','arm':'A1','target_ess':20.,'evaluation_id':'e','reaction_id':'R',
            'truth_pair_id':x['truth_pair_id'],'truth_tie':tie,'usefulness':float(x['correct_gain']>1e-12 and info>1e-12) if not tie else np.nan,
            'correct_gain':x['correct_gain'],'wrong_gain':x['wrong_gain'],'information_advantage':info,
            'correct_minus_wrong':x['correct_gain']-x['wrong_gain'] if not tie else np.nan})
    ev,_=cs.aggregate_evaluation_reaction_records(__import__('pandas').DataFrame(records))
    row=ev.iloc[0]
    assert row.n_distinct_truth_pairs==expected['n_distinct_pairs']
    assert row.n_non_tie_truth_pairs==expected['n_non_tie_pairs']
    assert row.n_truth_ties==expected['n_truth_tie_pairs']
    assert row.correct_gain_mean==pytest.approx(expected['mean_correct_gain'],abs=1e-15)
    assert row.information_advantage_mean==pytest.approx(expected['mean_information_advantage'],abs=1e-15)
    assert row.usefulness==pytest.approx(expected['directionally_useful_fraction'],abs=1e-15)


def test_development_confirmation_prediction_never_trains_on_confirmation():
    import pandas as pd
    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import Ridge
    scaler_rows = []
    ridge_fits = []
    original_fit = StandardScaler.fit
    original_ridge_fit = Ridge.fit
    def tracked_fit(self, X, y=None, **kwargs):
        scaler_rows.append(len(X))
        return original_fit(self, X, y, **kwargs)
    def tracked_ridge_fit(self, X, y, **kwargs):
        ridge_fits.append((len(X), self.alpha))
        return original_ridge_fit(self, X, y, **kwargs)
    StandardScaler.fit = tracked_fit
    Ridge.fit = tracked_ridge_fit
    try:
        ids=[f'R{i}' for i in range(40)]
        outcome=pd.DataFrame({'arm':'A1','reaction_id':ids,'usefulness':(np.arange(40)%5)/5})
        geometry=pd.DataFrame({'arm':'A1','reaction_id':ids,'eta2':np.arange(40)%7,'H_dir':np.arange(40)%9})
        registry=pd.DataFrame({'reaction_id':ids,'subsystem':[f'path{i%6}' for i in range(40)]})
        split=pd.DataFrame({'reaction_id':ids,'analysis_split':['DEVELOPMENT']*30+['CONFIRMATION_HOLDOUT']*10})
        got=cs.associations(outcome,geometry,geometry,registry,split)
    finally:
        StandardScaler.fit = original_fit
        Ridge.fit = original_ridge_fit
    assert set(got.evaluation)=={'development_pathway_grouped_cv','confirmation_from_development_fit'}
    assert 30 in scaler_rows
    assert 40 not in scaler_rows  # scaling never sees confirmation rows
    assert ridge_fits and 40 not in [n for n, _ in ridge_fits]
    assert all(n <= 30 and alpha == 1.0 for n, alpha in ridge_fits)


def test_pathway_group_folds_are_isolated():
    import pandas as pd
    from sklearn.model_selection import GroupKFold
    groups=pd.Series([f'path{i%6}' for i in range(30)])
    folds=GroupKFold(5).split(np.zeros((30,2)),np.zeros(30),groups)
    assert all(set(groups.iloc[tr]).isdisjoint(set(groups.iloc[te])) for tr,te in folds)


def test_truth_exclusion_removes_whole_vector_and_renormalizes_weights():
    weights = np.array([0.1, 0.2, 0.3, 0.4]); flux = np.arange(12).reshape(4,3)
    keep = np.array([0,1,3]); retained = flux[keep]
    retained_weights = weights[keep] / weights[keep].sum()
    assert retained.tolist() == flux[[0,1,3]].tolist()
    assert retained_weights.sum() == pytest.approx(1.)
    assert retained_weights.tolist() == pytest.approx([1/7,2/7,4/7])


def test_saved_frozen_baseline_and_ess20_parity_gates():
    """Use local corrected result tables when present; skip on a clean clone."""
    outputs = ROOT / 'reproduced/reviewer_checks/corrected'
    a3 = outputs / 'experiment_3/experiment_3a/experiment3a_recomputed_vs_frozen_case_audit.tsv'
    a4 = outputs / 'experiment_4/experiment4_ess20_outcome_parity.tsv'
    if not a3.is_file() or not a4.is_file():
        pytest.skip('corrected upstream result tables are local ignored outputs; run Experiments 3A and 4 to enable baseline parity integration checks')
    import pandas as pd
    audit3 = pd.read_csv(a3, sep='\t')
    endpoint = audit3.loc[audit3.field.eq('usefulness')]
    assert set(endpoint.arm) == {'A1', 'A2'}
    assert endpoint.n_different.eq(0).all()
    assert endpoint.reference.isin(['frozen PL2B', 'frozen A22']).all()
    audit4 = pd.read_csv(a4, sep='\t')
    assert audit4.status.eq('PASS').all()
    assert audit4.field.isin(['n_distinct_truth_pairs','n_non_tie_truth_pairs','n_truth_ties',
                              'usefulness','correct_gain_mean','wrong_gain_mean',
                              'information_advantage_mean']).all()
