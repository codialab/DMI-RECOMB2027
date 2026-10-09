#!/usr/bin/env python3
"""Post-hoc candidate-level sensitivity analyses against a read-only MRI checkout.

All source artifacts are verified before use. No source-root files are written.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, lzma, math, os, sys, subprocess
from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd

EXPECTED = {
 'outputs/secondary_anchor_exploratory/dmi_sa_six_arm_current_20260926/reconstructed/flux_cache.npz':'f9644da4ebd875cf5f0b34cc58b20faf0a0457d7c7a44635e8a010d8896c3dbb',
 '12_recomb_method_comparison/outputs_dmi_fraction_ensemble_v2/candidate_fraction_table.tsv.xz':'421c2c1d196fbd088a4e92e2ecc998ad578aab814a025fff8e1c843597fed1da',
 'outputs/dmi_bridge_bio0_audit_v1/BRIDGEBIO0_FOUR_ARM_WEIGHTS.tsv.xz':'843e17b7012a6fa212d1e9ac9e1dc6e309e854d94ab5c4b3f9e316b589768f4f',
 'outputs/secondary_anchor_exploratory/dmi_sa_six_arm_current_20260926/current_candidate_observables.csv':'946cd613d673c9aa9367bce705ad5d05562c6a174acc830fbbd36cb395aa5504',
 'outputs/dmi_bridge_a20_dual_anchor_qualification_v1/BRIDGEA20_DUAL_ANCHOR_WEIGHTS.tsv.xz':'b153e2fde45dca9573c71c548783e704313f056d0980f64af60180ad90ebc1f4',
 'outputs/dmi_bridge_a20_dual_anchor_qualification_v1/BRIDGEA20_MOUSE_COORDINATES.tsv':'f41d0fef703b52201de7ed8a86d888bdca28d5262ce52145bf8e53b00a81e944',
 'outputs/dmi_bridge_pl2a_sign_only_prepare_v1/BRIDGEPL2A_TRUTH_SELECTION.tsv':'a3906367a89da57e471f5c9e7733c363e4e71b4715d5dd9c9309a46d9f8f630f',
 'outputs/dmi_bridge_pl2a_sign_only_prepare_v1/BRIDGEPL2A_TRUTH_PAIR_REGISTRY.tsv':'2986c98deb6563e2c4f192ba1f77bbf0385976b42266363bcb968ca6d0cd3894',
 'outputs/dmi_bridge_pl1_predictability_landscape_v1_resumable_v1/BRIDGEPL1_EVALUATION_REGISTRY.tsv':'aa468dcec11dea764682bd3d7eb14c03c003063e74d92080d0bfaabb366290a3',
 'outputs/dmi_bridge_pl2c_geometry_utility_development_v1/BRIDGEPL2C_REACTION_SPLIT.tsv':'7dd8676508416a13f88af184110ef40e0f8d92782cd5d62b42e7006af5316e2f',
 'outputs/dmi_bridge_a23_a1_a2_synthesis_v1/BRIDGEA23_MATCHED_POPULATION_AUDIT.tsv':'9957482aec24f4b8b128e2990aa1f6aeaf507cce83600d556ec07fa4cb885d60',
 'outputs/dmi_bridge_pl1_predictability_landscape_v1_resumable_v1/BRIDGEPL1_MANIFEST.json':'b0f81620ff24f5791dd964a8063416784db86a358b0464b8c68a504e84e5f529',
}
A1_CURRENT = '13_recomb_identifiability/outputs/sa_u1_u2_followup_20260926/current_vmax_base_weights.tsv.xz'
PL1_MANIFEST='outputs/dmi_bridge_pl1_predictability_landscape_v1_resumable_v1/BRIDGEPL1_MANIFEST.json'
A21_MANIFEST='outputs/dmi_bridge_a21_dual_anchor_geometry_v1/BRIDGEA21_MANIFEST.json'
A21_MANIFEST_SHA='db9948d7b99447a695cf8aa21429dde6d642e1fd0632bcce52643deac59e6feb'
A1_CURRENT_SHA = '91121bf1251354505cde5f5ae783dc8d14bc9a22ecdf8daa06be1e5d82e3385d'
A1_MANIFEST = '13_recomb_identifiability/outputs/sa_u1_u2_followup_20260926/current_vmax_base_weights_manifest.json'
A1_MANIFEST_SHA = '2abf9e7c4d3c922bda8c81426001c1f27c6e43f1761e9c9073e0d8ca078bd5bc'
ARCHIVED_KERNEL='outputs/secondary_anchor_exploratory/dmi_sa_six_arm_current_20260926/reconstruct_flux_cache.py'
ARCHIVED_KERNEL_SHA='70c77a76a863a7ca800b7c3610175417b2720c5cbc3d4c1b6db331ec3cef1d32'
CODE_PATHS=['scripts/dmi_bridge_pl1_predictability_core_v1.py','scripts/dmi_bridge_pl2_sign_only_core_v1.py','scripts/dmi_bridge_pl2a_prepare_v1.py','scripts/dmi_bridge_pl2b_production_v1.py','scripts/dmi_bridge_a20_dual_anchor_core_v1.py','scripts/dmi_bridge_a21_dual_anchor_geometry_core_v1.py','scripts/dmi_bridge_a22_dual_anchor_sign_only_core_v1.py','scripts/dmi_bridge_a23_a1_a2_synthesis_core_v1.py']
TRUTH_SEEDS = (20271028,20271029,20271030,20271031)
LAM = .25
ATOL = 1e-12
KEY = ['mouse_id','tumor','algorithm','ensemble_hash','rna_context_key']
CAND_KEY = ['algorithm','tumor','ensemble_hash','sample_index']
IDENTITY = ['algorithm','tumor','ensemble_hash','projection_hash','rna_context_key','rna_training_samples','sample_index']


def sha256(p:Path)->str:
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): h.update(b)
 return h.hexdigest()

def checked(root:Path, rel:str, expected:str)->Path:
 p=root/rel
 if not p.is_file() or p.is_symlink(): raise RuntimeError(f'missing/non-regular frozen source: {rel}')
 got=sha256(p)
 if got!=expected: raise RuntimeError(f'input SHA256 mismatch for {rel}: {got} != {expected}')
 return p

def load(root:Path):
 actual=Path(subprocess.run(['git','-C',str(root),'rev-parse','--show-toplevel'],check=True,capture_output=True,text=True).stdout.strip()).resolve()
 if actual!=root.resolve(): raise RuntimeError('upstream root must be the named MRI repository root')
 paths={r:checked(root,r,h) for r,h in EXPECTED.items()}
 checked(root,A1_CURRENT,A1_CURRENT_SHA); checked(root,A1_MANIFEST,A1_MANIFEST_SHA); checked(root,A21_MANIFEST,A21_MANIFEST_SHA); checked(root,ARCHIVED_KERNEL,ARCHIVED_KERNEL_SHA)
 import dmi_bridge_pl1_predictability_core_v1 as pl1
 import dmi_bridge_pl2_sign_only_core_v1 as pl2
 import dmi_bridge_a20_dual_anchor_core_v1 as a20
 panel=pd.read_csv(paths['12_recomb_method_comparison/outputs_dmi_fraction_ensemble_v2/candidate_fraction_table.tsv.xz'],sep='\t',compression='xz',float_precision='round_trip')
 if len(panel)!=640 or panel.duplicated(CAND_KEY).any(): raise RuntimeError('candidate panel must be 640 unique candidates')
 if sorted(panel.algorithm.unique())!=['CORDA','GIMME','RIPTiDe','iMAT']: raise RuntimeError('method panel mismatch')
 strata=panel.groupby(['algorithm','tumor','ensemble_hash'],sort=False).size()
 if len(strata)!=32 or not strata.eq(20).all(): raise RuntimeError('candidate panel is not exactly 32 strata of 20 vectors')
 cache=np.load(paths['outputs/secondary_anchor_exploratory/dmi_sa_six_arm_current_20260926/reconstructed/flux_cache.npz'],allow_pickle=False)
 mats=cache['mats']
 if mats.shape!=(32,20,4181) or not np.isfinite(mats).all(): raise RuntimeError('frozen flux cache shape/finiteness mismatch')
 rxns=list(map(str,cache['rxns'].tolist()))
 if len(set(rxns))!=4181: raise RuntimeError('duplicate reaction identity')
 # Input raw glucose/lactate operators from the original candidate panel; cache is flux-only.
 panel['candidate_Vmax']=np.maximum(-panel.EX_glc__D_e.to_numpy(float),0.)
 panel['candidate_Vlac']=np.maximum(-panel.LDH_L.to_numpy(float),0.)
 observ=pd.read_csv(paths['outputs/secondary_anchor_exploratory/dmi_sa_six_arm_current_20260926/current_candidate_observables.csv'],float_precision='round_trip')
 j=panel[CAND_KEY+['candidate_Vmax','candidate_Vlac']].merge(observ[CAND_KEY+['Vmax','Vlac']],on=CAND_KEY,validate='one_to_one')
 if len(j)!=640 or np.max(np.abs(j.candidate_Vmax-j.Vmax))>ATOL or np.max(np.abs(j.candidate_Vlac-j.Vlac))>ATOL:
  raise RuntimeError('raw candidate-panel observables do not reproduce frozen observable table')
 a1=pd.read_csv(root/A1_CURRENT,sep='\t',compression='xz',float_precision='round_trip')
 a1bio=pd.read_csv(paths['outputs/dmi_bridge_bio0_audit_v1/BRIDGEBIO0_FOUR_ARM_WEIGHTS.tsv.xz'],sep='\t',compression='xz',float_precision='round_trip')
 a1bio=a1bio.loc[a1bio.arm.eq('strong_anchor_baseline')]
 a2=pd.read_csv(paths['outputs/dmi_bridge_a20_dual_anchor_qualification_v1/BRIDGEA20_DUAL_ANCHOR_WEIGHTS.tsv.xz'],sep='\t',compression='xz',float_precision='round_trip')
 # Use original authoritative A1 and A2 weights and verify against the tables consumed by PL1/A20.
 idw=['mouse_id','algorithm','tumor','ensemble_hash','rna_context_key','sample_index']
 a1b=a1bio[idw+['weight']].merge(a1[idw+['weight']],on=idw,suffixes=('_bio0','_current'),validate='one_to_one')
 if len(a1b)!=3200 or np.max(np.abs(a1b.weight_bio0-a1b.weight_current))>1e-12: raise RuntimeError('A1 baseline weights fail current Vmax parity')
 a1=a1bio.rename(columns={'weight':'weight'})
 a2=a2.loc[a2.arm.eq('A2-L')].copy()
 evals=pd.read_csv(paths['outputs/dmi_bridge_pl1_predictability_landscape_v1_resumable_v1/BRIDGEPL1_EVALUATION_REGISTRY.tsv'],sep='\t',dtype=str)
 truth=pd.read_csv(paths['outputs/dmi_bridge_pl2a_sign_only_prepare_v1/BRIDGEPL2A_TRUTH_SELECTION.tsv'],sep='\t',dtype=str)
 pairs=pd.read_csv(paths['outputs/dmi_bridge_pl2a_sign_only_prepare_v1/BRIDGEPL2A_TRUTH_PAIR_REGISTRY.tsv'],sep='\t',dtype=str)
 if len(evals)!=400 or len(truth)!=480 or len(pairs)!=1200 or int(pairs.pair_evaluable.eq('true').sum())!=1110: raise RuntimeError('frozen registry accounting mismatch')
 return locals()

def weight_map(table:pd.DataFrame, arm:str)->dict:
 d=table.copy()
 if 'arm' in d: d=d[d.arm.eq(arm)]
 out={}
 for key,g in d.groupby(KEY,sort=False):
  g=g.sort_values('sample_index',kind='stable')
  w=g.weight.to_numpy(float)
  if len(w)!=20 or not np.isfinite(w).all() or np.any(w<0) or w.sum()<=0: raise RuntimeError(f'invalid 20-vector support {key}')
  out[tuple(map(str,key))]=(g.sample_index.to_numpy(int),w/w.sum())
 return out

def scalar_outcome(delta, weights, truth, core):
 """Independent check wrapper around the frozen PL2 scalar endpoint."""
 return core.evaluate_truth_case(delta_b=delta,baseline_weights=weights,truth_delta_b=float(truth),lambda_grid=(LAM,),reliability_grid=(.5,))

def vector_outcomes(left:np.ndarray,right:np.ndarray,wl:np.ndarray,wr:np.ndarray,truth_delta:np.ndarray,block:int=64)->dict:
 """Vectorized frozen PL2 absolute-error gains for one case and reaction blocks."""
 n=left.shape[0]; m=right.shape[0]
 joint=(wl[:,None]*wr[None,:]).reshape(-1)
 d=(left[:,None,:]-right[None,:,:]).transpose(2,0,1).reshape(left.shape[1],-1)
 neg=d < -ATOL; pos=d > ATOL; tie=~(neg|pos); mag=np.abs(d)
 masses=[neg@joint,tie@joint,pos@joint]
 abs_sums=[(neg*mag)@joint,(tie*mag)@joint,(pos*mag)@joint]
 sq=[neg@(joint*joint),tie@(joint*joint),pos@(joint*joint)]
 bmean=sum(abs_sums); terr=np.abs(bmean-np.abs(truth_delta))
 correct=np.full(len(truth_delta),np.nan); wrong=np.full_like(correct,np.nan); info=np.full_like(correct,np.nan)
 non=np.abs(truth_delta)>ATOL
 if not non.any(): return {'truth_direction':np.sign(truth_delta),'correct_gain':correct,'wrong_gain':wrong,'information_gain':info,'usefulness':np.full(len(truth_delta),np.nan)}
 hi=math.exp(LAM); lo=math.exp(-LAM)
 for is_correct,dest in [(True,correct),(False,wrong)]:
  obs=np.sign(truth_delta) if is_correct else -np.sign(truth_delta)
  fneg=np.where(obs>0,lo,hi); fpos=np.where(obs>0,hi,lo)
  z=fneg*masses[0]+masses[1]+fpos*masses[2]
  mean=(fneg*abs_sums[0]+abs_sums[1]+fpos*abs_sums[2])/z
  dest[non]=terr[non]-np.abs(mean[non]-np.abs(truth_delta[non]))
 info[non]=(correct[non]-wrong[non])/2.
 util=non&(correct>ATOL)&(info>ATOL)
 return {'truth_direction':np.sign(truth_delta),'correct_gain':correct,'wrong_gain':wrong,'information_gain':info,'usefulness':util.astype(float)}

def candidate_index(data):
 return {(str(a),str(t),str(e)):i for i,(a,t,e) in enumerate(zip(data['cache']['alg'],data['cache']['tumor'],data['cache']['eh']))}

def descriptor_block(left:np.ndarray,right:np.ndarray,wl:np.ndarray,wr:np.ndarray,tol:float=ATOL)->dict[str,np.ndarray]:
 """Vectorized exact sign/magnitude descriptors for a reaction block."""
 d=(left.T[:,:,None]-right.T[:,None,:]).reshape(left.shape[1],-1)
 joint=(wl[:,None]*wr[None,:]).reshape(-1)
 mag=np.abs(d); pos=d>tol; neg=d < -tol; tie=~(pos|neg)
 masks=(neg,tie,pos); masses=[]; sums=[]
 for mask in masks:
  masses.append(mask@joint); sums.append((mask*mag)@joint)
 overall=mag@joint
 total=((mag-overall[:,None])**2)@joint
 between=np.zeros(left.shape[1],float)
 for mass,sm in zip(masses,sums):
  mu=np.divide(sm,mass,out=np.zeros_like(sm),where=mass>0)
  between+=mass*(mu-overall)**2
 eta=np.divide(between,total,out=np.full_like(total,np.nan),where=total>tol)
 eta=np.clip(eta,0.,1.)
 m=np.column_stack(masses); nz=m>0
 ent=np.zeros(len(m),float)
 ent[nz.any(axis=1)]=-np.sum(np.where(nz,m*np.log(np.where(nz,m,1.)),0.),axis=1)[nz.any(axis=1)]/math.log(3.)
 return {'eta2':eta,'H_dir':ent,'D_dir':np.max(m,axis=1),'K_dir':nz.sum(axis=1),'non_tie_coverage':m[:,0]+m[:,2]}

def original_weights(table,mouse,tumor,alg,rna):
 d=table.loc[table.mouse_id.astype(str).eq(str(mouse)) & table.tumor.astype(str).eq(str(tumor)) & table.algorithm.astype(str).eq(str(alg)) & table.rna_context_key.astype(str).eq(str(rna))].sort_values('sample_index',kind='stable')
 if len(d)!=20: return None
 w=d.weight.to_numpy(float)
 return (w/w.sum()) if w.sum()>0 else None

def run3a(data,out,pilot=False):
 """Recompute full geometry for parity, then truth-excluded geometry on matched pairs."""
 root=data['root']; cache=data['cache']; evals=data['evals']; panel=data['panel']; rxns=data['rxns']; rindex={r:i for i,r in enumerate(rxns)}; pl1=data['pl1']
 import dmi_bridge_pl2a_prepare_v1 as pl2a
 target=[r for r in rxns if r not in {'HEX1','LDH_L'}]; tidx=np.asarray([rxns.index(r) for r in target])
 cachemap=candidate_index(data); panelidx={(str(r.algorithm),str(r.tumor),str(r.ensemble_hash),int(r.sample_index)):i for i,r in enumerate(panel.itertuples())}
 # Full original weighted 20x20 geometry parity against PL1 and A21 reaction summaries.
 summaries={
  'A1':pd.read_csv(root/'outputs/dmi_bridge_pl1_predictability_landscape_v1_resumable_v1/BRIDGEPL1_REACTION_SUMMARY.tsv',sep='\t'),
  'A2':pd.read_csv(root/'outputs/dmi_bridge_a21_dual_anchor_geometry_v1/BRIDGEA21_REACTION_SUMMARY.tsv.xz',sep='\t',compression='xz')}
 saved_hashes={}
 for arm,rel,manifest in [('A1','outputs/dmi_bridge_pl1_predictability_landscape_v1_resumable_v1/BRIDGEPL1_REACTION_SUMMARY.tsv',PL1_MANIFEST),('A2','outputs/dmi_bridge_a21_dual_anchor_geometry_v1/BRIDGEA21_REACTION_SUMMARY.tsv.xz',A21_MANIFEST)]:
  mf=json.loads((root/manifest).read_text()); expected=mf['artifact_sha256'][Path(rel).name]
  p=root/rel
  if sha256(p)!=expected: raise RuntimeError(f'{arm} frozen geometry summary hash mismatch')
  saved_hashes[arm]=expected
 baseline={a:{k:np.full((400,4179),np.nan) for k in ('eta2','H_dir','D_dir')} for a in ('A1','A2')}
 rowpos={str(r.evaluation_id):i for i,r in enumerate(evals.itertuples())}
 for ev in evals.itertuples(index=False):
  eidx=rowpos[str(ev.evaluation_id)]; alg=str(ev.algorithm); rna=str(ev.rna_context_key)
  for arm,wt in [('A1',data['a1bio']),('A2',data['a2'])]:
   wc=original_weights(wt,str(ev.ct2a_mouse),'CT2A',alg,rna); wg=original_weights(wt,str(ev.gl261_mouse),'GL261',alg,rna)
   if wc is None or wg is None: raise RuntimeError('full geometry support absent')
   ensC=str(panel.loc[panel.algorithm.eq(alg)&panel.tumor.eq('CT2A')&panel.rna_context_key.eq(rna),'ensemble_hash'].iloc[0])
   ensG=str(panel.loc[panel.algorithm.eq(alg)&panel.tumor.eq('GL261')&panel.rna_context_key.eq(rna),'ensemble_hash'].iloc[0])
   lc=cachemap[(alg,'CT2A',ensC)]; rg=cachemap[(alg,'GL261',ensG)]
   # A1 inventory excludes HEX1 only; A2's dual-anchor target inventory excludes both anchors.
   ids=[r for r in rxns if r!='HEX1' and (arm=='A1' or r!='LDH_L')]
   ix=np.asarray([rindex[x] for x in ids])
   vals=descriptor_block(np.take(data['mats'][lc],ix,axis=1),np.take(data['mats'][rg],ix,axis=1),wc,wg)
   # Store the common 4,179-reaction target rows for paired sensitivity and use all 4,180 for A1 parity separately.
   for k in ('eta2','H_dir','D_dir'):
    arr=np.asarray(vals[k],float)
    if arm=='A1': arr=arr[[i for i,r in enumerate(ids) if r!='LDH_L']]
    baseline[arm][k][eidx,:]=arr
 parity=[]
 for arm in ('A1','A2'):
  s=summaries[arm].set_index('reaction_id')
  for metric,col in [('eta2','sign_magnitude_eta2'),('H_dir','directional_entropy3'),('D_dir','dominant_sign_mass')]:
   finite_count=np.isfinite(baseline[arm][metric]).sum(axis=0)
   got=np.asarray([np.median(v[np.isfinite(v)]) if np.isfinite(v).any() else np.nan for v in baseline[arm][metric].T],float)
   for i,r in enumerate(target):
    expn=int(s.loc[r,col+'_finite_count']); gotn=int(finite_count[i]); expval=s.loc[r,col+'_median']
    exp=float(expval) if pd.notna(expval) else float('nan'); obs=float(got[i]) if gotn else float('nan')
    if gotn!=expn: raise RuntimeError(f'{arm} {metric} finite-count parity failed for {r}: {gotn} != {expn}')
    if gotn==0:
     if math.isfinite(exp): raise RuntimeError(f'{arm} {metric} NaN parity failed for {r}')
     diff=0.0
    else:
     if not math.isfinite(exp): raise RuntimeError(f'{arm} {metric} frozen value is nonfinite for {r}')
     diff=abs(obs-exp)
    parity.append({'arm':arm,'metric':metric,'reaction_id':r,'recomputed_median':obs,'frozen_median':exp,'finite_count':gotn,'abs_difference':diff})
 maxdiff=max(r['abs_difference'] for r in parity)
 if not math.isfinite(maxdiff) or maxdiff>1e-10: raise RuntimeError(f'20x20 geometry parity failed ({maxdiff}); stopping Experiment 3A')
 # Shared candidate identities; matched support is checked under both arms before any descriptor/outcome aggregation.
 pairs=data['pairs'].loc[data['pairs'].pair_evaluable.eq('true')].drop_duplicates(['evaluation_id','ct2a_candidate_id','gl261_candidate_id']).copy()
 if len(pairs)!=836: raise RuntimeError('original distinct evaluable truth-pair denominator changed')
 eid=evals.set_index('evaluation_id').to_dict('index'); truth_ids={}
 for row in panel.to_dict('records'): truth_ids[pl2a.candidate_identity(row)]=row
 pairrows=[]
 for p in pairs.itertuples(index=False):
  ev=eid[str(p.evaluation_id)]; key=(str(p.evaluation_id),str(p.ct2a_candidate_id),str(p.gl261_candidate_id))
  valid=True
  for arm,wt in [('A1',data['a1bio']),('A2',data['a2'])]:
   for mouse,tumor,cid in [(str(p.ct2a_mouse),'CT2A',str(p.ct2a_candidate_id)),(str(p.gl261_mouse),'GL261',str(p.gl261_candidate_id))]:
    c=truth_ids.get(cid)
    if c is None: raise RuntimeError('frozen candidate identity missing')
    w=original_weights(wt,mouse,tumor,str(ev['algorithm']),str(ev['rna_context_key']))
    if w is None or w.sum()-w[int(c['sample_index'])]<=0: valid=False
  if valid: pairrows.append(key)
 common=set(pairrows)
 def wf(mouse,tumor,alg,rna,target,arm): return original_weights(data['a1bio'] if arm=='A1' else data['a2'],mouse,tumor,alg,rna)
 # Four-method deterministic pilot: one common truth pair and one weak-cue reaction per method.
 pilot_rows=[]; eval_map=evals.set_index('evaluation_id').to_dict('index'); rx=target[0]
 for alg in ('CORDA','GIMME','RIPTiDe','iMAT'):
  candidates=pairs.loc[pairs.evaluation_id.map(lambda x:str(eval_map[str(x)]['algorithm'])==alg)]
  chosen=None
  for p in candidates.itertuples(index=False):
   key=(str(p.evaluation_id),str(p.ct2a_candidate_id),str(p.gl261_candidate_id))
   if key in common: chosen=p; break
  if chosen is None: raise RuntimeError(f'3A pilot has no matched truth pair for {alg}')
  ev=eval_map[str(chosen.evaluation_id)]; rowC=truth_ids[str(chosen.ct2a_candidate_id)]; rowG=truth_ids[str(chosen.gl261_candidate_id)]
  li=cachemap[(alg,'CT2A',str(rowC['ensemble_hash']))]; ri=cachemap[(alg,'GL261',str(rowG['ensemble_hash']))]; ix=rindex[rx]
  for arm,wt in [('A1',data['a1bio']),('A2',data['a2'])]:
   wc=wf(str(chosen.ct2a_mouse),'CT2A',alg,str(ev['rna_context_key']),20.,arm); wg=wf(str(chosen.gl261_mouse),'GL261',alg,str(ev['rna_context_key']),20.,arm)
   sc=int(chosen.ct2a_sample_index); sg=int(chosen.gl261_sample_index); kc=np.arange(20)!=sc; kg=np.arange(20)!=sg
   wc[sc]=0.; wg[sg]=0.; wc/=wc.sum(); wg/=wg.sum()
   left=np.take(data['mats'][li],[ix],axis=1)[kc]; right=np.take(data['mats'][ri],[ix],axis=1)[kg]
   z=descriptor_block(left,right,wc[kc],wg[kg])
   expected=pl1.target_landscape_metrics(anchor_ct2a=data['mats'][li,kc,rxns.index('HEX1')],anchor_gl261=data['mats'][ri,kg,rxns.index('HEX1')],target_ct2a=left[:,0],target_gl261=right[:,0],strong_weights_ct2a=wc[kc],strong_weights_gl261=wg[kg])
   expected_values={'eta2':expected['sign_magnitude_eta2'],'H_dir':expected['directional_entropy3'],'D_dir':expected['dominant_sign_mass']}
   diffs={}; undefined=[]
   for metric,ex in expected_values.items():
    ob=float(z[metric][0]); ex=float(ex)
    if math.isfinite(ob) and math.isfinite(ex): diffs[metric]=abs(ob-ex)
    elif math.isnan(ob) and math.isnan(ex): diffs[metric]=0.0; undefined.append(metric)
    else: raise RuntimeError(f'3A pilot defined/undefined parity mismatch for {alg}/{arm}/{metric}')
   md=max(diffs.values())
   if md>1e-12 or len(wc[kc])!=19 or len(wg[kg])!=19 or abs(wc[kc].sum()-1)>1e-12 or abs(wg[kg].sum()-1)>1e-12: raise RuntimeError('3A deterministic 19x19 pilot failed')
   pilot_rows.append({'algorithm':alg,'arm':arm,'evaluation_id':str(chosen.evaluation_id),'ct2a_candidate_id':str(chosen.ct2a_candidate_id),'gl261_candidate_id':str(chosen.gl261_candidate_id),'reaction_id':rx,'left_support':19,'right_support':19,'undefined_descriptor_fields':','.join(undefined),'max_descriptor_abs_diff':md,'status':'PASS'})
 out.mkdir(parents=True,exist_ok=True)
 (out/'experiment3a_pilot.json').write_text(json.dumps({'status':'PASS','selected_evaluations':4,'method_panel':['CORDA','GIMME','RIPTiDe','iMAT'],'reaction_id':rx,'cases':pilot_rows},indent=2)+'\n')
 if pilot:
  pd.DataFrame(parity).to_csv(out/'experiment3a_20x20_parity.tsv.gz',sep='\t',index=False,compression='gzip')
  return
 # Use the same distinct pair set for both arms; utility and descriptors are reaction-level aggregates.
 outcome=[]; supports=[]; geo=defaultdict(lambda:{k:0. for k in ('eta2','H_dir','D_dir','K_dir','non_tie_coverage','eta2_20x20','H_dir_20x20','D_dir_20x20','n','n_eta2','n_H_dir','n_D_dir','n_eta2_20x20','n_H_dir_20x20','n_D_dir_20x20')})
 p2=pairs.copy(); p2['pair_evaluable']='false'
 for idx,p in p2.iterrows():
  k=(str(p.evaluation_id),str(p.ct2a_candidate_id),str(p.gl261_candidate_id)); p2.loc[idx,'pair_evaluable']='true' if k in common else 'false'
 data['pairs']=p2
 for arm in ('A1','A2'):
  rows,supp=score_frozen_truth_pairs(data,wf,common,scenario='truth_excluded_q10_q50_q90',arm=arm,target=20.,outdir=out)
  outcome.extend(rows); supports.extend(supp)
  for p in pairs.loc[pairs.apply(lambda x:(str(x.evaluation_id),str(x.ct2a_candidate_id),str(x.gl261_candidate_id)) in common,axis=1)].itertuples(index=False):
   ev=eid[str(p.evaluation_id)]; alg=str(ev['algorithm']); rna=str(ev['rna_context_key'])
   cl=truth_ids[str(p.ct2a_candidate_id)]; gr=truth_ids[str(p.gl261_candidate_id)]
   ensC=str(cl['ensemble_hash']); ensG=str(gr['ensemble_hash']); li=cachemap[(alg,'CT2A',ensC)]; ri=cachemap[(alg,'GL261',ensG)]
   wc_full=wf(str(p.ct2a_mouse),'CT2A',alg,rna,20.,arm); wg_full=wf(str(p.gl261_mouse),'GL261',alg,rna,20.,arm)
   wc=wc_full.copy(); wg=wg_full.copy()
   def hold(w,sample):
    x=w.copy(); x[int(sample)]=0.; return x/x.sum()
   wc=hold(wc,int(p.ct2a_sample_index)); wg=hold(wg,int(p.gl261_sample_index))
   ix=np.asarray([rindex[x] for x in target])
   keepC=np.arange(20)!=int(p.ct2a_sample_index); keepG=np.arange(20)!=int(p.gl261_sample_index)
   z=descriptor_block(np.take(data['mats'][li],ix,axis=1)[keepC],np.take(data['mats'][ri],ix,axis=1)[keepG],wc[keepC],wg[keepG])
   z20=descriptor_block(np.take(data['mats'][li],ix,axis=1),np.take(data['mats'][ri],ix,axis=1),wc_full,wg_full)
   truths=data['mats'][li,int(p.ct2a_sample_index),ix]-data['mats'][ri,int(p.gl261_sample_index),ix]
   for j,r in enumerate(target):
    if abs(truths[j])<=ATOL: continue
    a=geo[(arm,r)]
    for k in ('eta2','H_dir','D_dir','K_dir','non_tie_coverage'):
     if math.isfinite(float(z[k][j])):
      a[k]+=float(z[k][j])
      if k in ('eta2','H_dir','D_dir'): a[f'n_{k}']+=1
    for k in ('eta2','H_dir','D_dir'):
     if math.isfinite(float(z20[k][j])):
      a[f'{k}_20x20']+=float(z20[k][j]); a[f'n_{k}_20x20']+=1
    a['n']+=1
  # Produce per-arm reaction-level tables after equal distinct-pair aggregation.
 georows=[]
 for (arm,r),z in geo.items():
  n=int(z['n']); georows.append({'arm':arm,'reaction_id':r,'n_distinct_non_tie_truth_pairs':n,**{f'n_{k}_truth_pairs':int(z[f'n_{k}']) for k in ('eta2','H_dir','D_dir')},**{f'n_{k}_20x20_truth_pairs':int(z[f'n_{k}_20x20']) for k in ('eta2','H_dir','D_dir')},**{k:(z[k]/z[f'n_{k}'] if z[f'n_{k}'] else float('nan')) for k in ('eta2','H_dir','D_dir')},**{k:(z[k]/n if n else float('nan')) for k in ('K_dir','non_tie_coverage')},**{f'{k}_20x20':(z[f'{k}_20x20']/z[f'n_{k}_20x20'] if z[f'n_{k}_20x20'] else float('nan')) for k in ('eta2','H_dir','D_dir')}})
 out.mkdir(parents=True,exist_ok=True)
 pd.DataFrame(parity).to_csv(out/'experiment3a_20x20_parity.tsv.gz',sep='\t',index=False,compression='gzip')
 pd.DataFrame(georows).to_csv(out/'experiment3a_truth_excluded_geometry.tsv.gz',sep='\t',index=False,compression='gzip')
 geometry=pd.DataFrame(georows)
 for metric in ('eta2','H_dir','D_dir'):
  geometry[f'{metric}_19x19']=geometry[metric]
  geometry[f'{metric}_change_19x19_minus_20x20']=geometry[f'{metric}_19x19']-geometry[f'{metric}_20x20']
 geometry.to_csv(out/'experiment3a_geometry_comparison.tsv.gz',sep='\t',index=False,compression='gzip')
 pd.DataFrame(outcome).to_csv(out/'experiment3a_usefulness.tsv.gz',sep='\t',index=False,compression='gzip')
 pd.DataFrame(supports).to_csv(out/'experiment3a_matched_pair_support.tsv.gz',sep='\t',index=False,compression='gzip')
 geomdf=pd.DataFrame(georows); outdf=pd.DataFrame(outcome)
 aligned=outdf.merge(geomdf,on=['arm','reaction_id'],how='outer',validate='one_to_one',suffixes=('_outcome','_geometry'))
 if aligned[['non_tie_cases','n_distinct_non_tie_truth_pairs']].isna().any().any() or not np.array_equal(aligned.non_tie_cases.to_numpy(int),aligned.n_distinct_non_tie_truth_pairs.to_numpy(int)): raise RuntimeError('truth-pair geometry/outcome aggregation denominator mismatch')
 rxreg=pd.read_csv(root/'outputs/dmi_bridge_pl2a_sign_only_prepare_v1/BRIDGEPL2A_REACTION_REGISTRY.tsv',sep='\t')
 split=pd.read_csv(root/'outputs/dmi_bridge_pl2c_geometry_utility_development_v1/BRIDGEPL2C_REACTION_SPLIT.tsv',sep='\t')
 assoc=[]; common_eta2_counts={}
 for arm in ('A1','A2'):
  armgeom=geomdf.loc[geomdf.arm.eq(arm)].copy()
  common_reactions=set(armgeom.loc[armgeom.eta2.notna()&armgeom.eta2_20x20.notna(),'reaction_id'])
  common_eta2_counts[arm]=len(common_reactions)
  sharedout=outdf.loc[outdf.arm.eq(arm)&outdf.reaction_id.isin(common_reactions)].copy()
  if not len(common_reactions): raise RuntimeError(f'no common finite eta2 reaction population for {arm}')
  for spec,cols in [('original_20x20',{'eta2':'eta2_20x20','H_dir':'H_dir_20x20','D_dir':'D_dir_20x20'}),('truth_excluded_19x19',{'eta2':'eta2','H_dir':'H_dir','D_dir':'D_dir'})]:
   g=armgeom[['arm','reaction_id',*cols.values()]].rename(columns={v:k for k,v in cols.items()})
   sub=associations(sharedout,g,g,rxreg,split); sub=sub.loc[sub.arm.eq(arm)].copy(); sub.insert(1,'geometry_spec',spec); assoc.append(sub)
 assoc=pd.concat(assoc,ignore_index=True)
 assoc.to_csv(out/'experiment3a_eta2_usefulness.tsv',sep='\t',index=False)
 (out/'experiment3a_manifest.json').write_text(json.dumps({'status':'20x20_PARITY_PASS','geometry_label':'truth-exclusion sensitivity descriptor; not a pre-cue predictor','full_descriptor_reaction_inventory':4180,'weak_cue_reaction_inventory':4179,'frozen_evaluable_distinct_pairs':836,'matched_a1_a2_pairs':len(common),'common_finite_eta2_reactions_by_arm':common_eta2_counts,'weight_calibration':'frozen original A1 and A2 ESS20 weights; no recalibration','max_parity_abs_diff':maxdiff,'source_hashes':EXPECTED|{A1_CURRENT:A1_CURRENT_SHA},'source_root':str(root),'code_sha256':{p:sha256(root/p) for p in CODE_PATHS},'archived_weight_kernel_sha256':ARCHIVED_KERNEL_SHA},indent=2)+'\n')

def score_frozen_truth_pairs(data, weights_for, common_keys, *, scenario:str, arm:str, target:float, outdir:Path, truth_pairs=None):
 """Re-evaluate frozen distinct truth pairs using PL2's exact sign-only endpoint."""
 cache=data['cache']; evals=data['evals']; pairs=data['pairs']; rxns=data['rxns']
 ids={r:i for i,r in enumerate(rxns)}
 # PL2 uses 4,179 weak-cue reactions, excluding both anchor reactions.
 targets=[r for r in rxns if r not in {'HEX1','LDH_L'}]
 if len(targets)!=4179: raise RuntimeError('weak-cue reaction universe must be 4,179')
 ridx=np.asarray([ids[r] for r in targets],dtype=int)
 cache_map=candidate_index(data)
 source_pairs=(pairs if truth_pairs is None else truth_pairs)
 truth=source_pairs.loc[source_pairs.pair_evaluable.astype(str).eq('true')].drop_duplicates(['evaluation_id','ct2a_candidate_id','gl261_candidate_id']).copy()
 if len(truth)==0: raise RuntimeError('scenario has no distinct evaluable truth pairs')
 if scenario=='original_q10_q50_q90' and len(truth)!=836: raise RuntimeError('original distinct evaluable truth-pair population must be exactly 836')
 e_map=evals.set_index('evaluation_id').to_dict('index')
 total={r:{'n':0,'non_tie':0,'use':0.,'cg':0.,'wg':0.,'ig':0.,'sep':0.} for r in targets}
 support=[]
 for pr in truth.itertuples(index=False):
  ev=e_map[str(pr.evaluation_id)]
  casekey=(str(pr.evaluation_id),str(pr.ct2a_candidate_id),str(pr.gl261_candidate_id))
  if casekey not in common_keys:
   support.append({'scenario':scenario,'arm':arm,'target_ess':target,'evaluation_id':casekey[0],'ct2a_candidate_id':casekey[1],'gl261_candidate_id':casekey[2],'eligible':False,'reason':'outside paired A1/A2 support intersection'})
   continue
  conditions=[]
  valid=True
  for tumor,mouse,cid,sample in [('CT2A',str(pr.ct2a_mouse),str(pr.ct2a_candidate_id),int(pr.ct2a_sample_index)),('GL261',str(pr.gl261_mouse),str(pr.gl261_candidate_id),int(pr.gl261_sample_index))]:
   key=(mouse,tumor,str(ev['algorithm']),str(ev['ensemble_hash']) if 'ensemble_hash' in ev else '',str(ev['rna_context_key']))
   # Candidate ensemble identity is recovered from the evaluation's method/context and tumor.
   cands=data['panel'].loc[data['panel'].algorithm.astype(str).eq(str(ev['algorithm'])) & data['panel'].tumor.astype(str).eq(tumor) & data['panel'].rna_context_key.astype(str).eq(str(ev['rna_context_key']))].copy()
   if len(cands)!=20: raise RuntimeError('evaluation-to-condition candidate pool mismatch')
   ensemble=str(cands.ensemble_hash.iloc[0]); cache_i=cache_map[(str(ev['algorithm']),tumor,ensemble)]
   sample_to_pos={int(x):i for i,x in enumerate(range(20))}
   samples=np.arange(20,dtype=int)
   w=weights_for(mouse,tumor,str(ev['algorithm']),str(ev['rna_context_key']),target,arm)
   if w is None or w.sum()<=0 or sample not in sample_to_pos or w[sample_to_pos[sample]]<0:
    valid=False; break
   w=w/w.sum(); held=w.copy(); held[sample_to_pos[sample]]=0.; mass=float(held.sum())
   if mass<=0: valid=False; break
   held/=mass
   kept=np.delete(samples,sample_to_pos[sample])
   conditions.append((cache_i,sample,kept,held[kept]))
  if not valid:
   support.append({'scenario':scenario,'arm':arm,'target_ess':target,'evaluation_id':str(pr.evaluation_id),'ct2a_candidate_id':casekey[1],'gl261_candidate_id':casekey[2],'eligible':False,'reason':'zero conditional support after truth exclusion'})
   continue
  (li,ls,lkeep,lw),(ri,rs,rkeep,rw)=conditions
  left=data['mats'][li,lkeep,:][:,ridx]
  right=data['mats'][ri,rkeep,:][:,ridx]
  truth_delta=data['mats'][li,ls,ridx]-data['mats'][ri,rs,ridx]
  # Cross-check a deterministic reaction against the frozen scalar operator before accepting this case.
  vec=vector_outcomes(left[:,[0]],right[:,[0]],lw,rw,truth_delta[[0]])
  delta=(left[:,0,None]-right[None,:,0]).reshape(-1); joint=(lw[:,None]*rw[None,:]).reshape(-1)
  scalar=scalar_outcome(delta,joint,float(truth_delta[0]),data['pl2'])
  if int(vec['truth_direction'][0])!=int(scalar['truth_direction']): raise RuntimeError('vector/scalar truth-direction parity failed')
  lam=scalar.get('lambda_results',[])
  if lam and abs(float(vec['correct_gain'][0])-float(lam[0]['correct_abs_error_gain']))>1e-10: raise RuntimeError('vector/scalar gain parity failed')
  values=vector_outcomes(left,right,lw,rw,truth_delta)
  non=np.abs(truth_delta)>ATOL
  for j,r in enumerate(targets):
   if not non[j]: continue
   z=total[r]; z['n']+=1; z['non_tie']+=1; z['use']+=float(values['usefulness'][j]); z['cg']+=float(values['correct_gain'][j]); z['wg']+=float(values['wrong_gain'][j]); z['ig']+=float(values['information_gain'][j]); z['sep']+=float(values['correct_gain'][j]-values['wrong_gain'][j])
  support.append({'scenario':scenario,'arm':arm,'target_ess':target,'evaluation_id':str(pr.evaluation_id),'ct2a_candidate_id':casekey[1],'gl261_candidate_id':casekey[2],'eligible':True,'reason':''})
 rows=[]
 for r,z in total.items():
  if z['non_tie']:
   rows.append({'scenario':scenario,'arm':arm,'target_ess':target,'reaction_id':r,'n_distinct_truth_cases':z['n'],'non_tie_cases':z['non_tie'],'usefulness':z['use']/z['non_tie'],'correct_gain_mean':z['cg']/z['non_tie'],'wrong_gain_mean':z['wg']/z['non_tie'],'information_gain_mean':z['ig']/z['non_tie'],'correct_minus_wrong_mean':z['sep']/z['non_tie']})
 return rows,support

def build_alternative_truths(data):
 """Outcome-blind weighted HEX1 quantiles and seeded weighted whole-vector draws."""
 import dmi_bridge_pl2a_prepare_v1 as pl2a
 panel=data['panel']; cache=data['cache']; rxns=data['rxns']; hidx=rxns.index('HEX1')
 original=data['truth'].copy()
 original['scenario']='original_q10_q50_q90'; original['slot']=original['requested_quantile'].map({'0.10':'q10','0.50':'q50','0.90':'q90'}); original['seed']=''
 original=original.rename(columns={'requested_quantile':'quantile_original'})
 cache_map=candidate_index(data); wtab=data['a1bio']; records=[]
 pools=panel.groupby(['algorithm','tumor','ensemble_hash','rna_context_key'],sort=True)
 mice=sorted(data['a1bio'][['mouse_id','tumor']].drop_duplicates().itertuples(index=False,name=None))
 for mouse,tumor in mice:
  for (alg,ptumor,ensemble,rna),g in pools:
   if ptumor!=tumor: continue
   if len(g)!=20: raise RuntimeError('truth selection candidate pool is not 20')
   cache_i=cache_map[(str(alg),str(tumor),str(ensemble))]
   rows=g.to_dict('records')
   rows.sort(key=lambda r:(float(data['mats'][cache_i,int(r['sample_index']),hidx]),pl2a.candidate_identity(r)))
   w=original_weights(wtab,str(mouse),str(tumor),str(alg),str(rna))
   if w is None: raise RuntimeError('original A1 weighting support missing for truth selection')
   bysample={int(r.sample_index):float(r.weight) for r in wtab.loc[wtab.mouse_id.astype(str).eq(str(mouse)) & wtab.tumor.astype(str).eq(str(tumor)) & wtab.algorithm.astype(str).eq(str(alg)) & wtab.rna_context_key.astype(str).eq(str(rna))].itertuples()}
   weights=np.asarray([bysample[int(r['sample_index'])] for r in rows],float); weights/=weights.sum(); cum=np.cumsum(weights)
   for q in (.20,.50,.80):
    j=min(int(np.searchsorted(cum,q,side='left')),19); r=rows[j]
    records.append({'scenario':'weighted_q20_q50_q80','slot':f'q{int(q*100)}','seed':'','mouse_id':str(mouse),'tumor':str(tumor),'algorithm':str(alg),'rna_context_key':str(rna),'candidate_id':pl2a.candidate_identity(r),'sample_index':int(r['sample_index']),'HEX1':float(data['mats'][cache_i,int(r['sample_index']),hidx])})
 for seed in TRUTH_SEEDS:
  rng=np.random.Generator(np.random.PCG64(seed))
  for mouse,tumor in mice:
   for (alg,ptumor,ensemble,rna),g in pools:
    if ptumor!=tumor: continue
    cache_i=cache_map[(str(alg),str(tumor),str(ensemble))]
    rows=g.sort_values('sample_index',kind='stable').to_dict('records')
    w=original_weights(wtab,str(mouse),str(tumor),str(alg),str(rna))
    if w is None: raise RuntimeError('random truth selection support missing')
    for draw,sample in enumerate(rng.choice(20,size=3,replace=True,p=w/w.sum()),start=1):
     r=rows[int(sample)]
     records.append({'scenario':f'weighted_random_seed_{seed}','slot':f'draw{draw}','seed':seed,'mouse_id':str(mouse),'tumor':str(tumor),'algorithm':str(alg),'rna_context_key':str(rna),'candidate_id':pl2a.candidate_identity(r),'sample_index':int(r['sample_index']),'HEX1':float(data['mats'][cache_i,int(r['sample_index']),hidx])})
 alt=pd.DataFrame(records)
 # The frozen primary truth set remains in the same table as the alternative scenarios.
 cols=['scenario','slot','seed','mouse_id','tumor','algorithm','rna_context_key','candidate_id','sample_index','HEX1']
 return pd.concat([original[cols],alt[cols]],ignore_index=True)

def pair_alternative_truths(data,selections,weight_fn):
 cand={}
 for r in data['panel'].to_dict('records'):
  import dmi_bridge_pl2a_prepare_v1 as pl2a
  cand[pl2a.candidate_identity(r)]=r
 evals=data['evals']; pairs=[]; duplicate_rows=[]
 for scenario,sel in selections.groupby('scenario',sort=True):
  slot_values=sorted(sel.slot.unique())
  bypool={(str(r.mouse_id),str(r.tumor),str(r.algorithm),str(r.rna_context_key),str(r.slot)):r for r in sel.itertuples(index=False)}
  for ev in evals.itertuples(index=False):
   for slot in slot_values:
    l= bypool.get((str(ev.ct2a_mouse),'CT2A',str(ev.algorithm),str(ev.rna_context_key),str(slot)))
    r= bypool.get((str(ev.gl261_mouse),'GL261',str(ev.algorithm),str(ev.rna_context_key),str(slot)))
    if l is None or r is None: continue
    lc=cand[str(l.candidate_id)]; rc=cand[str(r.candidate_id)]
    common=True
    for arm in ('A1','A2'):
     for mouse,tumor,context,row in [(str(ev.ct2a_mouse),'CT2A',str(ev.rna_context_key),lc),(str(ev.gl261_mouse),'GL261',str(ev.rna_context_key),rc)]:
      w=weight_fn(mouse,tumor,str(ev.algorithm),context,20.,arm)
      sample=int(row['sample_index'])
      if w is None or w.sum()-w[sample]<=0: common=False
    pairs.append({'scenario':scenario,'slot':str(slot),'evaluation_id':str(ev.evaluation_id),'algorithm':str(ev.algorithm),'rna_context_key':str(ev.rna_context_key),'ct2a_mouse':str(ev.ct2a_mouse),'gl261_mouse':str(ev.gl261_mouse),'ct2a_candidate_id':str(l.candidate_id),'gl261_candidate_id':str(r.candidate_id),'ct2a_sample_index':int(l.sample_index),'gl261_sample_index':int(r.sample_index),'pair_evaluable':'true' if common else 'false'})
  sc=sel.copy()
  dup=sc.groupby(['mouse_id','tumor','algorithm','rna_context_key']).candidate_id.agg(lambda x:len(x)-x.nunique()).reset_index(name='duplicate_selection_count')
  dup.insert(0,'scenario',scenario); duplicate_rows.extend(dup.to_dict('records'))
 return pd.DataFrame(pairs),pd.DataFrame(duplicate_rows)

def associations(outcomes:pd.DataFrame, geom_a1:pd.DataFrame, geom_a2:pd.DataFrame, reaction_registry:pd.DataFrame, split:pd.DataFrame):
 from scipy.stats import spearmanr,rankdata
 from sklearn.model_selection import GroupKFold
 from sklearn.pipeline import make_pipeline
 from sklearn.preprocessing import StandardScaler
 from sklearn.linear_model import Ridge
 results=[]
 for arm,geom in [('A1',geom_a1),('A2',geom_a2)]:
  gm=geom.loc[geom.arm.eq(arm)].drop(columns=['arm'],errors='ignore') if 'arm' in geom.columns else geom
  dat=outcomes.loc[outcomes.arm.eq(arm)].merge(gm,on='reaction_id',validate='one_to_one').merge(reaction_registry[['reaction_id','subsystem']],on='reaction_id',validate='one_to_one').merge(split[['reaction_id','analysis_split']],on='reaction_id',validate='one_to_one')
  for split_name,sub in [('ALL',dat),*[(x,dat.loc[dat.analysis_split.eq(x)]) for x in sorted(dat.analysis_split.unique())]]:
   for subset in ('eta2','H_dir'):
    z=sub[[subset,'usefulness']].dropna(); rho=float(spearmanr(z[subset],z.usefulness).statistic) if len(z)>2 else float('nan')
    results.append({'arm':arm,'reaction_split':split_name,'descriptor':subset,'n_reactions':len(z),'spearman_rho':rho})
   z=sub[['eta2','H_dir','usefulness','subsystem']].dropna()
   if len(z)>4 and z.subsystem.nunique()>=5:
    xr=rankdata(z.eta2); yr=rankdata(z.usefulness); hr=rankdata(z.H_dir); X=np.c_[np.ones(len(z)),hr]; residx=xr-X@np.linalg.lstsq(X,xr,rcond=None)[0]; residy=yr-X@np.linalg.lstsq(X,yr,rcond=None)[0]; partial=float(np.corrcoef(residx,residy)[0,1])
    folds=list(GroupKFold(n_splits=5).split(z[['H_dir','eta2']],z.usefulness,z.subsystem.astype(str)))
    score={}
    for name,cols in [('entropy', ['H_dir']),('entropy_plus_eta2',['H_dir','eta2'])]:
     pred=np.full(len(z),np.nan); x=z[cols].to_numpy(float); y=z.usefulness.to_numpy(float)
     for tr,te in folds:
      model=make_pipeline(StandardScaler(),Ridge(alpha=1.0)); model.fit(x[tr],y[tr]); pred[te]=model.predict(x[te])
     score[name]=float(np.mean(np.abs(pred-y)))
    results.append({'arm':arm,'reaction_split':split_name,'descriptor':'eta2_increment','n_reactions':len(z),'partial_rank_eta2_given_entropy':partial,'mae_entropy':score['entropy'],'mae_entropy_plus_eta2':score['entropy_plus_eta2'],'mae_improvement':score['entropy']-score['entropy_plus_eta2']})
 return pd.DataFrame(results)

def run3b(data,out):
 root=data['root']; out.mkdir(parents=True,exist_ok=True)
 selections=build_alternative_truths(data)
 # Calibration-20 conditional weights are the exact frozen A1/A2 vectors.
 def wf(mouse,tumor,alg,rna,target,arm): return original_weights(data['a1bio'] if arm=='A1' else data['a2'],mouse,tumor,alg,rna)
 pairs,duplicates=pair_alternative_truths(data,selections,wf)
 allout=[]; alls=[]; saved=data['pairs'].copy()
 for scenario,pp in pairs.groupby('scenario',sort=True):
  # The shared support intersection is fixed before computing either arm.
  common=set((str(r.evaluation_id),str(r.ct2a_candidate_id),str(r.gl261_candidate_id)) for r in pp.loc[pp.pair_evaluable.eq('true')].itertuples())
  pp=pp.copy(); pp['pair_evaluable']=pp.apply(lambda r:'true' if (str(r.evaluation_id),str(r.ct2a_candidate_id),str(r.gl261_candidate_id)) in common else 'false',axis=1)
  data['pairs']=pp
  for arm in ('A1','A2'):
   rows,supp=score_frozen_truth_pairs(data,wf,common,scenario=scenario,arm=arm,target=20.,outdir=out,truth_pairs=pp)
   allout.extend(rows); alls.extend(supp)
 data['pairs']=saved
 # Frozen 20x20 geometry is the pre-cue descriptor for alternative truth populations.
 geom=[]
 for arm,path in [('A1',root/'outputs/dmi_bridge_pl1_predictability_landscape_v1_resumable_v1/BRIDGEPL1_REACTION_SUMMARY.tsv'),('A2',root/'outputs/dmi_bridge_a21_dual_anchor_geometry_v1/BRIDGEA21_REACTION_SUMMARY.tsv.xz')]:
  f=pd.read_csv(path,sep='\t',compression='infer')
  geom.append((arm,f[['reaction_id','sign_magnitude_eta2_median','directional_entropy3_median','dominant_sign_mass_median']].rename(columns={'sign_magnitude_eta2_median':'eta2','directional_entropy3_median':'H_dir','dominant_sign_mass_median':'D_dir'})))
 rxreg=pd.read_csv(root/'outputs/dmi_bridge_pl2a_sign_only_prepare_v1/BRIDGEPL2A_REACTION_REGISTRY.tsv',sep='\t')
 split=pd.read_csv(root/'outputs/dmi_bridge_pl2c_geometry_utility_development_v1/BRIDGEPL2C_REACTION_SPLIT.tsv',sep='\t')
 outdf=pd.DataFrame(allout); outdf['scenario']=outdf.scenario.astype(str)
 metrics=[]
 for (scenario,arm),sub in outdf.groupby(['scenario','arm'],sort=True):
  for scope,ids in [('ALL',set(sub.reaction_id)),*[(x,set(split.loc[split.analysis_split.eq(x),'reaction_id'])) for x in sorted(split.analysis_split.unique())]]:
   p=sub.loc[sub.reaction_id.isin(ids)]
   for garm,g in geom:
    if garm!=arm: continue
    m=p.merge(g,on='reaction_id',validate='one_to_one')
    assoc=associations(m.assign(arm=arm),g, g,rxreg,split) if False else None
    for desc in ('eta2','H_dir'):
     from scipy.stats import spearmanr
     z=m[[desc,'usefulness']].dropna(); metrics.append({'scenario':scenario,'arm':arm,'reaction_split':scope,'descriptor':desc,'n_reactions':len(z),'spearman_rho':float(spearmanr(z[desc],z.usefulness).statistic) if len(z)>2 else float('nan')})
    from scipy.stats import rankdata
    from sklearn.model_selection import GroupKFold
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import Ridge
    z=m[['eta2','H_dir','usefulness','reaction_id']].merge(rxreg[['reaction_id','subsystem']],on='reaction_id',validate='one_to_one').dropna()
    if len(z)>4 and z.subsystem.nunique()>=5:
     xr=rankdata(z.eta2); yr=rankdata(z.usefulness); hr=rankdata(z.H_dir); X0=np.c_[np.ones(len(z)),hr]; rx=xr-X0@np.linalg.lstsq(X0,xr,rcond=None)[0]; ry=yr-X0@np.linalg.lstsq(X0,yr,rcond=None)[0]; partial=float(np.corrcoef(rx,ry)[0,1])
     folds=list(GroupKFold(n_splits=5).split(z[['H_dir','eta2']],z.usefulness,z.subsystem.astype(str))); maes=[]
     for cols in (['H_dir'],['H_dir','eta2']):
      pred=np.full(len(z),np.nan); xx=z[cols].to_numpy(float); yy=z.usefulness.to_numpy(float)
      for tr,te in folds:
       model=make_pipeline(StandardScaler(),Ridge(alpha=1.0)); model.fit(xx[tr],yy[tr]); pred[te]=model.predict(xx[te])
      maes.append(float(np.mean(np.abs(pred-yy))))
     metrics.append({'scenario':scenario,'arm':arm,'reaction_split':scope,'descriptor':'eta2_increment','n_reactions':len(z),'partial_rank_eta2_given_entropy':partial,'mae_entropy':maes[0],'mae_entropy_plus_eta2':maes[1],'mae_improvement':maes[0]-maes[1]})
 outdf.to_csv(out/'experiment3b_usefulness.tsv.gz',sep='\t',index=False,compression='gzip')
 pd.DataFrame(alls).to_csv(out/'experiment3b_pair_support.tsv.gz',sep='\t',index=False,compression='gzip')
 pairs.to_csv(out/'experiment3b_pair_eligibility.tsv.gz',sep='\t',index=False,compression='gzip')
 selections.to_csv(out/'experiment3b_selected_candidates.tsv.gz',sep='\t',index=False,compression='gzip')
 duplicates.to_csv(out/'experiment3b_duplicates.tsv',sep='\t',index=False)
 pd.DataFrame(metrics).to_csv(out/'experiment3b_associations.tsv',sep='\t',index=False)
 scenario_denominators={str(s):{'selection_rows':int(len(g)),'evaluable_selection_rows':int(g.pair_evaluable.eq('true').sum()),'distinct_evaluable_truth_pairs':int(g.loc[g.pair_evaluable.eq('true')].drop_duplicates(['evaluation_id','ct2a_candidate_id','gl261_candidate_id']).shape[0])} for s,g in pairs.groupby('scenario',sort=True)}
 (out/'experiment3b_manifest.json').write_text(json.dumps({'status':'COMPLETED','quantiles':[.20,.50,.80],'random_seeds':list(TRUTH_SEEDS),'random_sampling':'weighted categorical with replacement; three selections per pool','paired_arms':'same selected candidate IDs; scenarios use the common support intersection','scenario_denominators':scenario_denominators,'lambda':LAM,'weight_calibration':'frozen original A1 and A2 ESS20 weights; no recalibration','source_hashes':EXPECTED|{A1_CURRENT:A1_CURRENT_SHA},'source_root':str(root),'code_sha256':{p:sha256(root/p) for p in CODE_PATHS},'archived_weight_kernel_sha256':ARCHIVED_KERNEL_SHA},indent=2)+'\n')

def weighted_geometry_summary(data,weights_for,target,arm):
 """Reaction-level pre-cue geometry from each calibrated 20x20 candidate-pair distribution."""
 cache=data['cache']; evals=data['evals']; rxns=data['rxns']; rindex={r:i for i,r in enumerate(rxns)}; cachemap=candidate_index(data)
 targets=[r for r in rxns if r not in {'HEX1','LDH_L'}]; sums={k:np.zeros(len(targets),float) for k in ('eta2','H_dir','D_dir')}; counts={k:np.zeros(len(targets),int) for k in sums}; eval_counts=0
 for ev in evals.itertuples(index=False):
  alg=str(ev.algorithm); rna=str(ev.rna_context_key)
  wc=weights_for(str(ev.ct2a_mouse),'CT2A',alg,rna,target,arm); wg=weights_for(str(ev.gl261_mouse),'GL261',alg,rna,target,arm)
  if wc is None or wg is None or wc.sum()<=0 or wg.sum()<=0: continue
  wc=wc/wc.sum(); wg=wg/wg.sum()
  ec=str(data['panel'].loc[data['panel'].algorithm.eq(alg)&data['panel'].tumor.eq('CT2A')&data['panel'].rna_context_key.eq(rna),'ensemble_hash'].iloc[0])
  eg=str(data['panel'].loc[data['panel'].algorithm.eq(alg)&data['panel'].tumor.eq('GL261')&data['panel'].rna_context_key.eq(rna),'ensemble_hash'].iloc[0])
  li=cachemap[(alg,'CT2A',ec)]; ri=cachemap[(alg,'GL261',eg)]
  ix=np.asarray([rindex[r] for r in targets])
  z=descriptor_block(np.take(data['mats'][li],ix,axis=1),np.take(data['mats'][ri],ix,axis=1),wc,wg)
  eval_counts+=1
  for k in sums:
   finite=np.isfinite(z[k]); sums[k][finite]+=z[k][finite]; counts[k][finite]+=1
 return pd.DataFrame({'reaction_id':targets,'arm':arm,'target_ess':target,'n_geometry_evaluations':eval_counts,**{f'n_{k}_finite':counts[k] for k in sums},**{k:np.divide(v,counts[k],out=np.full_like(v,np.nan),where=counts[k]>0) for k,v in sums.items()}})

def run4(data,out):
 """Recalibrate A1/A2 globally at 320 candidates, gated by exact ESS20 parity."""
 root=data['root']; panel=data['panel']; a20=data['a20mod']
 a1base=pd.read_csv(root/A1_CURRENT,sep='\t',compression='xz',float_precision='round_trip')
 a2base=data['a2'].copy()
 mouse=pd.read_csv(root/'outputs/dmi_bridge_a20_dual_anchor_qualification_v1/BRIDGEA20_MOUSE_COORDINATES.tsv',sep='\t')
 coords={(str(r.mouse_id),str(r.observable)):float(r.coordinate) for r in mouse.itertuples()}
 idc=['algorithm','tumor','ensemble_hash','projection_hash','rna_context_key','sample_index']
 # Merge raw panel observables into the authoritative baseline row order (320 per mouse).
 raw=panel[idc+['EX_glc__D_e','LDH_L']].copy()
 raw['Vmax_raw']=np.maximum(-raw.EX_glc__D_e.to_numpy(float),0.)
 raw['Vlac_raw']=np.maximum(-raw.LDH_L.to_numpy(float),0.)
 a1base=a1base.merge(raw,on=idc,how='left',validate='many_to_one',sort=False)
 if len(a1base)!=3200 or a1base[['Vmax_raw','Vlac_raw']].isna().any().any(): raise RuntimeError('raw observable join failed')
 outputs=[]; parity=[]; all_weights={}
 for mouse_id,group in a1base.groupby('mouse_id',sort=True):
  g=group.copy(); tumor=str(g.tumor.iloc[0])
  if g.tumor.nunique()!=1 or len(g)!=320: raise RuntimeError('ESS pool is not exactly one global 320-candidate tumor-mouse pool')
  qv=a20.candidate_rank_coordinate(g.Vmax_raw.to_numpy(float))
  ql=a20.candidate_rank_coordinate(g.Vlac_raw.to_numpy(float))
  tv=coords[(str(mouse_id),'Vmax')]; tl=coords[(str(mouse_id),'Vlac')]
  d1=(qv-tv)**2; d2=d1+(ql-tl)**2
  frozen1=g.weight.to_numpy(float)
  a2g=a2base.loc[a2base.mouse_id.astype(str).eq(str(mouse_id))].merge(g[idc],on=idc,how='right',validate='one_to_one',sort=False)
  if len(a2g)!=320 or a2g.weight.isna().any(): raise RuntimeError('A2 weights do not join to the exact A1 320-candidate order')
  frozen2=a2g.weight.to_numpy(float)
  for target in (10.,20.,40.):
   r1=a20.solve_temperature_for_ess(d1,target_ess=target)
   r2=a20.solve_temperature_for_ess(d2,target_ess=target)
   if target==20.:
    diff1=float(np.max(np.abs(r1.weights-frozen1))); diff2=float(np.max(np.abs(r2.weights-frozen2)))
    ess1=float(1/np.dot(frozen1,frozen1)); ess2=float(1/np.dot(frozen2,frozen2))
    parity.append({'mouse_id':str(mouse_id),'tumor':tumor,'a1_max_abs_weight_diff':diff1,'a2_max_abs_weight_diff':diff2,'a1_frozen_ess':ess1,'a1_reproduced_ess':r1.achieved_ess,'a2_frozen_ess':ess2,'a2_reproduced_ess':r2.achieved_ess})
    if diff1>1e-12 or diff2>1e-12 or abs(ess1-20)>1e-10 or abs(ess2-20)>1e-10: raise RuntimeError('ESS-20 original numerical parity failed; stopping Experiment 4')
   for arm,res in (('A1',r1),('A2',r2)):
    all_weights[(str(mouse_id),arm,target)]=(g[idc].copy(),res.weights.copy())
    for (alg,rna),idx in g.groupby(['algorithm','rna_context_key'],sort=True).groups.items():
     positions=g.index.get_indexer(idx); cw=res.weights[positions]
     # Conditional weights are normalized within each original 20-vector stratum.
     if len(cw)!=20: raise RuntimeError('conditional stratum does not contain 20 vectors')
     cw_sum=float(cw.sum()); has_support=bool(cw_sum>0); cw=(cw/cw_sum) if has_support else np.zeros_like(cw); cess=float(1/np.dot(cw,cw)) if has_support else float('nan')
     outputs.append({'mouse_id':str(mouse_id),'tumor':tumor,'arm':arm,'algorithm':str(alg),'rna_context_key':str(rna),'target_global_ess':target,'temperature':res.temperature,'global_ess':res.achieved_ess,'conditional_ess':cess,'conditional_positive_n':int(np.count_nonzero(cw>0)),'conditional_support':has_support})
 out.mkdir(parents=True,exist_ok=True)
 pd.DataFrame(parity).to_csv(out/'experiment4_ess20_parity.tsv',sep='\t',index=False)
 pd.DataFrame(outputs).to_csv(out/'experiment4_weight_support.tsv.gz',sep='\t',index=False,compression='gzip')
 # Product-weight ESS is the product of the two condition-specific conditional ESS values.
 ws=pd.DataFrame(outputs)
 prows=[]
 for target in (10.,20.,40.):
  sub=ws.loc[ws.target_global_ess.eq(target)]
  for r in data['evals'].itertuples():
   for arm in ('A1','A2'):
    left=sub.loc[sub.mouse_id.eq(str(r.ct2a_mouse)) & sub.arm.eq(arm) & sub.algorithm.eq(str(r.algorithm)) & sub.rna_context_key.eq(str(r.rna_context_key))]
    right=sub.loc[sub.mouse_id.eq(str(r.gl261_mouse)) & sub.arm.eq(arm) & sub.algorithm.eq(str(r.algorithm)) & sub.rna_context_key.eq(str(r.rna_context_key))]
    if len(left)!=1 or len(right)!=1: raise RuntimeError('evaluation conditional ESS join failed')
    prows.append({'evaluation_id':str(r.evaluation_id),'arm':arm,'target_global_ess':target,'ct2a_conditional_ess':float(left.conditional_ess.iloc[0]),'gl261_conditional_ess':float(right.conditional_ess.iloc[0]),'product_weight_ess':float(left.conditional_ess.iloc[0]*right.conditional_ess.iloc[0])})
 pd.DataFrame(prows).to_csv(out/'experiment4_product_ess.tsv.gz',sep='\t',index=False,compression='gzip')
 def weights_for(mouse_id,tumor,alg,rna,target,arm):
  ids,w=all_weights[(str(mouse_id),arm,float(target))]
  mask=ids.algorithm.astype(str).eq(str(alg)).to_numpy() & ids.rna_context_key.astype(str).eq(str(rna)).to_numpy()
  group=ids.loc[mask].sort_values('sample_index',kind='stable')
  vals=w[mask][np.argsort(group.sample_index.to_numpy(int),kind='stable')]
  return vals if len(vals)==20 else None
 truth=data['pairs'].loc[data['pairs'].pair_evaluable.eq('true')].drop_duplicates(['evaluation_id','ct2a_candidate_id','gl261_candidate_id'])
 common_by_target={}
 for target in (10.,20.,40.):
  keys=set()
  for pr in truth.itertuples(index=False):
   ev=data['evals'].set_index('evaluation_id').loc[str(pr.evaluation_id)]
   ok=True
   for arm in ('A1','A2'):
    for tumor,mouse,sample in [('CT2A',str(pr.ct2a_mouse),int(pr.ct2a_sample_index)),('GL261',str(pr.gl261_mouse),int(pr.gl261_sample_index))]:
     w=weights_for(mouse,tumor,str(ev.algorithm),str(ev.rna_context_key),target,arm)
     if w is None or w.sum()<=0 or (w.sum()-w[sample])<=0: ok=False
   if ok: keys.add((str(pr.evaluation_id),str(pr.ct2a_candidate_id),str(pr.gl261_candidate_id)))
  common_by_target[target]=keys
 all_case=[]; all_support=[]
 for target in (10.,20.,40.):
  for arm in ('A1','A2'):
   rows,supp=score_frozen_truth_pairs(data,weights_for,common_by_target[target],scenario='ESS_'+str(int(target)),arm=arm,target=target,outdir=out)
   all_case.extend(rows); all_support.extend(supp)
 pd.DataFrame(all_case).to_csv(out/'experiment4_gain_usefulness.tsv.gz',sep='\t',index=False,compression='gzip')
 pd.DataFrame(all_support).to_csv(out/'experiment4_pair_support.tsv.gz',sep='\t',index=False,compression='gzip')
 if all_case:
  cf=pd.DataFrame(all_case)
  paired=cf.pivot(index=['target_ess','reaction_id'],columns='arm',values=['usefulness','correct_gain_mean','wrong_gain_mean','correct_minus_wrong_mean'])
  paired.columns=['_'.join(col) for col in paired.columns]; paired=paired.reset_index()
  paired.to_csv(out/'experiment4_a1_a2_paired.tsv.gz',sep='\t',index=False,compression='gzip')
 # Recompute eta2/entropy from the same recalibrated 20x20 weights for each target.
 geoms=[]
 for target in (10.,20.,40.):
  for arm in ('A1','A2'):
   geoms.append(weighted_geometry_summary(data,weights_for,target,arm))
 geom=pd.concat(geoms,ignore_index=True)
 geom.to_csv(out/'experiment4_geometry_summary.tsv.gz',sep='\t',index=False,compression='gzip')
 rxreg=pd.read_csv(root/'outputs/dmi_bridge_pl2a_sign_only_prepare_v1/BRIDGEPL2A_REACTION_REGISTRY.tsv',sep='\t')
 split=pd.read_csv(root/'outputs/dmi_bridge_pl2c_geometry_utility_development_v1/BRIDGEPL2C_REACTION_SPLIT.tsv',sep='\t')
 assoc=[]
 casedf=pd.DataFrame(all_case)
 for target in (10.,20.,40.):
  for arm in ('A1','A2'):
   od=casedf.loc[casedf.target_ess.eq(target)&casedf.arm.eq(arm)]
   gm=geom.loc[geom.target_ess.eq(target)&geom.arm.eq(arm)]
   if len(od) and len(gm):
    temp=associations(od,gm.rename(columns={'eta2':'eta2','H_dir':'H_dir','D_dir':'D_dir'}),gm,rxreg,split) if False else None
    from scipy.stats import spearmanr
    joined=od.merge(gm,on='reaction_id',validate='one_to_one').merge(rxreg[['reaction_id','subsystem']],on='reaction_id',validate='one_to_one').merge(split[['reaction_id','analysis_split']],on='reaction_id',validate='one_to_one')
    for scope,sub in [('ALL',joined),*[(x,joined.loc[joined.analysis_split.eq(x)]) for x in sorted(joined.analysis_split.unique())]]:
     for desc in ('eta2','H_dir'):
      z=sub[[desc,'usefulness']].dropna(); assoc.append({'target_ess':target,'arm':arm,'reaction_split':scope,'descriptor':desc,'n_reactions':len(z),'spearman_rho':float(spearmanr(z[desc],z.usefulness).statistic) if len(z)>2 else float('nan')})
     z=sub[['eta2','H_dir','usefulness']].dropna()
     if len(z)>2:
      from scipy.stats import rankdata
      xr=rankdata(z.eta2); yr=rankdata(z.usefulness); hr=rankdata(z.H_dir); X0=np.c_[np.ones(len(z)),hr]; rx=xr-X0@np.linalg.lstsq(X0,xr,rcond=None)[0]; ry=yr-X0@np.linalg.lstsq(X0,yr,rcond=None)[0]
      assoc.append({'target_ess':target,'arm':arm,'reaction_split':scope,'descriptor':'eta2_increment','n_reactions':len(z),'partial_rank_eta2_given_entropy':float(np.corrcoef(rx,ry)[0,1])})
 pd.DataFrame(assoc).to_csv(out/'experiment4_eta2_usefulness_associations.tsv',sep='\t',index=False)
 manifest={'status':'ESS20_PARITY_PASS_OUTCOME_RECOMPUTED','targets':[10,20,40],'lambda':LAM,'truth_registry_sha256':EXPECTED['outputs/dmi_bridge_pl2a_sign_only_prepare_v1/BRIDGEPL2A_TRUTH_PAIR_REGISTRY.tsv'],'distinct_original_truth_pairs':int(len(truth)),'paired_eligible_truth_pairs_by_target':{str(int(t)):len(k) for t,k in common_by_target.items()},'source_hashes':EXPECTED|{A1_CURRENT:A1_CURRENT_SHA,A1_MANIFEST:A1_MANIFEST_SHA,ARCHIVED_KERNEL:ARCHIVED_KERNEL_SHA},'code_sha256':{p:sha256(root/p) for p in CODE_PATHS},'calibration':{'operator':'dmi_bridge_a20_dual_anchor_core_v1.solve_temperature_for_ess','global_pool_n':320,'targets':[10,20,40],'low':1e-12,'initial_high':1.0,'max_high':1e6,'geometric_bisection_iterations':80},'source_root':str(root),'note':'A1-A2 contrasts are calibration sensitivities, not isolated effects of adding a measurement.'}
 (out/'experiment4_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')

def main():
 ap=argparse.ArgumentParser(); ap.add_argument('experiment',choices=['3a','3b','4']); ap.add_argument('--upstream-root',type=Path,default=Path('/home/pty/work/project_MRI_brain_tumor')); ap.add_argument('--output-root',type=Path,required=True); ap.add_argument('--pilot',action='store_true'); a=ap.parse_args()
 sys.path.insert(0,str(a.upstream_root/'scripts'))
 data=load(a.upstream_root); data['root']=a.upstream_root
 # Source core modules imported by load; bind explicitly for numerical calls.
 import dmi_bridge_pl1_predictability_core_v1 as pl1
 import dmi_bridge_pl2_sign_only_core_v1 as pl2
 import dmi_bridge_a20_dual_anchor_core_v1 as a20
 data.update(pl1=pl1,pl2=pl2,a20mod=a20)
 if a.experiment=='4': run4(data,a.output_root)
 elif a.experiment=='3a': run3a(data,a.output_root,pilot=a.pilot)
 elif a.experiment=='3b': run3b(data,a.output_root)
 else: raise RuntimeError(f'Experiment {a.experiment} runner is not yet staged')

if __name__=='__main__': main()
