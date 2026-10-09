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

def vector_outcomes(left:np.ndarray,right:np.ndarray,wl:np.ndarray,wr:np.ndarray,truth_delta:np.ndarray,block:int=64,core=None)->dict:
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
 # The frozen scalar implementation is authoritative at the 1e-12 decision
 # boundary; vectorized reductions can differ by a few ulps after cancellation.
 if core is not None:
  borderline=non & ((np.abs(correct-ATOL)<1e-10)|(np.abs(wrong-ATOL)<1e-10)|(np.abs(info-ATOL)<1e-10))
  for j in np.flatnonzero(borderline):
   delta=(left[:,j,None]-right[None,:,j]).reshape(-1)
   ref=scalar_outcome(delta,joint,float(truth_delta[j]),core)['lambda_results'][0]
   correct[j]=float(ref['correct_abs_error_gain']);wrong[j]=float(ref['wrong_abs_error_gain'])
   info[j]=float(correct[j]-ref['random_sign_expected_abs_error_gain'])
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
  rows,supp,_=score_frozen_truth_pairs(data,wf,common,scenario='truth_excluded_q10_q50_q90',arm=arm,target=20.,outdir=out)
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
    a=geo[(arm,r,str(p.evaluation_id))]
    for k in ('eta2','H_dir','D_dir','K_dir','non_tie_coverage'):
     if math.isfinite(float(z[k][j])):
      a[k]+=float(z[k][j])
      if k in ('eta2','H_dir','D_dir'): a[f'n_{k}']+=1
    for k in ('eta2','H_dir','D_dir'):
     if math.isfinite(float(z20[k][j])):
      a[f'{k}_20x20']+=float(z20[k][j]); a[f'n_{k}_20x20']+=1
    a['n']+=1
 # First average descriptors within each evaluation over distinct non-tie pairs;
 # then average evaluations equally at the reaction level.
 eval_geo=[]
 for (arm,r,evaluation_id),z in geo.items():
  n=int(z['n']); eval_geo.append({'arm':arm,'reaction_id':r,'evaluation_id':evaluation_id,'n_non_tie_truth_pairs':n,**{k:(z[k]/z[f'n_{k}'] if z[f'n_{k}'] else np.nan) for k in ('eta2','H_dir','D_dir')},**{f'{k}_20x20':(z[f'{k}_20x20']/z[f'n_{k}_20x20'] if z[f'n_{k}_20x20'] else np.nan) for k in ('eta2','H_dir','D_dir')}})
 eg=pd.DataFrame(eval_geo); georows=[]
 for (arm,r),g in eg.groupby(['arm','reaction_id'],sort=True):
  row={'arm':arm,'reaction_id':r,'n_geometry_evaluations':int(g.evaluation_id.nunique()),'n_distinct_non_tie_truth_pairs':int(g.n_non_tie_truth_pairs.sum())}
  for k in ('eta2','H_dir','D_dir','eta2_20x20','H_dir_20x20','D_dir_20x20'):
   row[k]=float(g[k].mean()) if g[k].notna().any() else np.nan;row[f'n_finite_{k}_evaluations']=int(g[k].notna().sum())
  georows.append(row)
 eg.to_csv(out/'experiment3a_evaluation_geometry.tsv.gz',sep='\t',index=False,compression='gzip')
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
 out_n=aligned.n_non_tie_truth_pairs.fillna(0).to_numpy(int); geo_n=aligned.n_distinct_non_tie_truth_pairs.fillna(0).to_numpy(int)
 if not np.array_equal(out_n,geo_n): raise RuntimeError('truth-pair geometry/outcome aggregation denominator mismatch')
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

def aggregate_evaluation_reaction_records(edf:pd.DataFrame):
 """Frozen PL2-style distinct-pair → evaluation → equal-evaluation reaction means."""
 if not len(edf): return pd.DataFrame(),[]
 edf=edf.drop_duplicates(['evaluation_id','reaction_id','truth_pair_id'],keep='first')
 ev=edf.groupby(['scenario','arm','target_ess','evaluation_id','reaction_id'],sort=True,as_index=False).agg(n_distinct_truth_pairs=('truth_pair_id','nunique'),n_non_tie_truth_pairs=('truth_tie',lambda x:int((~x).sum())),n_truth_ties=('truth_tie','sum'),usefulness=('usefulness','mean'),correct_gain_mean=('correct_gain','mean'),wrong_gain_mean=('wrong_gain','mean'),information_advantage_mean=('information_advantage','mean'),correct_minus_wrong_mean=('correct_minus_wrong','mean'))
 ev['response_defined']=ev.n_non_tie_truth_pairs.gt(0); rows=[]
 endpoints=['usefulness','correct_gain_mean','wrong_gain_mean','information_advantage_mean','correct_minus_wrong_mean']
 for (sc,a,t,r),g in ev.groupby(['scenario','arm','target_ess','reaction_id'],sort=True):
  row={'scenario':sc,'arm':a,'target_ess':t,'reaction_id':r,'n_eligible_evaluations':int(g.evaluation_id.nunique()),'n_defined_response_evaluations':int(g.response_defined.sum()),'n_distinct_truth_pairs':int(g.n_distinct_truth_pairs.sum()),'n_non_tie_truth_pairs':int(g.n_non_tie_truth_pairs.sum()),'n_truth_ties':int(g.n_truth_ties.sum())}
  for c in endpoints: row[c]=float(g[c].mean()) if g[c].notna().any() else np.nan; row[f'n_finite_{c}_evaluations']=int(g[c].notna().sum())
  rows.append(row)
 return ev,rows

def reaction_summary_from_evaluations(ev:pd.DataFrame):
 """Pandas equal-evaluation PL2D-style reaction means and finite denominators."""
 if not len(ev):return pd.DataFrame()
 grouped=ev.groupby(['scenario','arm','target_ess','reaction_id'],sort=True)
 out=grouped.agg(n_eligible_evaluations=('evaluation_id','nunique'),n_defined_response_evaluations=('response_defined','sum'),n_distinct_truth_pairs=('n_distinct_truth_pairs','sum'),n_non_tie_truth_pairs=('n_non_tie_truth_pairs','sum'),n_truth_ties=('n_truth_ties','sum'),usefulness=('usefulness','mean'),correct_gain_mean=('correct_gain_mean','mean'),wrong_gain_mean=('wrong_gain_mean','mean'),information_advantage_mean=('information_advantage_mean','mean'),correct_minus_wrong_mean=('correct_minus_wrong_mean','mean'),n_finite_usefulness_evaluations=('usefulness','count'),n_finite_correct_gain_mean_evaluations=('correct_gain_mean','count'),n_finite_wrong_gain_mean_evaluations=('wrong_gain_mean','count'),n_finite_information_advantage_mean_evaluations=('information_advantage_mean','count'),n_finite_correct_minus_wrong_mean_evaluations=('correct_minus_wrong_mean','count')).reset_index()
 return out

def frozen_case_evaluations(root:Path,pairs:pd.DataFrame,arm:str,scenario:str,target:float,keys=None):
 """Reaggregate frozen PL2B/A22 pair cases, preserving their exact tol decisions."""
 sys.path.insert(0,str(root/'scripts'))
 import dmi_bridge_pl2b_production_v1 as pl2b
 import json as _json
 rxns=np.load(root/'outputs/secondary_anchor_exploratory/dmi_sa_six_arm_current_20260926/reconstructed/flux_cache.npz',allow_pickle=False)['rxns'].astype(str).tolist()
 targets=[r for r in rxns if r not in {'HEX1','LDH_L'}];nr=len(targets)
 if nr!=4179:raise RuntimeError('frozen case aggregation requires the 4,179-reaction weak-cue universe')
 p=pairs.loc[pairs.pair_evaluable.astype(str).eq('true')].drop_duplicates(['evaluation_id','ct2a_candidate_id','gl261_candidate_id']).copy()
 tuples=list(zip(p.evaluation_id.astype(str),p.ct2a_candidate_id.astype(str),p.gl261_candidate_id.astype(str)))
 if keys is not None:
  allowed=set(keys);mask=[k in allowed for k in tuples];p=p.loc[mask].copy();tuples=[k for k,m in zip(tuples,mask) if m]
 pair_ids=[pl2b.truth_pair_id(*k) for k in tuples]
 if len(set(pair_ids))!=len(pair_ids):raise RuntimeError('frozen truth-pair IDs are not unique')
 pairset=set(pair_ids);eids=sorted(p.evaluation_id.astype(str).unique());emap={x:i for i,x in enumerate(eids)};rmap={x:i for i,x in enumerate(targets)};size=len(eids)*nr
 sums={k:np.zeros(size,dtype=np.float64) for k in ('n_distinct','n_non','n_tie','correct','wrong','info','use')}
 folder='outputs/dmi_bridge_pl2b_sign_only_utility_v1' if arm=='A1' else 'outputs/dmi_bridge_a22_dual_anchor_sign_only_utility_v1'
 files=sorted((root/folder).glob('BRIDGEPL2B_CASE_OUTCOMES.part-*.tsv.xz' if arm=='A1' else 'BRIDGEA22_CASE_OUTCOMES.part-*.tsv.xz'))
 if not files:raise RuntimeError(f'frozen PL2/A22 case parts missing for {arm}')
 cols=['truth_pair_id','evaluation_id','reaction_id','truth_status','lambda_0_25_correct_absolute_error_gain','lambda_0_25_wrong_absolute_error_gain','lambda_0_25_random_sign_expected_gain']
 for part in files:
  for d in pd.read_csv(part,sep='\t',compression='xz',usecols=cols,chunksize=120000,float_precision='round_trip'):
   d=d.loc[d.truth_pair_id.isin(pairset)&d.reaction_id.isin(rmap)]
   if d.empty:continue
   tie=d.truth_status.eq('TRUTH_TIE').to_numpy();c=d.lambda_0_25_correct_absolute_error_gain.to_numpy(float);w=d.lambda_0_25_wrong_absolute_error_gain.to_numpy(float);random=d.lambda_0_25_random_sign_expected_gain.to_numpy(float);info=c-random;use=(~tie)&(c>ATOL)&(info>ATOL)
   idx=d.evaluation_id.map(emap).to_numpy(int)*nr+d.reaction_id.map(rmap).to_numpy(int)
   for name,val in [('n_distinct',np.ones(len(d))),('n_non',(~tie).astype(float)),('n_tie',tie.astype(float)),('correct',np.where(tie,0.,c)),('wrong',np.where(tie,0.,w)),('info',np.where(tie,0.,info)),('use',use.astype(float))]:sums[name]+=np.bincount(idx,weights=val,minlength=size)
 total=float(sums['n_distinct'].sum());expected=float(len(pairset)*nr)
 if total!=expected:raise RuntimeError(f'frozen {arm} case support mismatch: {total} != {expected}')
 ix=np.flatnonzero(sums['n_distinct']>0);ei=ix//nr;ri=ix%nr;non=sums['n_non'][ix]
 def mean(name):return np.divide(sums[name][ix],non,out=np.full(len(ix),np.nan),where=non>0)
 return pd.DataFrame({'scenario':scenario,'arm':arm,'target_ess':target,'evaluation_id':[eids[i] for i in ei],'reaction_id':[targets[i] for i in ri],'n_distinct_truth_pairs':sums['n_distinct'][ix].astype(int),'n_non_tie_truth_pairs':non.astype(int),'n_truth_ties':sums['n_tie'][ix].astype(int),'response_defined':non>0,'usefulness':mean('use'),'correct_gain_mean':mean('correct'),'wrong_gain_mean':mean('wrong'),'information_advantage_mean':mean('info'),'correct_minus_wrong_mean':mean('correct')-mean('wrong')})

def frozen_case_evaluations_groups(root:Path,pairs:pd.DataFrame,arm:str,target:float,key_groups:dict):
 """Aggregate multiple overlapping frozen truth-pair populations in one source scan."""
 import dmi_bridge_pl2b_production_v1 as pl2b
 valid=pairs.loc[pairs.pair_evaluable.astype(str).eq('true')].drop_duplicates(['evaluation_id','ct2a_candidate_id','gl261_candidate_id']).copy()
 tuples=list(zip(valid.evaluation_id.astype(str),valid.ct2a_candidate_id.astype(str),valid.gl261_candidate_id.astype(str)))
 pid_groups={}
 for group,keys in key_groups.items():
  for key in keys:
   pid=pl2b.truth_pair_id(*key);pid_groups.setdefault(pid,[]).append(group)
 pid_set=set(pid_groups)
 rxns=np.load(root/'outputs/secondary_anchor_exploratory/dmi_sa_six_arm_current_20260926/reconstructed/flux_cache.npz',allow_pickle=False)['rxns'].astype(str).tolist();targets=[r for r in rxns if r not in {'HEX1','LDH_L'}];nr=len(targets)
 if nr!=4179:raise RuntimeError('frozen grouped aggregation requires 4,179 weak-cue reactions')
 eids=sorted(valid.evaluation_id.astype(str).unique());emap={x:i for i,x in enumerate(eids)};rmap={x:i for i,x in enumerate(targets)};size=len(eids)*nr
 names=('n_distinct','n_non','n_tie','correct','wrong','info','use');accs={g:{n:np.zeros(size,dtype=np.float64) for n in names} for g in key_groups}
 folder='outputs/dmi_bridge_pl2b_sign_only_utility_v1' if arm=='A1' else 'outputs/dmi_bridge_a22_dual_anchor_sign_only_utility_v1';pattern='BRIDGEPL2B_CASE_OUTCOMES.part-*.tsv.xz' if arm=='A1' else 'BRIDGEA22_CASE_OUTCOMES.part-*.tsv.xz'
 cols=['truth_pair_id','evaluation_id','reaction_id','truth_status','lambda_0_25_correct_absolute_error_gain','lambda_0_25_wrong_absolute_error_gain','lambda_0_25_random_sign_expected_gain']
 for part in sorted((root/folder).glob(pattern)):
  for d in pd.read_csv(part,sep='\t',compression='xz',usecols=cols,chunksize=120000,float_precision='round_trip'):
   d=d.loc[d.truth_pair_id.isin(pid_set)&d.reaction_id.isin(rmap)]
   if d.empty:continue
   tie=d.truth_status.eq('TRUTH_TIE').to_numpy();c=d.lambda_0_25_correct_absolute_error_gain.to_numpy(float);w=d.lambda_0_25_wrong_absolute_error_gain.to_numpy(float);info=c-d.lambda_0_25_random_sign_expected_gain.to_numpy(float);use=(~tie)&(c>ATOL)&(info>ATOL);idx=d.evaluation_id.map(emap).to_numpy(int)*nr+d.reaction_id.map(rmap).to_numpy(int)
   for pid,positions in d.groupby('truth_pair_id',sort=False).indices.items():
    ix=idx[np.asarray(positions,dtype=int)];ti=tie[np.asarray(positions,dtype=int)];cc=c[np.asarray(positions,dtype=int)];ww=w[np.asarray(positions,dtype=int)];ii=info[np.asarray(positions,dtype=int)];uu=use[np.asarray(positions,dtype=int)]
    for group in pid_groups[pid]:
     a=accs[group]
     for n,v in [('n_distinct',np.ones(len(ix))),('n_non',(~ti).astype(float)),('n_tie',ti.astype(float)),('correct',np.where(ti,0.,cc)),('wrong',np.where(ti,0.,ww)),('info',np.where(ti,0.,ii)),('use',uu.astype(float))]:a[n]+=np.bincount(ix,weights=v,minlength=size)
 result={}
 for group,a in accs.items():
  expected=len(key_groups[group])*nr
  if a['n_distinct'].sum()!=expected:raise RuntimeError(f'frozen support mismatch for {group}')
  ix=np.flatnonzero(a['n_distinct']>0);ei=ix//nr;ri=ix%nr;nn=a['n_non'][ix]
  def mean(n):return np.divide(a[n][ix],nn,out=np.full(len(ix),np.nan),where=nn>0)
  result[group]=pd.DataFrame({'scenario':group,'arm':arm,'target_ess':target,'evaluation_id':[eids[i] for i in ei],'reaction_id':[targets[i] for i in ri],'n_distinct_truth_pairs':a['n_distinct'][ix].astype(int),'n_non_tie_truth_pairs':nn.astype(int),'n_truth_ties':a['n_tie'][ix].astype(int),'response_defined':nn>0,'usefulness':mean('use'),'correct_gain_mean':mean('correct'),'wrong_gain_mean':mean('wrong'),'information_advantage_mean':mean('info'),'correct_minus_wrong_mean':mean('correct')-mean('wrong')})
 return result

def restore3a_frozen_outcomes(data,out):
 out=Path(out);pairs=data['pairs'].loc[data['pairs'].pair_evaluable.astype(str).eq('true')].drop_duplicates(['evaluation_id','ct2a_candidate_id','gl261_candidate_id'])
 if len(pairs)!=836:raise RuntimeError('3A frozen matched population is not 836 distinct pairs')
 computed=[];restored=[]
 for arm in ('A1','A2'):
  current=pd.read_csv(out/f'truth_excluded_q10_q50_q90_{arm}_evaluation_reaction.tsv.gz',sep='\t')
  frozen=frozen_case_evaluations(data['root'],pairs,arm,'truth_excluded_q10_q50_q90',20.)
  joined=current.merge(frozen,on=['evaluation_id','reaction_id'],suffixes=('_computed','_frozen'),validate='one_to_one')
  for a,b in [('n_distinct_truth_pairs','n_distinct_truth_pairs'),('n_non_tie_truth_pairs','n_non_tie_truth_pairs'),('n_truth_ties','n_truth_ties'),('correct_gain_mean','correct_gain_mean'),('wrong_gain_mean','wrong_gain_mean'),('information_advantage_mean','information_advantage_mean'),('usefulness','usefulness')]:
   x=joined[f'{a}_computed'].to_numpy(float);y=joined[f'{b}_frozen'].to_numpy(float);d=np.abs(x-y);computed.append({'arm':arm,'field':a,'n_rows':len(joined),'max_abs_difference':float(np.nanmax(d)),'n_different':int(np.sum(d>ATOL)),'tolerance':ATOL,'reference':'frozen PL2B' if arm=='A1' else 'frozen A22'})
  restored.append(frozen)
  frozen.to_csv(out/f'truth_excluded_q10_q50_q90_{arm}_evaluation_reaction.tsv.gz',sep='\t',index=False,compression='gzip')
 ev=pd.concat(restored,ignore_index=True);rx=reaction_summary_from_evaluations(ev)
 rx.to_csv(out/'experiment3a_usefulness.tsv.gz',sep='\t',index=False,compression='gzip')
 pd.DataFrame(computed).to_csv(out/'experiment3a_recomputed_vs_frozen_case_audit.tsv',sep='\t',index=False)
 if not np.array_equal(ev.n_distinct_truth_pairs.to_numpy(int),ev.n_non_tie_truth_pairs.to_numpy(int)+ev.n_truth_ties.to_numpy(int)):raise RuntimeError('frozen case pair/tie accounting identity failed')

def restore3b_frozen_baseline(data,out):
 """Restore the unchanged q10/q50/q90 scenario from frozen PL2B/A22 cases."""
 out=Path(out);pairs=data['pairs'].loc[data['pairs'].pair_evaluable.astype(str).eq('true')].drop_duplicates(['evaluation_id','ct2a_candidate_id','gl261_candidate_id'])
 if len(pairs)!=836:raise RuntimeError('3B original frozen baseline must contain 836 pairs')
 restored=[];audit=[]
 for arm in ('A1','A2'):
  path=out/f'original_q10_q50_q90_{arm}_evaluation_reaction.tsv.gz';current=pd.read_csv(path,sep='\t')
  frozen=frozen_case_evaluations(data['root'],pairs,arm,'original_q10_q50_q90',20.)
  j=current.merge(frozen,on=['evaluation_id','reaction_id'],suffixes=('_prior','_frozen'),validate='one_to_one')
  for col in ('n_distinct_truth_pairs','n_non_tie_truth_pairs','n_truth_ties','usefulness','correct_gain_mean','wrong_gain_mean','information_advantage_mean'):
   delta=np.abs(j[f'{col}_prior'].to_numpy(float)-j[f'{col}_frozen'].to_numpy(float));audit.append({'arm':arm,'endpoint':col,'n_evaluation_reaction_rows':len(j),'n_changed_over_1e-12':int(np.sum(delta>ATOL)),'max_abs_difference':float(np.nanmax(delta)),'reference':'frozen PL2B' if arm=='A1' else 'frozen A22'})
  frozen.to_csv(path,sep='\t',index=False,compression='gzip');restored.append(frozen)
 ev=pd.concat(restored,ignore_index=True);rx=reaction_summary_from_evaluations(ev)
 full=pd.read_csv(out/'experiment3b_usefulness.tsv.gz',sep='\t')
 full=full.loc[full.scenario.ne('original_q10_q50_q90')]
 full=pd.concat([full,rx],ignore_index=True)
 full.to_csv(out/'experiment3b_usefulness.tsv.gz',sep='\t',index=False,compression='gzip')
 pd.DataFrame(audit).to_csv(out/'experiment3b_recomputed_vs_frozen_case_audit.tsv',sep='\t',index=False)
 # Refit the paired-arm association summaries on identical reaction IDs.
 rxreg=pd.read_csv(data['root']/'outputs/dmi_bridge_pl2a_sign_only_prepare_v1/BRIDGEPL2A_REACTION_REGISTRY.tsv',sep='\t')
 split=pd.read_csv(data['root']/'outputs/dmi_bridge_pl2c_geometry_utility_development_v1/BRIDGEPL2C_REACTION_SPLIT.tsv',sep='\t')
 geom=[]
 for arm,path in [('A1',data['root']/'outputs/dmi_bridge_pl1_predictability_landscape_v1_resumable_v1/BRIDGEPL1_REACTION_SUMMARY.tsv'),('A2',data['root']/'outputs/dmi_bridge_a21_dual_anchor_geometry_v1/BRIDGEA21_REACTION_SUMMARY.tsv.xz')]:
  f=pd.read_csv(path,sep='\t',compression='infer');geom.append(f[['reaction_id','sign_magnitude_eta2_median','directional_entropy3_median','dominant_sign_mass_median']].assign(arm=arm).rename(columns={'sign_magnitude_eta2_median':'eta2','directional_entropy3_median':'H_dir','dominant_sign_mass_median':'D_dir'}))
 geo=pd.concat(geom,ignore_index=True);metrics=[]
 for scenario,g in full.groupby('scenario',sort=True):
  m=paired_arm_associations(g,geo,rxreg,split)
  if len(m):
   if 'scenario' not in m.columns:m.insert(0,'scenario',scenario)
   metrics.extend(m.to_dict('records'))
 pd.DataFrame(metrics).to_csv(out/'experiment3b_associations.tsv',sep='\t',index=False)

def finalize3b_associations(data,out):
 out=Path(out);full=pd.read_csv(out/'experiment3b_usefulness.tsv.gz',sep='\t');root=data['root'];rxreg=pd.read_csv(root/'outputs/dmi_bridge_pl2a_sign_only_prepare_v1/BRIDGEPL2A_REACTION_REGISTRY.tsv',sep='\t');split=pd.read_csv(root/'outputs/dmi_bridge_pl2c_geometry_utility_development_v1/BRIDGEPL2C_REACTION_SPLIT.tsv',sep='\t');geom=[]
 for arm,path in [('A1',root/'outputs/dmi_bridge_pl1_predictability_landscape_v1_resumable_v1/BRIDGEPL1_REACTION_SUMMARY.tsv'),('A2',root/'outputs/dmi_bridge_a21_dual_anchor_geometry_v1/BRIDGEA21_REACTION_SUMMARY.tsv.xz')]:
  f=pd.read_csv(path,sep='\t',compression='infer');geom.append(f[['reaction_id','sign_magnitude_eta2_median','directional_entropy3_median','dominant_sign_mass_median']].assign(arm=arm).rename(columns={'sign_magnitude_eta2_median':'eta2','directional_entropy3_median':'H_dir','dominant_sign_mass_median':'D_dir'}))
 geo=pd.concat(geom,ignore_index=True);metrics=[]
 for scenario,g in full.groupby('scenario',sort=True):
  m=paired_arm_associations(g,geo,rxreg,split)
  if len(m):
   if 'scenario' not in m.columns:m.insert(0,'scenario',scenario)
   metrics.extend(m.to_dict('records'))
 pd.DataFrame(metrics).to_csv(out/'experiment3b_associations.tsv',sep='\t',index=False)

def restore4_ess20_frozen(data,out):
 """Use official frozen case outcomes at ESS20 for both populations and arms."""
 out=Path(out);support=pd.read_csv(out/'experiment4_pair_support.tsv.gz',sep='\t');all_eval=[];audit=[]
 current=pd.read_csv(out/'experiment4_gain_usefulness.tsv.gz',sep='\t')
 for arm in ('A1','A2'):
  groups={}
  for scenario in ('ESS_20_eligible','ESS_20_fixed_common'):
   s=support.loc[support.scenario.eq(scenario)&support.arm.eq(arm)&support.eligible]
   groups[scenario]=set(zip(s.evaluation_id.astype(str),s.ct2a_candidate_id.astype(str),s.gl261_candidate_id.astype(str)))
  frozen_groups=frozen_case_evaluations_groups(data['root'],data['pairs'],arm,20.,groups)
  for scenario,frozen in frozen_groups.items():
   s=support.loc[support.scenario.eq(scenario)&support.arm.eq(arm)&support.eligible]
   path=out/f'{scenario}_{arm}_evaluation_reaction.tsv.gz';prior=pd.read_csv(path,sep='\t')
   j=prior.merge(frozen,on=['evaluation_id','reaction_id'],suffixes=('_prior','_frozen'),validate='one_to_one')
   for col in ('n_distinct_truth_pairs','n_non_tie_truth_pairs','n_truth_ties','usefulness','correct_gain_mean','wrong_gain_mean','information_advantage_mean'):
    d=np.abs(j[f'{col}_prior'].to_numpy(float)-j[f'{col}_frozen'].to_numpy(float));audit.append({'scenario':scenario,'arm':arm,'endpoint':col,'n_rows':len(j),'n_changed_over_1e-12':int(np.sum(d>ATOL)),'max_abs_difference':float(np.nanmax(d))})
   frozen.to_csv(path,sep='\t',index=False,compression='gzip');all_eval.append(frozen)
   new=reaction_summary_from_evaluations(frozen)
   current=current.loc[~(current.scenario.eq(scenario)&current.arm.eq(arm))]
   current=pd.concat([current,new],ignore_index=True)
 current.to_csv(out/'experiment4_gain_usefulness.tsv.gz',sep='\t',index=False,compression='gzip')
 pd.DataFrame(audit).to_csv(out/'experiment4_recomputed_vs_frozen_case_audit.tsv',sep='\t',index=False)

def fixed_common_population(common_by_target:dict):
 """Exact identity intersection over all ESS targets and anchor arms."""
 if not common_by_target: return set()
 return set.intersection(*(set(v) for v in common_by_target.values()))

def score_frozen_truth_pairs(data, weights_for, common_keys, *, scenario:str, arm:str, target:float, outdir:Path, truth_pairs=None):
 """Score distinct pairs per evaluation, then weight defined evaluations equally."""
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
 accum={}; nr=len(targets)
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
  vec=vector_outcomes(left[:,[0]],right[:,[0]],lw,rw,truth_delta[[0]],core=data['pl2'])
  delta=(left[:,0,None]-right[None,:,0]).reshape(-1); joint=(lw[:,None]*rw[None,:]).reshape(-1)
  scalar=scalar_outcome(delta,joint,float(truth_delta[0]),data['pl2'])
  if int(vec['truth_direction'][0])!=int(scalar['truth_direction']): raise RuntimeError('vector/scalar truth-direction parity failed')
  lam=scalar.get('lambda_results',[])
  if lam and abs(float(vec['correct_gain'][0])-float(lam[0]['correct_abs_error_gain']))>1e-10: raise RuntimeError('vector/scalar gain parity failed')
  values=vector_outcomes(left,right,lw,rw,truth_delta,core=data['pl2'])
  non=np.abs(truth_delta)>ATOL
  z=accum.setdefault(casekey[0],{'n_distinct':np.zeros(nr,dtype=np.int32),'n_non_tie':np.zeros(nr,dtype=np.int32),'n_ties':np.zeros(nr,dtype=np.int32),**{k:np.zeros(nr,dtype=np.float64) for k in ('usefulness','correct_gain_mean','wrong_gain_mean','information_advantage_mean','correct_minus_wrong_mean')}})
  z['n_distinct']+=1;z['n_non_tie'][non]+=1;z['n_ties'][~non]+=1
  z['usefulness'][non]+=values['usefulness'][non];z['correct_gain_mean'][non]+=values['correct_gain'][non];z['wrong_gain_mean'][non]+=values['wrong_gain'][non];z['information_advantage_mean'][non]+=values['information_gain'][non];z['correct_minus_wrong_mean'][non]+=(values['correct_gain'][non]-values['wrong_gain'][non])
  support.append({'scenario':scenario,'arm':arm,'target_ess':target,'evaluation_id':str(pr.evaluation_id),'ct2a_candidate_id':casekey[1],'gl261_candidate_id':casekey[2],'eligible':True,'reason':''})
 frames=[]
 for eid,z in accum.items():
  frame={'scenario':scenario,'arm':arm,'target_ess':target,'evaluation_id':eid,'reaction_id':targets,'n_distinct_truth_pairs':z['n_distinct'],'n_non_tie_truth_pairs':z['n_non_tie'],'n_truth_ties':z['n_ties']}
  for k in ('usefulness','correct_gain_mean','wrong_gain_mean','information_advantage_mean','correct_minus_wrong_mean'):
   frame[k]=np.divide(z[k],z['n_non_tie'],out=np.full(nr,np.nan),where=z['n_non_tie']>0)
  frames.append(pd.DataFrame(frame))
 ev=pd.concat(frames,ignore_index=True) if frames else pd.DataFrame(); ev['response_defined']=ev.n_non_tie_truth_pairs.gt(0) if len(ev) else pd.Series(dtype=bool)
 if len(ev):
  grouped=ev.groupby(['scenario','arm','target_ess','reaction_id'],sort=True)
  summary=grouped.agg(n_eligible_evaluations=('evaluation_id','nunique'),n_defined_response_evaluations=('response_defined','sum'),n_distinct_truth_pairs=('n_distinct_truth_pairs','sum'),n_non_tie_truth_pairs=('n_non_tie_truth_pairs','sum'),n_truth_ties=('n_truth_ties','sum'),usefulness=('usefulness','mean'),correct_gain_mean=('correct_gain_mean','mean'),wrong_gain_mean=('wrong_gain_mean','mean'),information_advantage_mean=('information_advantage_mean','mean'),correct_minus_wrong_mean=('correct_minus_wrong_mean','mean'),n_finite_usefulness_evaluations=('usefulness','count'),n_finite_correct_gain_mean_evaluations=('correct_gain_mean','count'),n_finite_wrong_gain_mean_evaluations=('wrong_gain_mean','count'),n_finite_information_advantage_mean_evaluations=('information_advantage_mean','count'),n_finite_correct_minus_wrong_mean_evaluations=('correct_minus_wrong_mean','count')).reset_index();rows=summary.to_dict('records')
 else: rows=[]
 outdir.mkdir(parents=True,exist_ok=True)
 ev.to_csv(outdir/f'{scenario}_{arm}_evaluation_reaction.tsv.gz',sep='\t',index=False,compression='gzip')
 return rows,support,ev

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
 split=split.copy()
 split['analysis_split']=split.analysis_split.astype(str).str.upper().map({'DEVELOPMENT':'development','CONFIRMATION_HOLDOUT':'confirmation'})
 if split.analysis_split.isna().any(): raise RuntimeError('unrecognized frozen reaction split label')
 results=[]
 for arm,geom in [('A1',geom_a1),('A2',geom_a2)]:
  gm=geom.loc[geom.arm.eq(arm)].drop(columns=['arm'],errors='ignore') if 'arm' in geom.columns else geom
  dat=outcomes.loc[outcomes.arm.eq(arm)].merge(gm,on='reaction_id',validate='one_to_one').merge(reaction_registry[['reaction_id','subsystem']],on='reaction_id',validate='one_to_one').merge(split[['reaction_id','analysis_split']],on='reaction_id',validate='one_to_one')
  dat=dat[['reaction_id','eta2','H_dir','usefulness','subsystem','analysis_split']].dropna().copy()
  dev=dat.loc[dat.analysis_split.eq('development')].reset_index(drop=True); conf=dat.loc[dat.analysis_split.eq('confirmation')].reset_index(drop=True)
  def partial(z):
   if len(z)<4:return float('nan')
   xr=rankdata(z.eta2);yr=rankdata(z.usefulness);hr=rankdata(z.H_dir);X=np.c_[np.ones(len(z)),hr]
   return float(np.corrcoef(xr-X@np.linalg.lstsq(X,xr,rcond=None)[0],yr-X@np.linalg.lstsq(X,yr,rcond=None)[0])[0,1])
  if len(dev)>4 and dev.subsystem.nunique()>=5:
   folds=list(GroupKFold(n_splits=5).split(dev[['H_dir','eta2']],dev.usefulness,dev.subsystem.astype(str)))
   cvpred={name:np.full(len(dev),np.nan) for name in ('entropy','entropy_plus_eta2')}
   for fold,(tr,te) in enumerate(folds):
    if set(dev.subsystem.iloc[tr]) & set(dev.subsystem.iloc[te]): raise RuntimeError('pathway GroupKFold leakage')
    for name,cols in [('entropy',['H_dir']),('entropy_plus_eta2',['H_dir','eta2'])]:
     model=make_pipeline(StandardScaler(),Ridge(alpha=1.0));model.fit(dev.loc[tr,cols],dev.usefulness.iloc[tr]);cvpred[name][te]=model.predict(dev.loc[te,cols])
   for name,pred in cvpred.items(): results.append({'arm':arm,'descriptor':name,'evaluation':'development_pathway_grouped_cv','n_development_reactions':len(dev),'n_confirmation_reactions':len(conf),'mae':float(np.mean(np.abs(pred-dev.usefulness))),'spearman':float(spearmanr(pred,dev.usefulness).statistic),'partial_rank_eta2_given_entropy_development':partial(dev),'partial_rank_eta2_given_entropy_confirmation':partial(conf)})
   for name,cols in [('entropy',['H_dir']),('entropy_plus_eta2',['H_dir','eta2'])]:
    model=make_pipeline(StandardScaler(),Ridge(alpha=1.0));model.fit(dev[cols],dev.usefulness);pred=model.predict(conf[cols]) if len(conf) else np.array([])
    results.append({'arm':arm,'descriptor':name,'evaluation':'confirmation_from_development_fit','n_development_reactions':len(dev),'n_confirmation_reactions':len(conf),'mae':float(np.mean(np.abs(pred-conf.usefulness))) if len(conf) else np.nan,'spearman':float(spearmanr(pred,conf.usefulness).statistic) if len(conf)>2 else np.nan,'partial_rank_eta2_given_entropy_development':partial(dev),'partial_rank_eta2_given_entropy_confirmation':partial(conf)})
   for label,z in [('development_pathway_grouped_cv',dev),('confirmation_from_development_fit',conf)]:
    for desc,col in [('eta2_association','eta2'),('entropy_association','H_dir')]:
     results.append({'arm':arm,'descriptor':desc,'evaluation':label,'n_reactions':len(z),'spearman':float(spearmanr(z[col],z.usefulness).statistic) if len(z)>2 else np.nan,'partial_rank_eta2_given_entropy_development':partial(dev),'partial_rank_eta2_given_entropy_confirmation':partial(conf)})
 out=pd.DataFrame(results);increments=[]
 for (arm,ev),g in out.loc[out.descriptor.isin(['entropy','entropy_plus_eta2'])].groupby(['arm','evaluation'],sort=True):
  base=g.loc[g.descriptor.eq('entropy')];aug=g.loc[g.descriptor.eq('entropy_plus_eta2')]
  if len(base)==1 and len(aug)==1:
   b=base.iloc[0];a=aug.iloc[0];increments.append({'arm':arm,'descriptor':'eta2_increment','evaluation':ev,'n_development_reactions':b.get('n_development_reactions',np.nan),'n_confirmation_reactions':b.get('n_confirmation_reactions',np.nan),'mae_improvement':float(b.mae-a.mae) if pd.notna(b.get('mae')) and pd.notna(a.get('mae')) else np.nan,'prediction_spearman_change':float(a.spearman-b.spearman) if pd.notna(b.get('spearman')) and pd.notna(a.get('spearman')) else np.nan,'partial_rank_eta2_given_entropy_development':a.partial_rank_eta2_given_entropy_development,'partial_rank_eta2_given_entropy_confirmation':a.partial_rank_eta2_given_entropy_confirmation})
 return pd.concat([out,pd.DataFrame(increments)],ignore_index=True)

def paired_arm_associations(outcomes,geometry,reaction_registry,split,*,target_column=None):
 """Calculate A1/A2 descriptor models on the exact shared finite reaction set."""
 results=[]; groups=['scenario']+([target_column] if target_column else [])
 grouper=groups[0] if len(groups)==1 else groups
 for keys,sub in outcomes.groupby(grouper,sort=True):
  keyvals=(keys,) if not isinstance(keys,tuple) else keys; metadata=dict(zip(groups,keyvals)); shared=None; arms={}
  for arm in ('A1','A2'):
   od=sub.loc[sub.arm.eq(arm)];gm=geometry.loc[geometry.arm.eq(arm)]
   joined=od[['reaction_id','usefulness']].merge(gm[['reaction_id','eta2','H_dir']],on='reaction_id',validate='one_to_one').dropna()
   arms[arm]=(od,gm); ids=set(joined.reaction_id.astype(str));shared=ids if shared is None else shared&ids
  for arm,(od,gm) in arms.items():
   if not shared:continue
   aa=associations(od.loc[od.reaction_id.astype(str).isin(shared)],gm.loc[gm.reaction_id.astype(str).isin(shared)],gm, reaction_registry,split)
   for k,v in metadata.items():
    if k not in aa.columns: aa.insert(0 if k=='scenario' else 1,k,v)
   aa['n_shared_a1_a2_reactions']=len(shared);results.extend(aa.to_dict('records'))
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
   rows,supp,_=score_frozen_truth_pairs(data,wf,common,scenario=scenario,arm=arm,target=20.,outdir=out,truth_pairs=pp)
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
 for scenario,sub in outdf.groupby('scenario',sort=True):
  g=pd.concat([g.assign(arm=arm) for arm,g in geom],ignore_index=True)
  assoc=paired_arm_associations(sub,g,rxreg,split)
  if len(assoc):assoc.insert(0,'scenario',scenario);metrics.extend(assoc.to_dict('records'))
 outdf.to_csv(out/'experiment3b_usefulness.tsv.gz',sep='\t',index=False,compression='gzip')
 pd.DataFrame(alls).to_csv(out/'experiment3b_pair_support.tsv.gz',sep='\t',index=False,compression='gzip')
 pairs.to_csv(out/'experiment3b_pair_eligibility.tsv.gz',sep='\t',index=False,compression='gzip')
 selections.to_csv(out/'experiment3b_selected_candidates.tsv.gz',sep='\t',index=False,compression='gzip')
 duplicates.to_csv(out/'experiment3b_duplicates.tsv',sep='\t',index=False)
 pd.DataFrame(metrics).to_csv(out/'experiment3b_associations.tsv',sep='\t',index=False)
 scenario_denominators={str(s):{'selection_rows':int(len(g)),'evaluable_selection_rows':int(g.pair_evaluable.eq('true').sum()),'distinct_evaluable_truth_pairs':int(g.loc[g.pair_evaluable.eq('true')].drop_duplicates(['evaluation_id','ct2a_candidate_id','gl261_candidate_id']).shape[0])} for s,g in pairs.groupby('scenario',sort=True)}
 (out/'experiment3b_manifest.json').write_text(json.dumps({'status':'COMPLETED','quantiles':[.20,.50,.80],'random_seeds':list(TRUTH_SEEDS),'random_sampling':'weighted categorical with replacement; three selections per pool','paired_arms':'same selected candidate IDs; scenarios use the common support intersection','scenario_denominators':scenario_denominators,'lambda':LAM,'weight_calibration':'frozen original A1 and A2 ESS20 weights; no recalibration','source_hashes':EXPECTED|{A1_CURRENT:A1_CURRENT_SHA},'source_root':str(root),'code_sha256':{p:sha256(root/p) for p in CODE_PATHS},'analysis_code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'archived_weight_kernel_sha256':ARCHIVED_KERNEL_SHA},indent=2)+'\n')

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
 fixed_common=fixed_common_population(common_by_target)
 all_case=[]; all_support=[]
 for target in (10.,20.,40.):
  for arm in ('A1','A2'):
   rows,supp,_=score_frozen_truth_pairs(data,weights_for,common_by_target[target],scenario='ESS_'+str(int(target))+'_eligible',arm=arm,target=target,outdir=out)
   all_case.extend(rows); all_support.extend(supp)
   rows,supp,_=score_frozen_truth_pairs(data,weights_for,fixed_common,scenario='ESS_'+str(int(target))+'_fixed_common',arm=arm,target=target,outdir=out)
   all_case.extend(rows); all_support.extend(supp)
 pd.DataFrame(all_case).to_csv(out/'experiment4_gain_usefulness.tsv.gz',sep='\t',index=False,compression='gzip')
 pd.DataFrame(all_support).to_csv(out/'experiment4_pair_support.tsv.gz',sep='\t',index=False,compression='gzip')
 if all_case:
  cf=pd.DataFrame(all_case)
  paired=cf.pivot(index=['scenario','target_ess','reaction_id'],columns='arm',values=['usefulness','correct_gain_mean','wrong_gain_mean','correct_minus_wrong_mean'])
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
 casedf=pd.DataFrame(all_case);assoc=[]
 for scenario in sorted(casedf.scenario.unique()):
  for target in (10.,20.,40.):
   od=casedf.loc[casedf.scenario.eq(scenario)&casedf.target_ess.eq(target)]
   gm=geom.loc[geom.target_ess.eq(target)]
   if len(od) and len(gm):
    aa=paired_arm_associations(od.assign(scenario=scenario),gm,rxreg,split)
    aa.insert(1,'target_ess',target);assoc.extend(aa.to_dict('records'))
 pd.DataFrame(assoc).to_csv(out/'experiment4_eta2_usefulness_associations.tsv',sep='\t',index=False)
 manifest={'status':'ESS20_PARITY_PASS_OUTCOME_RECOMPUTED','targets':[10,20,40],'lambda':LAM,'truth_registry_sha256':EXPECTED['outputs/dmi_bridge_pl2a_sign_only_prepare_v1/BRIDGEPL2A_TRUTH_PAIR_REGISTRY.tsv'],'distinct_original_truth_pairs':int(len(truth)),'paired_eligible_truth_pairs_by_target':{str(int(t)):len(k) for t,k in common_by_target.items()},'fixed_common_truth_pairs':len(fixed_common),'fixed_common_identity_sha256':hashlib.sha256('\n'.join(sorted('|'.join(k) for k in fixed_common)).encode()).hexdigest(),'source_hashes':EXPECTED|{A1_CURRENT:A1_CURRENT_SHA,A1_MANIFEST:A1_MANIFEST_SHA,ARCHIVED_KERNEL:ARCHIVED_KERNEL_SHA},'code_sha256':{p:sha256(root/p) for p in CODE_PATHS},'analysis_code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'calibration':{'operator':'dmi_bridge_a20_dual_anchor_core_v1.solve_temperature_for_ess','global_pool_n':320,'targets':[10,20,40],'low':1e-12,'initial_high':1.0,'max_high':1e6,'geometric_bisection_iterations':80},'source_root':str(root),'note':'A1-A2 contrasts are calibration sensitivities, not isolated effects of adding a measurement.'}
 (out/'experiment4_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')

def finalize3a(data,out):
 """Complete association outputs from a scored checkpoint after parity validation."""
 root=data['root'];out=Path(out)
 parity=pd.read_csv(out/'experiment3a_20x20_parity.tsv.gz',sep='\t'); maxdiff=float(parity.abs_difference.max())
 if not math.isfinite(maxdiff) or maxdiff>1e-10: raise RuntimeError('frozen geometry parity failed; refusing to finalize')
 geometry=pd.read_csv(out/'experiment3a_geometry_comparison.tsv.gz',sep='\t');outdf=pd.read_csv(out/'experiment3a_usefulness.tsv.gz',sep='\t')
 rxreg=pd.read_csv(root/'outputs/dmi_bridge_pl2a_sign_only_prepare_v1/BRIDGEPL2A_REACTION_REGISTRY.tsv',sep='\t')
 split=pd.read_csv(root/'outputs/dmi_bridge_pl2c_geometry_utility_development_v1/BRIDGEPL2C_REACTION_SPLIT.tsv',sep='\t');assoc=[]
 for arm in ('A1','A2'):
  ag=geometry.loc[geometry.arm.eq(arm)];common=set(ag.loc[ag.eta2_19x19.notna()&ag.eta2_20x20.notna(),'reaction_id']);oo=outdf.loc[outdf.arm.eq(arm)&outdf.reaction_id.isin(common)]
  for spec,metrics in [('original_20x20',('eta2_20x20','H_dir_20x20','D_dir_20x20')),('truth_excluded_19x19',('eta2_19x19','H_dir_19x19','D_dir_19x19'))]:
   gm=ag[['arm','reaction_id',*metrics]].rename(columns=dict(zip(metrics,('eta2','H_dir','D_dir'))));result=associations(oo,gm,gm,rxreg,split);result=result.loc[result.arm.eq(arm)].copy();result.insert(1,'geometry_spec',spec);assoc.append(result)
 pd.concat(assoc,ignore_index=True).to_csv(out/'experiment3a_eta2_usefulness.tsv',sep='\t',index=False)
 manifest={'status':'20x20_PARITY_PASS','geometry_label':'truth-exclusion sensitivity descriptor; not a pre-cue predictor','full_descriptor_reaction_inventory':4180,'weak_cue_reaction_inventory':4179,'frozen_evaluable_distinct_pairs':836,'max_parity_abs_diff':maxdiff,'source_hashes':EXPECTED|{A1_CURRENT:A1_CURRENT_SHA},'source_root':str(root),'code_sha256':{p:sha256(root/p) for p in CODE_PATHS},'archived_weight_kernel_sha256':ARCHIVED_KERNEL_SHA,'analysis_code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 (out/'experiment3a_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')

def rescore3a(data,out):
 """Rerun weak-cue outcomes from frozen identities while reusing parity-passed geometry."""
 out=Path(out);pairs=data['pairs'].loc[data['pairs'].pair_evaluable.eq('true')].drop_duplicates(['evaluation_id','ct2a_candidate_id','gl261_candidate_id']).copy()
 if len(pairs)!=836:raise RuntimeError('original matched truth-pair denominator changed')
 evals=data['evals'].set_index('evaluation_id').to_dict('index');import dmi_bridge_pl2a_prepare_v1 as pl2a
 candidates={pl2a.candidate_identity(r):r for r in data['panel'].to_dict('records')};keys=[]
 for p in pairs.itertuples(index=False):
  ev=evals[str(p.evaluation_id)];ok=True
  for arm,wt in [('A1',data['a1bio']),('A2',data['a2'])]:
   for mouse,tumor,cid in [(str(p.ct2a_mouse),'CT2A',str(p.ct2a_candidate_id)),(str(p.gl261_mouse),'GL261',str(p.gl261_candidate_id))]:
    c=candidates[cid];w=original_weights(wt,mouse,tumor,str(ev['algorithm']),str(ev['rna_context_key']))
    if w is None or w.sum()-w[int(c['sample_index'])]<=0:ok=False
  if ok:keys.append((str(p.evaluation_id),str(p.ct2a_candidate_id),str(p.gl261_candidate_id)))
 common=set(keys)
 if len(common)!=836:raise RuntimeError(f'matched support changed during 3A rescore: {len(common)}')
 pp=data['pairs'].copy();pp['pair_evaluable']=[('true' if (str(r.evaluation_id),str(r.ct2a_candidate_id),str(r.gl261_candidate_id)) in common else 'false') for r in pp.itertuples(index=False)];data['pairs']=pp
 def wf(mouse,tumor,alg,rna,target,arm):return original_weights(data['a1bio'] if arm=='A1' else data['a2'],mouse,tumor,alg,rna)
 outcomes=[];supports=[]
 for arm in ('A1','A2'):
  rows,supp,_=score_frozen_truth_pairs(data,wf,common,scenario='truth_excluded_q10_q50_q90',arm=arm,target=20.,outdir=out);outcomes.extend(rows);supports.extend(supp)
 pd.DataFrame(outcomes).to_csv(out/'experiment3a_usefulness.tsv.gz',sep='\t',index=False,compression='gzip')
 pd.DataFrame(supports).to_csv(out/'experiment3a_matched_pair_support.tsv.gz',sep='\t',index=False,compression='gzip')
 # Frozen PL2C evaluation-level parity is a hard gate before downstream associations.
 actual=pd.read_csv(out/'truth_excluded_q10_q50_q90_A1_evaluation_reaction.tsv.gz',sep='\t')
 frozen=pd.read_csv(data['root']/'outputs/dmi_bridge_pl2c_geometry_utility_development_v1/BRIDGEPL2C_DEVELOPMENT_EVALUATION_REACTION.tsv.xz',sep='\t',compression='xz',usecols=['evaluation_id','reaction_id','n_distinct_pairs','n_non_tie_pairs','n_truth_tie_pairs','mean_correct_gain','mean_information_advantage','directionally_useful_fraction'])
 joined=actual.merge(frozen,on=['evaluation_id','reaction_id'],suffixes=('_new','_frozen'),validate='one_to_one')
 checks=[]
 for a,b in [('n_distinct_truth_pairs','n_distinct_pairs'),('n_non_tie_truth_pairs','n_non_tie_pairs'),('n_truth_ties','n_truth_tie_pairs'),('correct_gain_mean','mean_correct_gain'),('information_advantage_mean','mean_information_advantage'),('usefulness','directionally_useful_fraction')]:
  d=np.abs(joined[a].to_numpy(float)-joined[b].to_numpy(float));maximum=float(np.nanmax(d))
  checks.append({'field':a,'reference':b,'n_rows':len(joined),'max_abs_difference':maximum,'tolerance':ATOL,'status':'PASS' if maximum<=ATOL else 'FAIL'})
 split=pd.read_csv(data['root']/'outputs/dmi_bridge_pl2c_geometry_utility_development_v1/BRIDGEPL2C_REACTION_SPLIT.tsv',sep='\t')
 dev_ids=set(split.loc[split.analysis_split.eq('DEVELOPMENT'),'reaction_id'].astype(str))-{'LDH_L'}
 dev=joined.loc[joined.reaction_id.astype(str).isin(dev_ids)]
 reaction_new=dev.groupby('reaction_id',sort=True)[['usefulness','correct_gain_mean','information_advantage_mean']].agg(['mean','count'])
 reaction_ref=frozen.loc[frozen.reaction_id.astype(str).isin(dev_ids)].groupby('reaction_id',sort=True)[['directionally_useful_fraction','mean_correct_gain','mean_information_advantage']].agg(['mean','count'])
 for a,b in [('usefulness','directionally_useful_fraction'),('correct_gain_mean','mean_correct_gain'),('information_advantage_mean','mean_information_advantage')]:
  for suffix in ('mean','count'):
   d=np.abs(reaction_new[(a,suffix)].to_numpy(float)-reaction_ref[(b,suffix)].to_numpy(float));maximum=float(np.nanmax(d))
   checks.append({'field':f'reaction_{a}_{suffix}','reference':f'PL2C_{b}_{suffix}','n_rows':len(d),'max_abs_difference':maximum,'tolerance':ATOL,'status':'PASS' if maximum<=ATOL else 'FAIL'})
 if any(x['status']!='PASS' for x in checks):
  pd.DataFrame(checks).to_csv(out/'experiment3a_frozen_aggregation_parity.tsv',sep='\t',index=False)
  raise RuntimeError('frozen PL2C evaluation-level numerical parity failed; stop affected analyses')
 pd.DataFrame(checks).to_csv(out/'experiment3a_frozen_aggregation_parity.tsv',sep='\t',index=False)

def finalize4(data,out):
 """Finalize paired A1/A2 associations and diagnostics from completed ESS tables."""
 root=data['root'];out=Path(out)
 casedf=pd.read_csv(out/'experiment4_gain_usefulness.tsv.gz',sep='\t');geom=pd.read_csv(out/'experiment4_geometry_summary.tsv.gz',sep='\t')
 paired=casedf.pivot(index=['scenario','target_ess','reaction_id'],columns='arm',values=['usefulness','correct_gain_mean','wrong_gain_mean','correct_minus_wrong_mean'])
 paired.columns=['_'.join(c) for c in paired.columns];paired.reset_index().to_csv(out/'experiment4_a1_a2_paired.tsv.gz',sep='\t',index=False,compression='gzip')
 rxreg=pd.read_csv(root/'outputs/dmi_bridge_pl2a_sign_only_prepare_v1/BRIDGEPL2A_REACTION_REGISTRY.tsv',sep='\t')
 split=pd.read_csv(root/'outputs/dmi_bridge_pl2c_geometry_utility_development_v1/BRIDGEPL2C_REACTION_SPLIT.tsv',sep='\t');assocs=[]
 for scenario in sorted(casedf.scenario.unique()):
  for target in (10.,20.,40.):
   od=casedf.loc[casedf.scenario.eq(scenario)&casedf.target_ess.eq(target)];gm=geom.loc[geom.target_ess.eq(target)]
   if len(od) and len(gm):
    aa=paired_arm_associations(od,gm,rxreg,split)
    if len(aa):aa.insert(1,'target_ess',target);assocs.extend(aa.to_dict('records'))
 pd.DataFrame(assocs).to_csv(out/'experiment4_eta2_usefulness_associations.tsv',sep='\t',index=False)
 support=pd.read_csv(out/'experiment4_pair_support.tsv.gz',sep='\t')
 weights=pd.read_csv(out/'experiment4_weight_support.tsv.gz',sep='\t');product=pd.read_csv(out/'experiment4_product_ess.tsv.gz',sep='\t')
 summary=[]
 for (scenario,target,arm),g in casedf.groupby(['scenario','target_ess','arm'],sort=True):
  row={'scenario':scenario,'target_ess':target,'arm':arm,'n_finite_reactions':int(g.usefulness.notna().sum()),'mean_usefulness_fraction_across_reactions':float(g.usefulness.mean()),'mean_correct_gain_across_reactions':float(g.correct_gain_mean.mean()),'mean_wrong_gain_across_reactions':float(g.wrong_gain_mean.mean()),'mean_correct_minus_wrong_across_reactions':float(g.correct_minus_wrong_mean.mean())}
  ss=support.loc[support.scenario.eq(scenario)&support.arm.eq(arm)]
  row['n_supported_truth_pairs']=int(ss.loc[ss.eligible].drop_duplicates(['evaluation_id','ct2a_candidate_id','gl261_candidate_id']).shape[0]);row['n_supported_evaluations']=int(ss.loc[ss.eligible,'evaluation_id'].nunique())
  w=weights.loc[weights.arm.eq(arm)&weights.target_global_ess.eq(target)];row['mean_conditional_ess']=float(w.conditional_ess.mean());row['min_conditional_ess']=float(w.conditional_ess.min());row['n_conditional_support_pools']=int(w.conditional_support.sum())
  pr=product.loc[product.arm.eq(arm)&product.target_global_ess.eq(target)];row['mean_product_weight_ess']=float(pr.product_weight_ess.mean());summary.append(row)
 pd.DataFrame(summary).to_csv(out/'experiment4_endpoint_summary.tsv',sep='\t',index=False)
 # ESS20 must reproduce the independently completed original 3A outcome table.
 baseline=[];threea=Path(data['root']).parent/'recomb2027-dmi/reproduced/reviewer_checks/corrected/experiment_3/experiment_3a'
 for arm in ('A1','A2'):
  got=pd.read_csv(out/f'ESS_20_eligible_{arm}_evaluation_reaction.tsv.gz',sep='\t')
  ref=pd.read_csv(threea/f'truth_excluded_q10_q50_q90_{arm}_evaluation_reaction.tsv.gz',sep='\t')
  j=got.merge(ref,on=['evaluation_id','reaction_id'],suffixes=('_ess20','_frozen'),validate='one_to_one')
  for a,b in [('n_distinct_truth_pairs','n_distinct_truth_pairs'),('n_non_tie_truth_pairs','n_non_tie_truth_pairs'),('n_truth_ties','n_truth_ties'),('usefulness','usefulness'),('correct_gain_mean','correct_gain_mean'),('wrong_gain_mean','wrong_gain_mean'),('information_advantage_mean','information_advantage_mean')]:
   d=np.abs(j[f'{a}_ess20'].to_numpy(float)-j[f'{b}_frozen'].to_numpy(float));mx=float(np.nanmax(d))
   baseline.append({'arm':arm,'field':a,'n_rows':len(j),'max_abs_difference':mx,'tolerance':ATOL,'status':'PASS' if mx<=ATOL else 'FAIL'})
 pd.DataFrame(baseline).to_csv(out/'experiment4_ess20_outcome_parity.tsv',sep='\t',index=False)
 if any(x['status']!='PASS' for x in baseline):raise RuntimeError('ESS20 outcome parity with the frozen original population failed; stopping Experiment 4 interpretation')
 common=support.loc[support.scenario.str.endswith('_fixed_common')&support.eligible].drop_duplicates(['evaluation_id','ct2a_candidate_id','gl261_candidate_id'])
 common_keys=set(zip(common.evaluation_id.astype(str),common.ct2a_candidate_id.astype(str),common.gl261_candidate_id.astype(str)))
 target_counts={}
 for target in (10,20,40):
  q=support.loc[support.scenario.eq(f'ESS_{target}_eligible')&support.eligible]
  target_counts[str(target)]=int(q.drop_duplicates(['evaluation_id','ct2a_candidate_id','gl261_candidate_id']).shape[0])
 manifest={'status':'ESS20_WEIGHT_AND_OUTCOME_PARITY_PASS','targets':[10,20,40],'lambda':LAM,'truth_registry_sha256':EXPECTED['outputs/dmi_bridge_pl2a_sign_only_prepare_v1/BRIDGEPL2A_TRUTH_PAIR_REGISTRY.tsv'],'distinct_original_truth_pairs':int(len(data['pairs'].loc[data['pairs'].pair_evaluable.eq('true')].drop_duplicates(['evaluation_id','ct2a_candidate_id','gl261_candidate_id']))),'paired_eligible_truth_pairs_by_target':target_counts,'fixed_common_truth_pairs':len(common_keys),'fixed_common_identity_sha256':hashlib.sha256('\n'.join('|'.join(x) for x in sorted(common_keys)).encode()).hexdigest(),'ess20_outcome_parity':'PASS within 1e-12 against corrected frozen 3A outcomes','paired_association_population':'intersection of finite reaction IDs across both arms for each scenario and ESS target','source_hashes':EXPECTED|{A1_CURRENT:A1_CURRENT_SHA,A1_MANIFEST:A1_MANIFEST_SHA,ARCHIVED_KERNEL:ARCHIVED_KERNEL_SHA},'source_root':str(root),'code_sha256':{p:sha256(root/p) for p in CODE_PATHS},'analysis_code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'calibration':{'operator':'dmi_bridge_a20_dual_anchor_core_v1.solve_temperature_for_ess','global_pool_n':320,'targets':[10,20,40],'low':1e-12,'initial_high':1.0,'max_high':1e6,'geometric_bisection_iterations':80},'note':'A1-A2 contrasts are calibration sensitivities, not isolated effects of adding a measurement.'}
 (out/'experiment4_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')

def rescore_a1_ess20_from_frozen(data,out):
 """Use the exact PL2 A1 baseline weights at ESS20 to avoid ulp-scale endpoint flips."""
 out=Path(out);support=pd.read_csv(out/'experiment4_pair_support.tsv.gz',sep='\t')
 selected={}
 for scenario in ('ESS_20_eligible','ESS_20_fixed_common'):
  s=support.loc[support.scenario.eq(scenario)&support.arm.eq('A1')&support.eligible]
  selected[scenario]=set(zip(s.evaluation_id.astype(str),s.ct2a_candidate_id.astype(str),s.gl261_candidate_id.astype(str)))
 def wf(mouse,tumor,alg,rna,target,arm):
  if arm!='A1':raise RuntimeError('A1 ESS20 rescore must not alter the paired A2 arm')
  return original_weights(data['a1bio'],mouse,tumor,alg,rna)
 all_rows=[];all_support=[]
 for scenario,keys in selected.items():
  rows,supp,_=score_frozen_truth_pairs(data,wf,keys,scenario=scenario,arm='A1',target=20.,outdir=out)
  all_rows.extend(rows);all_support.extend(supp)
 old=pd.read_csv(out/'experiment4_gain_usefulness.tsv.gz',sep='\t')
 replace={'ESS_20_eligible','ESS_20_fixed_common'}
 old=old.loc[~(old.scenario.isin(replace)&old.target_ess.eq(20.)&old.arm.eq('A1'))]
 current=pd.concat([old,pd.DataFrame(all_rows)],ignore_index=True);current.to_csv(out/'experiment4_gain_usefulness.tsv.gz',sep='\t',index=False,compression='gzip')
 oldsup=pd.read_csv(out/'experiment4_pair_support.tsv.gz',sep='\t')
 oldsup=oldsup.loc[~(oldsup.scenario.isin(replace)&oldsup.arm.eq('A1'))]
 pd.concat([oldsup,pd.DataFrame(all_support)],ignore_index=True).to_csv(out/'experiment4_pair_support.tsv.gz',sep='\t',index=False,compression='gzip')
 # Record the sub-ulp agreement between calibration weights and the exact PL2 A1 table.
 a1cur=pd.read_csv(data['root']/A1_CURRENT,sep='\t',compression='xz',float_precision='round_trip');key=['mouse_id','algorithm','tumor','ensemble_hash','rna_context_key','sample_index']
 b=data['a1bio'][key+['weight']];j=a1cur[key+['weight']].merge(b,on=key,suffixes=('_calibrated','_PL2_frozen'),validate='one_to_one')
 rows=[]
 for mouse,g in j.groupby('mouse_id',sort=True):
  ess_cal=1/np.dot(g.weight_calibrated,g.weight_calibrated);ess_frozen=1/np.dot(g.weight_PL2_frozen,g.weight_PL2_frozen)
  rows.append({'mouse_id':mouse,'max_abs_weight_difference':float(np.max(np.abs(g.weight_calibrated-g.weight_PL2_frozen))),'calibrated_ess':ess_cal,'PL2_frozen_ess':ess_frozen,'status':'PASS' if np.max(np.abs(g.weight_calibrated-g.weight_PL2_frozen))<=1e-12 and abs(ess_frozen-20.)<=1e-10 else 'FAIL'})
 pd.DataFrame(rows).to_csv(out/'experiment4_a1_ess20_frozen_weight_parity.tsv',sep='\t',index=False)
 if any(x['status']!='PASS' for x in rows):raise RuntimeError('frozen A1 ESS20 weight parity failed')

def main():
 ap=argparse.ArgumentParser(); ap.add_argument('experiment',choices=['3a','3b','4']); ap.add_argument('--upstream-root',type=Path,default=Path('/home/pty/work/project_MRI_brain_tumor')); ap.add_argument('--output-root',type=Path,required=True); ap.add_argument('--pilot',action='store_true'); ap.add_argument('--finalize-existing',action='store_true'); ap.add_argument('--rescore-existing',action='store_true'); ap.add_argument('--restore-frozen-baseline',action='store_true'); ap.add_argument('--restore-ess20-frozen',action='store_true'); ap.add_argument('--finalize-4-existing',action='store_true'); ap.add_argument('--rescore-a1-ess20',action='store_true'); a=ap.parse_args()
 sys.path.insert(0,str(a.upstream_root/'scripts'))
 data=load(a.upstream_root); data['root']=a.upstream_root
 # Source core modules imported by load; bind explicitly for numerical calls.
 import dmi_bridge_pl1_predictability_core_v1 as pl1
 import dmi_bridge_pl2_sign_only_core_v1 as pl2
 import dmi_bridge_a20_dual_anchor_core_v1 as a20
 data.update(pl1=pl1,pl2=pl2,a20mod=a20)
 if a.experiment=='4' and a.restore_ess20_frozen: restore4_ess20_frozen(data,a.output_root); finalize4(data,a.output_root)
 elif a.experiment=='4' and a.rescore_a1_ess20: rescore_a1_ess20_from_frozen(data,a.output_root)
 elif a.experiment=='4' and a.finalize_4_existing: finalize4(data,a.output_root)
 elif a.experiment=='4': run4(data,a.output_root)
 elif a.experiment=='3a' and a.restore_frozen_baseline: restore3a_frozen_outcomes(data,a.output_root); finalize3a(data,a.output_root)
 elif a.experiment=='3b' and a.restore_frozen_baseline: restore3b_frozen_baseline(data,a.output_root)
 elif a.experiment=='3b' and a.finalize_existing: finalize3b_associations(data,a.output_root)
 elif a.experiment=='3a' and a.rescore_existing: rescore3a(data,a.output_root)
 elif a.experiment=='3a' and a.finalize_existing: finalize3a(data,a.output_root)
 elif a.experiment=='3a': run3a(data,a.output_root,pilot=a.pilot)
 elif a.experiment=='3b': run3b(data,a.output_root)
 else: raise RuntimeError(f'Experiment {a.experiment} runner is not yet staged')

if __name__=='__main__': main()
