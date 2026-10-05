#!/usr/bin/env python3
"""Freeze outcome-blind A2-L residual geometry on the frozen current panel."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import lzma
import math
import os
from pathlib import Path
import shutil

import numpy as np

from scripts import dmi_bridge_pl1_predictability_landscape_v1 as pl1
from scripts import dmi_bridge_pl1_storage_v1 as storage
from scripts import dmi_bridge_pl1_predictability_resume_v1 as resume
from scripts.dmi_bridge_pl1_predictability_batch_v1 import target_landscape_metrics_block
from scripts.dmi_bridge_pl1_predictability_core_v1 import normalize_weights
from scripts import dmi_bridge_a21_dual_anchor_geometry_core_v1 as core

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'outputs/dmi_bridge_a21_dual_anchor_geometry_v1'
A20 = ROOT / 'outputs/dmi_bridge_a20_dual_anchor_qualification_v1'
A20_MANIFEST_SHA = '01f8004a0bf9d5358d7ce7281e827cc1a1fbe00a6750a5c063a1edad2eedf380'
A20_WEIGHTS_SHA = 'b153e2fde45dca9573c71c548783e704313f056d0980f64af60180ad90ebc1f4'
PANEL_ID = '9c27c8e75aa128a4a4cb15547f23557845f37d35a7323e875af5a1711a527f42'
PL1_MANIFEST_SHA = 'b0f81620ff24f5791dd964a8063416784db86a358b0464b8c68a504e84e5f529'
PL1_PARTS_SHA = '97c70016f53ca8497da5d7c34865f90e5c9dd3dbd32ffb531252492ba024cb21'
CONTEXTS, EVALS, REACTIONS = 16, 400, 4179
ROWS_PER_CONTEXT = 25 * REACTIONS
PROXIMAL = {'Glycolysis/gluconeogenesis', 'Pyruvate metabolism', 'Citric acid cycle', 'Glutamate metabolism'}
EXTRA_METRICS = ('delta_hex1_mean', 'delta_hex1_sd', 'delta_hex1_target_correlation', 'abs_delta_hex1_target_correlation', 'delta_lactate_mean', 'delta_lactate_sd', 'delta_lactate_target_correlation', 'abs_delta_lactate_target_correlation')
EXTRA_STATUSES = ('hex1_delta_degenerate', 'hex1_target_coupling_status', 'lactate_delta_degenerate', 'lactate_target_coupling_status')
SUMMARY_METRICS = pl1.REACTION_STATS + ('delta_lactate_mean', 'delta_lactate_sd', 'delta_lactate_target_correlation', 'abs_delta_lactate_target_correlation')


def require(ok: bool, message: str) -> None:
    if not ok:
        raise RuntimeError(message)


def json_write(path: Path, value: dict) -> None:
    tmp = path.with_name(path.name + '.tmp')
    tmp.write_text(json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + '\n')
    os.replace(tmp, path)


def safe_hash(relative: str, expected: str) -> str:
    path = ROOT / relative
    require(path.is_file() and not path.is_symlink() and path.resolve().is_relative_to(ROOT), f'unsafe or missing source: {relative}')
    got = pl1.sha256(path)
    require(got == expected, f'source SHA256 mismatch: {relative}')
    return got


def source_gate() -> dict[str, str]:
    sources = dict(pl1.source_hashes())
    mrel = 'outputs/dmi_bridge_a20_dual_anchor_qualification_v1/BRIDGEA20_MANIFEST.json'
    wrel = 'outputs/dmi_bridge_a20_dual_anchor_qualification_v1/BRIDGEA20_DUAL_ANCHOR_WEIGHTS.tsv.xz'
    sources[mrel] = safe_hash(mrel, A20_MANIFEST_SHA)
    sources[wrel] = safe_hash(wrel, A20_WEIGHTS_SHA)
    manifest = json.loads((ROOT / mrel).read_text())
    require(manifest['status'] == 'BRIDGE_A20_PARTIAL_QUALIFICATION' and manifest['panel_id'] == PANEL_ID, 'A2.0 admission failed')
    require(manifest['qualification'] == {'A2-G': 'NOT_QUALIFIED', 'A2-L': 'QUALIFIED'}, 'A2 arm status mismatch')
    require(manifest['scientific_outcomes_computed'] is False and manifest['a2_g_outcome_accessed'] is False, 'A2 outcome flag mismatch')
    require(manifest['row_counts']['A2-L_weights'] == 3200 and manifest['row_counts']['A2-G_weights'] == 0 and manifest['row_counts']['weight_summary'] == 10, 'A2 row accounting mismatch')
    for name, digest in manifest['artifact_sha256'].items():
        relative = f'outputs/dmi_bridge_a20_dual_anchor_qualification_v1/{name}'
        sources[relative] = safe_hash(relative, digest)
    require(len(manifest['source_sha256']) == 19, 'A2 source audit cardinality mismatch')
    for rec in manifest['source_sha256'].values():
        relative, digest = rec['relative_path'], rec['sha256']
        sources[relative] = safe_hash(relative, digest)
    sources['outputs/dmi_bridge_pl1_predictability_landscape_v1_resumable_v1/BRIDGEPL1_MANIFEST.json'] = safe_hash('outputs/dmi_bridge_pl1_predictability_landscape_v1_resumable_v1/BRIDGEPL1_MANIFEST.json', PL1_MANIFEST_SHA)
    sources['outputs/dmi_bridge_pl1_predictability_landscape_v1_resumable_v1/BRIDGEPL1_REACTION_EVALUATION_FEATURES.parts.json'] = safe_hash('outputs/dmi_bridge_pl1_predictability_landscape_v1_resumable_v1/BRIDGEPL1_REACTION_EVALUATION_FEATURES.parts.json', PL1_PARTS_SHA)
    parts = json.loads((ROOT / 'outputs/dmi_bridge_pl1_predictability_landscape_v1_resumable_v1/BRIDGEPL1_REACTION_EVALUATION_FEATURES.parts.json').read_text())
    require(parts['logical_content_sha256'] == '46aa9f0410f89d03b020e01fefa517348faae1d55d9bfe8387ac04b39bb9ac93' and parts['total_data_rows'] == 1_672_000, 'PL1 logical content identity mismatch')
    for part in parts['parts']:
        relative = 'outputs/dmi_bridge_pl1_predictability_landscape_v1_resumable_v1/' + part['filename']
        sources[relative] = safe_hash(relative, part['sha256'])
    return sources


def code_hashes() -> dict[str, str]:
    paths = [Path(__file__), Path(core.__file__), Path(pl1.__file__), Path(resume.__file__), ROOT / 'scripts/dmi_bridge_pl1_storage_v1.py', ROOT / 'scripts/dmi_bridge_pl1_predictability_core_v1.py', ROOT / 'scripts/dmi_bridge_pl1_predictability_batch_v1.py']
    return {str(p.relative_to(ROOT)): pl1.sha256(p) for p in paths}


def load_weights(by_key: dict) -> dict:
    rows = pl1.read_tsv(A20 / 'BRIDGEA20_DUAL_ANCHOR_WEIGHTS.tsv.xz')
    require(len(rows) == 3200 and all(r['arm'] == 'A2-L' for r in rows), 'A2-L weight inventory mismatch')
    candidate_keys = set()
    for group in by_key.values():
        for r in group:
            candidate_keys.add(tuple(r[k] for k in ('algorithm', 'tumor', 'ensemble_hash', 'projection_hash', 'rna_context_key', 'sample_index')))
    wm: dict[tuple, dict[int, float]] = {}
    for r in rows:
        key = tuple(r[k] for k in ('algorithm', 'tumor', 'ensemble_hash', 'projection_hash', 'rna_context_key', 'sample_index'))
        require(key in candidate_keys, 'A2-L weight/candidate key mismatch')
        group = (r['mouse_id'], r['tumor'], r['algorithm'], r['ensemble_hash'], r['rna_context_key'])
        sample = int(r['sample_index']); weight = float(r['weight'])
        require(math.isfinite(weight) and weight > 0, 'nonpositive A2-L weight')
        wm.setdefault(group, {})
        require(sample not in wm[group], 'duplicate A2-L weight key')
        wm[group][sample] = weight
    require(len(wm) == 160 and all(set(v) == set(range(20)) for v in wm.values()), 'A2-L strata incomplete')
    return {k: [(i, v[i]) for i in range(20)] for k, v in wm.items()}


def candidate_lactate(by_key: dict, cache_index: dict) -> np.ndarray:
    arr = np.empty((32, 20), dtype=float)
    for key, rows in by_key.items():
        vals = np.array([float(r['lac_production']) for r in rows])
        raw = np.array([float(r['LDH_L']) for r in rows])
        require(np.all(np.isfinite(vals)) and np.allclose(vals, core.lactate_production(raw), atol=1e-14, rtol=0), 'raw panel Vlac observable mismatch')
        arr[cache_index[key]] = vals
    return arr


def layout():
    original, status, _ = resume._metric_layout()
    metrics = original + list(EXTRA_METRICS)
    statuses = status + list(EXTRA_STATUSES)
    return metrics, statuses, {m: i for i, m in enumerate(metrics)}


def feature_fields(metrics, statuses):
    return ['evaluation_id', 'algorithm', 'rna_context_key', 'ct2a_mouse', 'gl261_mouse', 'reaction_id', 'subsystem', 'same_subsystem_as_HEX1', 'same_subsystem_as_LDH_L', 'proximal_set_member'] + metrics + statuses


def parity(batch, pos, h_c, h_g, l_c, l_g, b_c, b_g, wc, wg):
    scalar = core.target_landscape_metrics_a21(hex1_ct2a=h_c, hex1_gl261=h_g, lactate_ct2a=l_c, lactate_gl261=l_g, target_ct2a=b_c, target_gl261=b_g, strong_weights_ct2a=wc, strong_weights_gl261=wg)
    for key, expected in scalar.items():
        got = batch[key][pos]
        if isinstance(expected, str):
            require(str(got) == expected, f'scalar status mismatch: {key}')
        elif isinstance(expected, bool):
            require(bool(got) == expected, f'scalar boolean mismatch: {key}')
        elif math.isnan(float(expected)):
            require(math.isnan(float(got)), f'scalar NaN mismatch: {key}')
        else:
            require(np.isclose(float(got), float(expected), rtol=5e-12, atol=5e-12), f'scalar numeric mismatch: {key}')


def context_compute(spec, wm, mats, lactate, cols, hex_idx, metrics, statuses, metric_pos, values, flags, states):
    ci, gi = spec['ct2a_cache_index'], spec['gl261_cache_index']
    hc, hg = mats[ci, :, hex_idx], mats[gi, :, hex_idx]
    lc, lg = lactate[ci], lactate[gi]
    e = spec['eval_start']
    for cmouse in spec['ct2a_mice']:
        wc = normalize_weights([w for _, w in wm[(cmouse, 'CT2A', spec['method'], spec['ct2a_ensemble'], spec['rna_context_key'])]])
        for gmouse in spec['gl261_mice']:
            wg = normalize_weights([w for _, w in wm[(gmouse, 'GL261', spec['method'], spec['gl261_ensemble'], spec['rna_context_key'])]])
            for start in range(0, REACTIONS, pl1.REACTION_BLOCK_SIZE):
                stop = min(start + pl1.REACTION_BLOCK_SIZE, REACTIONS)
                bc, bg = mats[ci][:, cols[start:stop]], mats[gi][:, cols[start:stop]]
                block = target_landscape_metrics_block(anchor_ct2a=hc, anchor_gl261=hg, target_ct2a=bc, target_gl261=bg, strong_weights_ct2a=wc, strong_weights_gl261=wg)
                lb = target_landscape_metrics_block(anchor_ct2a=lc, anchor_gl261=lg, target_ct2a=bc, target_gl261=bg, strong_weights_ct2a=wc, strong_weights_gl261=wg)
                block.update(delta_hex1_mean=block['delta_a_mean'], delta_hex1_sd=block['delta_a_sd'], delta_hex1_target_correlation=block['delta_anchor_target_correlation'], abs_delta_hex1_target_correlation=block['abs_delta_anchor_target_correlation'], delta_lactate_mean=lb['delta_a_mean'], delta_lactate_sd=lb['delta_a_sd'], delta_lactate_target_correlation=lb['delta_anchor_target_correlation'], abs_delta_lactate_target_correlation=lb['abs_delta_anchor_target_correlation'], hex1_delta_degenerate=block['anchor_delta_degenerate'], hex1_target_coupling_status=block['anchor_target_coupling_status'], lactate_delta_degenerate=lb['anchor_delta_degenerate'], lactate_target_coupling_status=lb['anchor_target_coupling_status'])
                for name in metrics:
                    values[e, start:stop, metric_pos[name]] = block[name]
                for j, name in enumerate(statuses):
                    if name.endswith('_status'):
                        states[e, start:stop, j] = np.asarray(block[name] == 'OK', dtype='u1')
                    else:
                        flags[e, start:stop, j] = np.asarray(block[name], dtype='u1')
                if (e - spec['eval_start']) in (0, 12, 24):
                    for pos in sorted({0, (stop-start)//2, stop-start-1}):
                        parity(block, pos, hc, hg, lc, lg, bc[:, pos], bg[:, pos], wc, wg)
            e += 1
    require(e == spec['eval_stop'], 'context evaluation count mismatch')
    values.flush(); flags.flush(); states.flush()


def write_context_part(path, spec, registry, rxns, panel, fields, metrics, statuses, metric_pos, values, flags, states):
    tmp = path.with_name(path.name + '.tmp')
    count = 0
    hex_sub, lac_sub = panel['HEX1'], panel['LDH_L']
    with lzma.open(tmp, 'wt', preset=0, newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for e in range(spec['eval_start'], spec['eval_stop']):
            reg = registry[e]
            for i, rid in enumerate(rxns):
                sub = panel.get(rid, 'UNMAPPED_PATHWAY')
                row = {k: reg[k] for k in ('evaluation_id', 'algorithm', 'rna_context_key', 'ct2a_mouse', 'gl261_mouse')}
                row.update(reaction_id=rid, subsystem=sub, same_subsystem_as_HEX1=sub == hex_sub, same_subsystem_as_LDH_L=sub == lac_sub, proximal_set_member=sub in PROXIMAL)
                row.update({m: values[e, i, metric_pos[m]] for m in metrics})
                for j, name in enumerate(statuses):
                    row[name] = ('OK' if states[e, i, j] else 'DEGENERATE') if name.endswith('_status') else bool(flags[e, i, j])
                writer.writerow({k: pl1.f(row[k]) for k in fields})
                count += 1
    require(count == ROWS_PER_CONTEXT, 'context feature count mismatch')
    os.replace(tmp, path)
    return count


def summary_arrays(values, eval_slice, metric_pos):
    out = {}
    for metric in SUMMARY_METRICS:
        x = np.asarray(values[eval_slice, :, metric_pos[metric]], dtype=float)
        finite = np.isfinite(x); count = finite.sum(axis=0)
        with np.errstate(invalid='ignore'), np.testing.suppress_warnings() as sup:
            sup.filter(RuntimeWarning)
            q10, med, q90 = np.nanquantile(x, (0.1, 0.5, 0.9), axis=0, method='linear')
        out[metric] = (med, q10, q90, count, x.shape[0] - count)
    return out


def summary_row(base, summaries, i, flags, states, sl):
    row = dict(base)
    for metric, arrays in summaries.items():
        for suffix, arr in zip(('median', 'q10', 'q90', 'finite_count', 'nonfinite_count'), arrays):
            row[f'{metric}_{suffix}'] = arr[i]
    for j, name in ((0, 'target_delta_degenerate'), (1, 'target_abs_magnitude_degenerate'), (2, 'anchor_delta_degenerate'), (5, 'hex1_delta_degenerate'), (7, 'lactate_delta_degenerate')):
        row[f'{name}_count'] = int(np.count_nonzero(flags[sl, i, j]))
    for j, name in ((3, 'sign_magnitude_degenerate'), (4, 'anchor_target_coupling_degenerate'), (6, 'hex1_target_coupling_degenerate'), (8, 'lactate_target_coupling_degenerate')):
        row[f'{name}_count'] = int(np.count_nonzero(states[sl, i, j] == 0))
    return row


def finish(out, work, work_manifest, specs, registry, rxns, panel, fields, metrics, statuses, metric_pos, values, flags, states, sources, codes):
    tmp = out.with_name(out.name + '.publish_tmp')
    require(not tmp.exists(), 'conflicting A2.1 publication staging directory')
    tmp.mkdir()
    try:
        pl1.write_tsv(tmp / 'BRIDGEA21_EVALUATION_REGISTRY.tsv', list(registry[0]), registry)
        logical = hashlib.sha256(); parts = []
        for i in range(CONTEXTS):
            src = work / f'features_context_{i:02d}.tsv.xz'
            name = f'BRIDGEA21_REACTION_EVALUATION_FEATURES.part-{i:03d}.tsv.xz'
            shutil.copyfile(src, tmp / name)
            with lzma.open(src, 'rt', newline='') as fh:
                header = fh.readline()
                require(header.rstrip('\n').split('\t') == fields, 'feature header mismatch')
                if i == 0: logical.update(header.encode())
                count = 0
                for line in fh:
                    logical.update(line.encode()); count += 1
            require(count == ROWS_PER_CONTEXT, 'feature part row mismatch')
            parts.append({'filename': name, 'rows': count, 'bytes': (tmp/name).stat().st_size, 'sha256': pl1.sha256(tmp/name)})
        parts_manifest = {'schema': 'bridge.a21.partitioned_tsv.v1', 'logical_artifact': 'BRIDGEA21_REACTION_EVALUATION_FEATURES', 'format': 'tsv.xz.parts', 'number_of_parts': 16, 'total_data_rows': EVALS*REACTIONS, 'column_names': fields, 'logical_content_sha256': logical.hexdigest(), 'parts': parts}
        json_write(tmp / 'BRIDGEA21_REACTION_EVALUATION_FEATURES.parts.json', parts_manifest)
        summary_cols = [f'{m}_{s}' for m in SUMMARY_METRICS for s in ('median','q10','q90','finite_count','nonfinite_count')]
        count_cols = [f'{s}_count' for s in ('target_delta_degenerate','target_abs_magnitude_degenerate','anchor_delta_degenerate','hex1_delta_degenerate','lactate_delta_degenerate','sign_magnitude_degenerate','anchor_target_coupling_degenerate','hex1_target_coupling_degenerate','lactate_target_coupling_degenerate')]
        cfields = ['algorithm','rna_context_key','reaction_id'] + summary_cols + count_cols
        fh, writer = pl1.open_xz_tsv(tmp / 'BRIDGEA21_CONTEXT_REACTION_SUMMARY.tsv.xz', cfields)
        try:
            for spec in specs:
                sl = slice(spec['eval_start'], spec['eval_stop'])
                summaries = summary_arrays(values, sl, metric_pos)
                for i, rid in enumerate(rxns):
                    row = summary_row({'algorithm':spec['method'],'rna_context_key':spec['rna_context_key'],'reaction_id':rid}, summaries, i, flags, states, sl)
                    writer.writerow({k:pl1.f(row[k]) for k in cfields})
        finally: fh.close()
        summaries = summary_arrays(values, slice(0,EVALS), metric_pos)
        rrows = [summary_row({'reaction_id':rid}, summaries, i, flags, states, slice(0,EVALS)) for i,rid in enumerate(rxns)]
        pl1.write_xz_tsv(tmp / 'BRIDGEA21_REACTION_SUMMARY.tsv.xz', ['reaction_id']+summary_cols+count_cols, rrows)
        pathway = {}
        for row in rrows:
            pathway.setdefault(panel.get(row['reaction_id'],'UNMAPPED_PATHWAY'),[]).append([row[f'{m}_median'] for m in SUMMARY_METRICS])
        prows=[]
        for sub, vectors in sorted(pathway.items()):
            arr=np.asarray(vectors,dtype=float); item={'subsystem':sub,'reaction_count':len(vectors)}
            for j,m in enumerate(SUMMARY_METRICS):
                item.update({f'{m}_{k}':v for k,v in pl1.qstats(arr[:,j].tolist()).items()})
            prows.append(item)
        pl1.write_tsv(tmp/'BRIDGEA21_PATHWAY_SUMMARY.tsv',list(prows[0]),prows)
        sensitivity=[]
        for scope,pred in (('primary',lambda s:True),('exclude_same_subsystem',lambda s:s != panel['HEX1']),('exclude_proximal_set',lambda s:s not in PROXIMAL)):
            subset=[r for r in rrows if pred(panel.get(r['reaction_id'],'UNMAPPED_PATHWAY'))]
            item={'scope':scope,'reaction_count':len(subset),'feature_count':EVALS*len(subset)}
            for m in SUMMARY_METRICS:
                item.update({f'{m}_{k}':v for k,v in pl1.qstats([r[f'{m}_median'] for r in subset]).items()})
            sensitivity.append(item)
        pl1.write_tsv(tmp/'BRIDGEA21_SENSITIVITY_SUMMARY.tsv',list(sensitivity[0]),sensitivity)
        counts={'contexts':CONTEXTS,'evaluation_registry':EVALS,'reaction_evaluation_features':EVALS*REACTIONS,'context_reaction_summary':CONTEXTS*REACTIONS,'reaction_summary':REACTIONS,'pathway_summary':len(prows)}
        require(len(prows)==50,'frozen pathway count mismatch')
        json_write(tmp/'BRIDGEA21_SOURCE_AUDIT.json',{'schema':'bridge.a21.source_audit.v1','source_sha256':sources,'code_sha256':codes,'a2_g_status':'NOT_QUALIFIED','pl1_consumed_logical_content_sha256':'46aa9f0410f89d03b020e01fefa517348faae1d55d9bfe8387ac04b39bb9ac93'})
        json_write(tmp/'BRIDGEA21_ANALYSIS_CONTRACT.json',core.analysis_contract())
        ess=np.asarray(values[:,:,metric_pos['joint_ess']])
        require(np.all(np.isfinite(ess)) and np.all(ess>0) and np.all(ess<=400+1e-9),'conditioned joint ESS invalid')
        qc={'status':'PASS','row_counts':counts,'conditioned_joint_ess_min':float(np.min(ess)),'conditioned_joint_ess_max':float(np.max(ess)),'scalar_parity':'PASS','completed_contexts':16,'a2_g_used':False,'scientific_outcomes_computed':False,'same_anchor_subsystem':panel['HEX1']==panel['LDH_L']}
        require(qc['same_anchor_subsystem'],'anchor subsystem mismatch')
        json_write(tmp/'BRIDGEA21_QC.json',qc)
        require(source_gate()==sources and code_hashes()==codes,'source/code changed during computation')
        hashes={p.name:pl1.sha256(p) for p in tmp.iterdir() if p.is_file()}
        manifest={'schema':core.SCHEMA,'status':core.COMPLETE_STATUS,'source_sha256':sources,'code_sha256':codes,'artifact_sha256':hashes,'row_counts':counts,'reaction_inventory':{'cache_reactions':4181,'primary_reactions':REACTIONS,'hex1_index':2793,'ldh_l_index':2909},'feature_logical_content_sha256':logical.hexdigest(),'a2_g_status':'NOT_QUALIFIED','a2_g_used':False,'scientific_outcomes_computed':False,'work_manifest_sha256':pl1.sha256(work/'WORK_MANIFEST.json')}
        json_write(tmp/'BRIDGEA21_MANIFEST.json',manifest)
        require(out.is_dir() and set(out.iterdir()) == {work}, 'conflicting A2.1 output files')
        for artifact in tmp.iterdir():
            os.replace(artifact, out / artifact.name)
        tmp.rmdir()
        return manifest
    except Exception:
        shutil.rmtree(tmp)
        raise


def build(out:Path=OUT)->dict:
    sources=source_gate(); codes=code_hashes()
    if out.exists():
        p=out/'BRIDGEA21_MANIFEST.json'
        require(p.is_file(),'conflicting incomplete A2.1 output')
        m=json.loads(p.read_text())
        require(m.get('status')==core.COMPLETE_STATUS and m.get('source_sha256')==sources and m.get('code_sha256')==codes,'conflicting A2.1 output identity')
        for name,digest in m['artifact_sha256'].items():
            require(pl1.sha256(out/name)==digest,f'completed artifact hash mismatch: {name}')
        metrics_raw=out/'_work'/'metrics.npy'; metrics_xz=out/'_work'/'metrics.npy.xz'
        if metrics_raw.is_file():
            arr=np.load(metrics_raw,mmap_mode='r',allow_pickle=False)
            require(arr.shape==(EVALS,REACTIONS,len(layout()[0])) and arr.dtype==np.dtype('f8'),'legacy metrics checkpoint metadata mismatch')
            del arr
            storage.compress_file_xz_verified(metrics_raw,metrics_xz,preset=9)
            metrics_raw.unlink()
        return m
    candidates,by_key,_,rxns,mats,cache_index,_=pl1.load_inputs()
    wm=load_weights(by_key); lactate=candidate_lactate(by_key,cache_index)
    require(np.count_nonzero(rxns=='HEX1')==1 and np.count_nonzero(rxns=='LDH_L')==1,'anchor reaction inventory mismatch')
    require(int(np.flatnonzero(rxns=='HEX1')[0])==2793 and int(np.flatnonzero(rxns=='LDH_L')[0])==2909,'anchor cache indices changed')
    primary=sorted(set(rxns.tolist())-{'HEX1','LDH_L'})
    require(len(primary)==REACTIONS,'primary reaction count mismatch')
    colmap={r:i for i,r in enumerate(rxns.tolist())}; cols=[colmap[r] for r in primary]
    panel={r['reaction_id']:r['subsystem'] for r in pl1.read_tsv(ROOT/'outputs/dmi_bridge_bio0_audit_v1/BRIDGEBIO0_REACTION_PANEL.tsv')}
    require(panel['HEX1']==panel['LDH_L']=='Glycolysis/gluconeogenesis','anchor subsystem annotation mismatch')
    specs,registry=resume._context_specs(by_key,wm,cache_index)
    for reg in registry:
        spec = next(s for s in specs if s['method'] == reg['algorithm'] and s['rna_context_key'] == reg['rna_context_key'])
        wc = normalize_weights([w for _, w in wm[(reg['ct2a_mouse'], 'CT2A', spec['method'], spec['ct2a_ensemble'], spec['rna_context_key'])]])
        wg = normalize_weights([w for _, w in wm[(reg['gl261_mouse'], 'GL261', spec['method'], spec['gl261_ensemble'], spec['rna_context_key'])]])
        reg['ct2a_conditioned_ess'] = float(1 / np.dot(wc, wc))
        reg['gl261_conditioned_ess'] = float(1 / np.dot(wg, wg))
        reg['conditioned_joint_ess'] = float(reg['ct2a_conditioned_ess'] * reg['gl261_conditioned_ess'])
    metrics,statuses,metric_pos=layout(); fields=feature_fields(metrics,statuses)
    work=out/'_work'; work.mkdir(parents=True,exist_ok=True)
    manifest_path=work/'WORK_MANIFEST.json'
    identity={'schema':'bridge.a21.context_work.v1','source_sha256':sources,'code_sha256':codes,'metric_fields':metrics,'status_fields':statuses,'completed_contexts':{}}
    if manifest_path.exists():
        work_manifest=json.loads(manifest_path.read_text())
        require({k:work_manifest[k] for k in identity if k!='completed_contexts'}=={k:identity[k] for k in identity if k!='completed_contexts'},'A2.1 work identity mismatch')
    else:
        require(set(work.iterdir())==set(),'unmanifested A2.1 work files')
        work_manifest=identity; json_write(manifest_path,work_manifest)
    shape=(EVALS,REACTIONS,len(metrics)); fshape=(EVALS,REACTIONS,len(statuses))
    metrics_raw=work/'metrics.npy'; metrics_xz=work/'metrics.npy.xz'
    def array(name,shape):
        path=work/name
        if path.exists():
            arr=np.load(path,mmap_mode='r+')
            require(arr.shape==shape,'checkpoint array shape mismatch')
            return arr
        if name=='metrics.npy' and metrics_xz.is_file():
            storage.decompress_file_xz(metrics_xz,path)
            arr=np.load(path,mmap_mode='r+',allow_pickle=False)
            require(arr.shape==shape and arr.dtype==np.dtype('f8'),'compressed checkpoint array metadata mismatch')
            return arr
        require(not work_manifest['completed_contexts'],'completed work missing array')
        return np.lib.format.open_memmap(path,mode='w+',dtype='f8' if name=='metrics.npy' else 'u1',shape=shape)
    values=array('metrics.npy',shape); flags=array('flags.npy',fshape); states=array('states.npy',fshape)
    for spec in specs:
        i=spec['context_index']; shard=work/f'features_context_{i:02d}.tsv.xz'
        rec=work_manifest['completed_contexts'].get(str(i))
        if rec:
            require(shard.is_file() and pl1.sha256(shard)==rec['sha256'] and rec['rows']==ROWS_PER_CONTEXT,'completed context shard conflict')
            continue
        shard.unlink(missing_ok=True)
        context_compute(spec,wm,mats,lactate,cols,2793,metrics,statuses,metric_pos,values,flags,states)
        rows=write_context_part(shard,spec,registry,primary,panel,fields,metrics,statuses,metric_pos,values,flags,states)
        work_manifest['completed_contexts'][str(i)]={'file':shard.name,'rows':rows,'sha256':pl1.sha256(shard)}
        json_write(manifest_path,work_manifest)
        print(f'A2.1 context {i+1}/{CONTEXTS} complete',flush=True)
    require(len(work_manifest['completed_contexts'])==CONTEXTS,'incomplete A2.1 work')
    for i in range(CONTEXTS):
        rec = work_manifest['completed_contexts'][str(i)]
        shard = work / rec['file']
        require(shard.is_file() and pl1.sha256(shard) == rec['sha256'] and rec['rows'] == ROWS_PER_CONTEXT,
                f'completed context {i} changed before publication')
    support = np.asarray(values[:, :, metric_pos['cartesian_support']])
    require(np.all(support == 400), 'finite Cartesian support mismatch')
    result=finish(out,work,work_manifest,specs,registry,primary,panel,fields,metrics,statuses,metric_pos,values,flags,states,sources,codes)
    values.flush()
    values._mmap.close()
    del values
    storage.compress_file_xz_verified(metrics_raw,metrics_xz,preset=9)
    metrics_raw.unlink()
    return result


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=OUT);args=parser.parse_args()
    print(json.dumps(build(args.output),sort_keys=True))


if __name__=='__main__': main()
