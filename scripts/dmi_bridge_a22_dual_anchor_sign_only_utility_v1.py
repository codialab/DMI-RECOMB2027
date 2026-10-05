#!/usr/bin/env python3
"""Matched PL2 sign-only utility on the qualified A2-L weights."""
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
from collections import Counter, defaultdict

import numpy as np

from scripts import dmi_bridge_a22_dual_anchor_sign_only_core_v1 as a22
from scripts import dmi_bridge_pl2_sign_only_core_v1 as scalar
from scripts import dmi_bridge_pl2b_batch_v1 as batch
from scripts import dmi_bridge_pl2b_production_v1 as pl2b
from scripts import dmi_bridge_pl2a_prepare_v1 as pl2a
from scripts import dmi_bridge_pl1_storage_v1 as storage

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'outputs/dmi_bridge_a22_dual_anchor_sign_only_utility_v1'
A20 = Path('outputs/dmi_bridge_a20_dual_anchor_qualification_v1')
A21 = Path('outputs/dmi_bridge_a21_dual_anchor_geometry_v1')
PL2A = Path('outputs/dmi_bridge_pl2a_sign_only_prepare_v1')
PART_COUNT = 128
BLOCK_SIZE = 64
MAX_PART_BYTES = 64 * 1024 * 1024
EXPECTED = {
    str(A20 / 'BRIDGEA20_MANIFEST.json'): '01f8004a0bf9d5358d7ce7281e827cc1a1fbe00a6750a5c063a1edad2eedf380',
    str(A21 / 'BRIDGEA21_MANIFEST.json'): 'db9948d7b99447a695cf8aa21429dde6d642e1fd0632bcce52643deac59e6feb',
    str(PL2A / 'BRIDGEPL2A_MANIFEST.json'): '8ae3eb60eb23614e379e38e3aa3cb72ec16c6bdb69aa9c74cc94499e31f1c312',
    'outputs/dmi_bridge_pl2b_sign_only_utility_v1/BRIDGEPL2B_MANIFEST.json': '7d9ca97a3f4ea5aa8b0b983e09ee83fc3c82f326d1b1cbf9c6137cfadf25539d',
    'dmi_bridge_a22_dual_anchor_sign_only_utility_foundation_v1.patch': '5b732959bb16410dd94f0bffa234cbdf7c49e242114b2d4fc7c588923c6fb2a5',
    'scripts/dmi_bridge_pl2_sign_only_core_v1.py': '9a32da4aa7c58ad7d9bfb29f38f37aa2a05813da4503675c3bf43c9542b0d4ce',
    'scripts/dmi_bridge_pl2b_batch_v1.py': '7da76f5e2d4f2f05a0c8f865014fa9a31954a8abf8ea0a03ea2700cc26b68c7b',
    'docs/DMI_BRIDGE_PL2_SIGN_ONLY_UTILITY.md': '4c6e8b717d3dc78cd0fc9803efff2826587c2a4de3d276eeb4c933e8841acea9',
    'docs/DMI_BRIDGE_PL2B_PRODUCTION_CONTRACT.md': '80e07f0d8081dfbdcda7f8441c7ad612dcbf1314b31664a754f0582dcfa5d84e',
}


def require(ok: bool, message: str) -> None:
    if not ok:
        raise RuntimeError(message)


def finalize_gain_storage(work: Path) -> None:
    raw=work/'summary_gains.dat'; compressed=work/'summary_gains.dat.xz'
    if raw.is_file():
        expected_bytes=3*a22.EXPECTED_CASE_ROWS*2*np.dtype('f8').itemsize
        require(raw.stat().st_size==expected_bytes,'legacy summary gains layout mismatch')
        storage.compress_file_xz_verified(raw,compressed,preset=9)
        raw.unlink()


def safe_file(relative: str) -> Path:
    path = ROOT / relative
    require(path.is_file() and not path.is_symlink() and path.resolve().is_relative_to(ROOT), f'unsafe source: {relative}')
    return path


def checked_hash(relative: str, expected: str) -> str:
    actual = pl2b.sha256_file(safe_file(relative))
    require(actual == expected, f'SHA256 mismatch: {relative}')
    return actual


def verify_bundle(base: Path, manifest: dict) -> None:
    for name, digest in manifest['artifact_sha256'].items():
        checked_hash(str(base / name), digest)


def admission() -> tuple[dict, dict, dict]:
    a22.validate_frozen_pl2_semantics()
    for relative, digest in EXPECTED.items():
        checked_hash(relative, digest)
    m20 = json.loads(safe_file(str(A20 / 'BRIDGEA20_MANIFEST.json')).read_text())
    m21 = json.loads(safe_file(str(A21 / 'BRIDGEA21_MANIFEST.json')).read_text())
    mp = json.loads(safe_file(str(PL2A / 'BRIDGEPL2A_MANIFEST.json')).read_text())
    require(m20['status'] == 'BRIDGE_A20_PARTIAL_QUALIFICATION' and m20['qualification'] == {'A2-L': 'QUALIFIED', 'A2-G': 'NOT_QUALIFIED'}, 'A2.0 qualification mismatch')
    require(m20['a2_g_outcome_accessed'] is False and m20['row_counts']['A2-L_weights'] == 3200 and m20['row_counts']['A2-G_weights'] == 0, 'A2.0 arm accounting mismatch')
    require(m21['status'] == 'BRIDGE_A21_DUAL_ANCHOR_GEOMETRY_FROZEN' and m21['scientific_outcomes_computed'] is False and m21['a2_g_used'] is False, 'A2.1 status mismatch')
    require(m21['feature_logical_content_sha256'] == '50bc3e47239e4b22d0931246f4fa0d5d126a5b69226d47cb342cb8a5b3ee99d2', 'A2.1 logical digest mismatch')
    require(m21['reaction_inventory'] == {'cache_reactions': 4181, 'primary_reactions': 4179, 'hex1_index': 2793, 'ldh_l_index': 2909}, 'A2.1 reaction inventory mismatch')
    require(m21['row_counts']['reaction_evaluation_features'] == 1671600 and m21['row_counts']['context_reaction_summary'] == 66864, 'A2.1 row count mismatch')
    require(mp['status'] == 'PL2A_READY_FOR_PL2B_SIGN_ONLY_UTILITY', 'PL2A status mismatch')
    verify_bundle(A20, m20)
    verify_bundle(A21, m21)
    verify_bundle(PL2A, mp)
    for base, manifest in ((A20, m20), (A21, m21)):
        for name, item in manifest['source_sha256'].items():
            if isinstance(item, dict):
                checked_hash(item['relative_path'], item['sha256'])
            else:
                checked_hash(name, item)
    for name, digest in m21['code_sha256'].items():
        checked_hash(name, digest)
    audit = json.loads(safe_file(str(A21 / 'BRIDGEA21_SOURCE_AUDIT.json')).read_text())
    for name, digest in audit['source_sha256'].items():
        checked_hash(name, digest)
    for name, digest in audit['code_sha256'].items():
        checked_hash(name, digest)
    qc = json.loads(safe_file(str(A21 / 'BRIDGEA21_QC.json')).read_text())
    require(qc['status'] == 'PASS' and qc['scalar_parity'] == 'PASS', 'A2.1 QC mismatch')
    old = json.loads(safe_file('outputs/dmi_bridge_pl2b_sign_only_utility_v1/BRIDGEPL2B_MANIFEST.json').read_text())
    require(old['case_logical_content_sha256'] == 'c3670992a0ad33698af6c7620429bed48e0ecaa98adb3c87e648eefb1ba70363', 'PL2B identity mismatch')
    return m20, m21, mp


def candidate_and_reaction_identity() -> tuple[list[dict], list[dict], dict, list[str]]:
    candidates = pl2a.read_tsv(safe_file(str(pl2a.PINNED_SOURCES['candidate_table'][0])))
    candidate_by_id = {pl2a.candidate_identity(row): row for row in candidates}
    require(len(candidates) == len(candidate_by_id) == 640, 'candidate identity mismatch')
    aliases = pl2b.read_tsv(safe_file(str(PL2A / 'BRIDGEPL2A_TRUTH_PAIR_REGISTRY.tsv')))
    distinct, alias_map = pl2b.collapse_truth_pairs(aliases)
    require(len(distinct) == 885 and len(alias_map) == 1200, 'truth pair count mismatch')
    old_registry = pl2b.read_tsv(safe_file('outputs/dmi_bridge_pl2b_sign_only_utility_v1/BRIDGEPL2B_DISTINCT_TRUTH_PAIR_REGISTRY.tsv'))
    require({row['truth_pair_id'] for row in distinct} == {row['truth_pair_id'] for row in old_registry}, 'canonical truth IDs changed')
    for row in distinct:
        for prefix, tumor in (('ct2a', 'CT2A'), ('gl261', 'GL261')):
            cid = str(row[f'{prefix}_candidate_id'])
            require(cid in candidate_by_id, 'frozen candidate missing')
            candidate = candidate_by_id[cid]
            require(candidate['tumor'] == tumor and candidate['algorithm'] == row['algorithm'] and candidate['rna_context_key'] == row['rna_context_key'] and int(candidate['sample_index']) == row[f'{prefix}_sample_index'], 'frozen candidate context mismatch')
    cache_path = safe_file(str(pl2a.PINNED_SOURCES['flux_cache'][0]))
    checked_hash(str(pl2a.PINNED_SOURCES['flux_cache'][0]), pl2a.PINNED_SOURCES['flux_cache'][1])
    with np.load(cache_path, allow_pickle=False) as cache:
        rxns = [str(x) for x in cache['rxns'].tolist()]
    require(len(rxns) == 4181 and len(set(rxns)) == 4181 and rxns.index('HEX1') == 2793 and rxns.index('LDH_L') == 2909, 'cache reaction inventory mismatch')
    targets = a22.common_target_reactions(rxns)
    a1 = pl2b.read_tsv(safe_file(str(PL2A / 'BRIDGEPL2A_REACTION_REGISTRY.tsv')))
    reactions = [r for r in a1 if r['reaction_id'] != 'LDH_L']
    require(len(reactions) == 4179 and {r['reaction_id'] for r in reactions} == set(targets), 'A1/A2 target set mismatch')
    with lzma.open(safe_file(str(A21 / 'BRIDGEA21_REACTION_SUMMARY.tsv.xz')), 'rt') as handle:
        next(handle)
        a21_ids = {line.split('\t', 1)[0] for line in handle}
    require(a21_ids == set(targets), 'A2.1 reaction identity mismatch')
    return distinct, alias_map, candidate_by_id, reactions


def support_specs(distinct: list[dict], candidate_by_id: dict) -> tuple[list[dict], list[dict], dict]:
    rows = pl2a.read_tsv(safe_file(str(A20 / 'BRIDGEA20_DUAL_ANCHOR_WEIGHTS.tsv.xz')))
    require(len(rows) == 3200 and all(r['arm'] == 'A2-L' for r in rows), 'A2-L weights mismatch')
    groups = defaultdict(dict)
    for r in rows:
        key = (r['mouse_id'], r['tumor'], r['algorithm'], r['ensemble_hash'], r['projection_hash'], r['rna_context_key'])
        sample = int(r['sample_index'])
        require(sample not in groups[key], 'duplicate A2-L sample')
        groups[key][sample] = float(r['weight'])
    require(len(groups) == 160 and all(set(g) == set(range(20)) for g in groups.values()), 'A2-L stratum mismatch')
    specs, audit = [], []
    for pair in distinct:
        condition = {}
        for prefix, tumor in (('ct2a', 'CT2A'), ('gl261', 'GL261')):
            cid = str(pair[f'{prefix}_candidate_id'])
            candidate = candidate_by_id[cid]
            key = (str(pair[f'{prefix}_mouse']), tumor, pair['algorithm'], candidate['ensemble_hash'], candidate['projection_hash'], candidate['rna_context_key'])
            pool = groups.get(key)
            require(pool is not None and len(pool) == 20, 'A2-L selected stratum missing')
            sample = int(candidate['sample_index'])
            ordered = sorted(pool)
            weights = np.asarray([pool[i] for i in ordered], dtype=float)
            location = ordered.index(sample)
            try:
                held = a22.remove_holdout_and_normalize(weights, location)
                status = 'EVALUABLE'
            except ValueError as exc:
                if 'zero post-holdout mass' not in str(exc):
                    raise
                held, status = None, 'ZERO_POST_HOLDOUT_MASS'
            condition[prefix] = None if held is None else {'held_samples': np.asarray(ordered[:location] + ordered[location + 1:], dtype=int), 'held_weights': held.weights, 'truth_sample': sample, 'cache_key': (candidate['algorithm'], tumor, candidate['ensemble_hash'])}
            audit.append({'truth_pair_id': pair['truth_pair_id'], 'evaluation_id': pair['evaluation_id'], 'primary_matched_pair': bool(pair['pair_evaluable']), 'condition': prefix.upper(), 'mouse_id': key[0], 'candidate_id': cid, 'sample_index': sample, 'pre_candidates': 20, 'post_candidates': 19, 'positive_weight_candidates_pre': 0 if held is None else held.positive_weight_candidates_pre, 'positive_weight_candidates_post': 0 if held is None else held.positive_weight_candidates_post, 'post_holdout_mass': 0.0 if held is None else held.post_holdout_mass, 'post_holdout_ess': float('nan') if held is None else held.post_holdout_ess, 'status': status})
        if pair['pair_evaluable']:
            require(condition['ct2a'] is not None and condition['gl261'] is not None, 'strict 836-pair support gate failed')
            specs.append({'pair': pair, **condition})
    specs.sort(key=lambda s: s['pair']['truth_pair_id'])
    audit.sort(key=lambda r: (r['truth_pair_id'], r['condition']))
    require(len(specs) == 836 and len(audit) == 1770, 'support accounting mismatch')
    masses = [float(r['post_holdout_mass']) for r in audit if r['primary_matched_pair']]
    esses = [float(r['post_holdout_ess']) for r in audit if r['primary_matched_pair']]
    qc = {'matched_evaluable': len(specs), 'nonprimary_a2l_evaluable': sum(all(r['status'] == 'EVALUABLE' for r in audit if r['truth_pair_id'] == pair['truth_pair_id']) for pair in distinct if not pair['pair_evaluable']), 'primary_post_mass_min': min(masses), 'primary_post_mass_max': max(masses), 'primary_post_ess_min': min(esses), 'primary_post_ess_max': max(esses), 'post_pool_size': 19}
    return specs, audit, qc


def part_name(index: int) -> str:
    return f'BRIDGEA22_CASE_OUTCOMES.part-{index:03d}.tsv.xz'


def scalar_parity(scalar_result: dict, block_result: dict, position: int) -> tuple[float, float, float, int]:
    """Permit breakpoint existence noise only when all endpoint gains are tied."""
    adjusted = dict(scalar_result)
    adjusted_rows = []
    degenerate = 0
    for row in scalar_result['lambda_results']:
        value = dict(row)
        lam = float(value['lambda_total'])
        block = block_result['lambda_results'][lam]
        left = float(value['break_even_reliability'])
        right = float(block['break_even_reliability'][position])
        if math.isfinite(left) != math.isfinite(right):
            gains = (float(value['correct_abs_error_gain']), float(value['wrong_abs_error_gain']),
                     float(block['correct_absolute_error_gain'][position]), float(block['wrong_absolute_error_gain'][position]))
            require(max(abs(gain) for gain in gains) <= scalar.GAIN_TOL,
                    'meaningful scalar/batch break-even existence mismatch')
            value['break_even_reliability'] = right
            degenerate += 1
        adjusted_rows.append(value)
    adjusted['lambda_results'] = adjusted_rows
    maximum, break_even, scaled = pl2b._scalar_parity(adjusted, block_result, position)
    return maximum, break_even, scaled, degenerate


def write_part(path: Path, specs: list[dict], reactions: list[dict], cols: np.ndarray, mats: np.ndarray) -> dict:
    tmp = path.with_name(path.name + '.tmp')
    tmp.unlink(missing_ok=True)
    count, zero_max, parity_max, be_scale_max, degenerate_break_even = 0, 0.0, 0.0, 0.0, 0
    with lzma.open(tmp, 'wt', newline='', preset=0) as handle:
        writer = csv.DictWriter(handle, fieldnames=pl2b.CASE_FIELDS, delimiter='\t', lineterminator='\n', extrasaction='raise')
        writer.writeheader()
        for spec in specs:
            left, right = spec['ct2a'], spec['gl261']
            for start in range(0, len(reactions), BLOCK_SIZE):
                stop = min(start + BLOCK_SIZE, len(reactions))
                block_cols = cols[start:stop]
                result = batch.evaluate_truth_block(ct2a_values=mats[left['cache_index']][left['held_samples']][:, block_cols], gl261_values=mats[right['cache_index']][right['held_samples']][:, block_cols], ct2a_weights=left['held_weights'], gl261_weights=right['held_weights'], truth_ct2a=mats[left['cache_index'], left['truth_sample'], block_cols], truth_gl261=mats[right['cache_index'], right['truth_sample'], block_cols])
                non_tie = np.asarray(result['truth_direction']) != 0
                zero = result['lambda_results'][0.0]
                for field, baseline in (('correct_magnitude_estimate', 'baseline_magnitude_estimate'), ('wrong_magnitude_estimate', 'baseline_magnitude_estimate'), ('correct_posterior_ess', 'baseline_joint_ess'), ('wrong_posterior_ess', 'baseline_joint_ess')):
                    if np.any(non_tie):
                        zero_max = max(zero_max, float(np.max(np.abs(zero[field][non_tie] - result[baseline][non_tie]))))
                if start in (0, (len(reactions) // 2 // BLOCK_SIZE) * BLOCK_SIZE, ((len(reactions) - 1) // BLOCK_SIZE) * BLOCK_SIZE):
                    for pos in sorted({0, len(block_cols) - 1}):
                        delta, joint = scalar.joint_delta_distribution(mats[left['cache_index']][left['held_samples'], block_cols[pos]], mats[right['cache_index']][right['held_samples'], block_cols[pos]], left['held_weights'], right['held_weights'])
                        check = a22.evaluate_truth_case(delta_b=delta, baseline_weights=joint, truth_delta_b=float(result['truth_delta_b'][pos]))
                        parity, _, scaled, degenerate = scalar_parity(check, result, pos)
                        parity_max, be_scale_max = max(parity_max, parity), max(be_scale_max, scaled)
                        degenerate_break_even += degenerate
                for pos, reaction in enumerate(reactions[start:stop]):
                    row = pl2b._case_row(spec['pair'], reaction, result, pos)
                    writer.writerow({field: pl2b.format_value(row[field]) for field in pl2b.CASE_FIELDS})
                    count += 1
    os.replace(tmp, path)
    require(path.stat().st_size < MAX_PART_BYTES and zero_max <= 1e-12 and parity_max <= 5e-12 and be_scale_max <= 5e-12, 'part QC failed')
    return {'filename': path.name, 'rows': count, 'bytes': path.stat().st_size, 'sha256': pl2b.sha256_file(path), 'lambda_zero_maximum_deviation': zero_max, 'scalar_batch_maximum_deviation': parity_max, 'scalar_batch_break_even_gain_scale_maximum_deviation': be_scale_max, 'degenerate_break_even_sentinel_mismatches': degenerate_break_even}


def checked_part(path: Path, record: dict, expected_rows: int) -> None:
    require(path.is_file() and not path.is_symlink() and path.stat().st_size == record['bytes'] and pl2b.sha256_file(path) == record['sha256'] and record['rows'] == expected_rows, f'checkpoint part mismatch: {path.name}')


def compute_parts(out: Path, specs: list[dict], reactions: list[dict], cols: np.ndarray, mats: np.ndarray, fingerprint: dict) -> list[dict]:
    work = out / '_work'
    work.mkdir(parents=True, exist_ok=True)
    checkpoint = work / 'WORK_MANIFEST.json'
    if checkpoint.exists():
        state = json.loads(checkpoint.read_text())
        require(state['fingerprint'] == fingerprint, 'checkpoint fingerprint mismatch')
    else:
        require(not any(work.glob('BRIDGEA22_CASE_OUTCOMES.part-*')), 'parts without checkpoint')
        state = {'schema': 'bridge.a22.work.v1', 'fingerprint': fingerprint, 'completed_parts': {}}
        pl2b.atomic_json(checkpoint, state)
    records = []
    for index in range(PART_COUNT):
        start, stop = len(specs) * index // PART_COUNT, len(specs) * (index + 1) // PART_COUNT
        name = part_name(index)
        expected_rows = (stop - start) * len(reactions)
        if name in state['completed_parts']:
            record = state['completed_parts'][name]
            checked_part(work / name, record, expected_rows)
        else:
            require(not (work / name).exists(), 'unregistered checkpoint part')
            record = write_part(work / name, specs[start:stop], reactions, cols, mats)
            require(record['rows'] == expected_rows, 'part row count mismatch')
            state['completed_parts'][name] = record
            pl2b.atomic_json(checkpoint, state)
        records.append(record)
        print(f'A2.2 part {index + 1}/{PART_COUNT} verified ({record["rows"]} rows)', flush=True)
    return records


def inspect_parts(directory: Path, records: list[dict]) -> tuple[dict, dict]:
    total = a22.EXPECTED_CASE_ROWS
    gain_path = directory / '_work' / 'summary_gains.dat'
    gains = np.memmap(gain_path, mode='w+', dtype='f8', shape=(3, total, 2))
    contexts = np.empty(total, dtype=np.uint8)
    non_tie = np.zeros(total, dtype=bool)
    context_to_code = {}
    logical = hashlib.sha256()
    previous = None
    offset = 0
    gain_fields = [(f'{prefix}_correct_absolute_error_gain', f'{prefix}_wrong_absolute_error_gain') for prefix in pl2b.LAMBDA_PREFIX.values()]
    for index, rec in enumerate(records):
        path = directory / '_work' / rec['filename']
        checked_part(path, rec, rec['rows'])
        with lzma.open(path, 'rb') as handle:
            header_line = handle.readline()
            header = header_line.decode().rstrip('\n').split('\t')
            require(header == pl2b.CASE_FIELDS, 'case schema mismatch')
            if index == 0:
                logical.update(header_line)
            pos = {field: header.index(field) for field in ('truth_pair_id', 'reaction_id', 'algorithm', 'rna_context_key', 'truth_status')}
            gain_pos = [(header.index(a), header.index(b)) for a, b in gain_fields]
            count = 0
            for line in handle:
                fields = line.decode().rstrip('\n').split('\t')
                require(len(fields) == len(header), 'case row width mismatch')
                key = (fields[pos['truth_pair_id']], fields[pos['reaction_id']])
                require(previous is None or previous < key, 'case order mismatch')
                previous = key
                label = fields[pos['algorithm']] + '\x1f' + fields[pos['rna_context_key']]
                if label not in context_to_code:
                    context_to_code[label] = len(context_to_code)
                row = offset + count
                contexts[row] = context_to_code[label]
                non_tie[row] = fields[pos['truth_status']] == 'NON_TIE'
                for j, (ca, wr) in enumerate(gain_pos):
                    gains[j, row, 0], gains[j, row, 1] = float(fields[ca]), float(fields[wr])
                logical.update(line)
                count += 1
        require(count == rec['rows'], 'decompressed part row count mismatch')
        offset += count
    require(offset == total, 'total case row count mismatch')
    gains.flush()
    parts = {'schema': 'bridge.a22.case_outcomes.partitioned_tsv.v1', 'logical_artifact': 'BRIDGEA22_CASE_OUTCOMES', 'format': 'tsv.xz.parts', 'partition_scheme': '128_contiguous_canonical_truth_pair_ranges', 'sort_key': ['truth_pair_id', 'reaction_id'], 'number_of_parts': PART_COUNT, 'column_names': pl2b.CASE_FIELDS, 'parts': records, 'total_data_rows': total, 'logical_content_sha256': logical.hexdigest(), 'maximum_part_bytes': max(r['bytes'] for r in records), 'part_size_ceiling_bytes': MAX_PART_BYTES}
    summary = {'gains': gains, 'context_codes': contexts, 'context_to_code': context_to_code, 'non_tie': non_tie}
    return parts, summary


def validate_final(out: Path, fingerprint: dict) -> dict:
    manifest = out / 'BRIDGEA22_MANIFEST.json'
    require(manifest.is_file(), 'A2.2 output exists without final manifest')
    value = json.loads(manifest.read_text())
    require(value['status'] == 'BRIDGE_A22_SIGN_ONLY_UTILITY_FROZEN' and value['fingerprint'] == fingerprint, 'final identity mismatch')
    for name, digest in value['artifact_sha256'].items():
        require(pl2b.sha256_file(safe_file(str(out.relative_to(ROOT) / name))) == digest, f'final artifact mismatch: {name}')
    matched = pl2b.read_tsv(safe_file(str(out.relative_to(ROOT) / 'BRIDGEA22_MATCHED_TRUTH_PAIR_REGISTRY.tsv')))
    require(len(matched) == 836 and all(row['pair_evaluable'] == 'true' for row in matched), 'matched registry population mismatch')
    require(len(pl2b.read_tsv(safe_file(str(out.relative_to(ROOT) / 'BRIDGEA22_ALIAS_MAP.tsv')))) == 1200, 'alias registry mismatch')
    return value


def produce(out: Path = OUT) -> dict:
    require(out == OUT, 'A2.2 output root is fixed')
    m20, m21, mp = admission()
    codes = {str(path.relative_to(ROOT)): pl2b.sha256_file(path) for path in (Path(__file__), Path(a22.__file__), Path(batch.__file__), Path(scalar.__file__), Path(pl2b.__file__), Path(storage.__file__))}
    fingerprint = {'schema': 'bridge.a22.fingerprint.v1', 'predecessor_manifest_sha256': {k: EXPECTED[k] for k in EXPECTED if k.endswith('MANIFEST.json')}, 'a21_feature_logical_content_sha256': m21['feature_logical_content_sha256'], 'a20_weights_sha256': m20['artifact_sha256']['BRIDGEA20_DUAL_ANCHOR_WEIGHTS.tsv.xz'], 'implementation_sha256': codes, 'lambda_grid': list(scalar.LAMBDA_GRID), 'reliability_grid': list(scalar.RELIABILITY_GRID), 'parts': PART_COUNT, 'block_size': BLOCK_SIZE}
    if out.exists() and (out / 'BRIDGEA22_MANIFEST.json').exists():
        manifest=validate_final(out, fingerprint)
        finalize_gain_storage(out/'_work')
        return {'status': 'NO_OP_EXISTING_IDENTICAL_A22', 'manifest': manifest}
    if out.exists():
        require(set(p.name for p in out.iterdir()) <= {'_work'}, 'conflicting final artifacts')
    distinct, aliases, candidate_by_id, reactions = candidate_and_reaction_identity()
    specs, support, support_qc = support_specs(distinct, candidate_by_id)
    out.mkdir(parents=True, exist_ok=True)
    work = out / '_work'
    work.mkdir(exist_ok=True)
    support_path = work / 'BRIDGEA22_TRUTH_SUPPORT_AUDIT.tsv'
    if support_path.exists():
        require(pl2b.read_tsv(support_path) == [{k: pl2b.format_value(v) for k, v in row.items()} for row in support], 'support checkpoint mismatch')
    else:
        pl2b.write_tsv(support_path, list(support[0]), support)
    # Reaction-B values are first opened only after all identity and support gates.
    _, _, _, _, mats, cache_info = pl2a.load_and_validate_inputs()
    cache_index = {(a, t, e): i for i, (a, t, e) in enumerate(zip(cache_info['alg'], cache_info['tumor'], cache_info['eh']))}
    for spec in specs:
        for prefix in ('ct2a', 'gl261'):
            spec[prefix]['cache_index'] = cache_index[spec[prefix]['cache_key']]
    colmap = {r: i for i, r in enumerate(cache_info['rxns'])}
    cols = np.asarray([colmap[r['reaction_id']] for r in reactions], dtype=int)
    records = compute_parts(out, specs, reactions, cols, mats, fingerprint)
    del mats
    parts, summary_data = inspect_parts(out, records)
    summaries = pl2b.build_context_summary(summary_data)
    gains = summary_data['gains']
    gains.flush()
    gains._mmap.close()
    del summary_data
    staging = work / 'publish'
    require(not staging.exists(), 'conflicting publication staging directory')
    staging.mkdir()
    shutil.copy2(support_path, staging / support_path.name)
    pl2b.write_tsv(staging / 'BRIDGEA22_MATCHED_TRUTH_PAIR_REGISTRY.tsv', pl2b._registry_fields(), [pair for pair in distinct if pair['pair_evaluable']])
    pl2b.write_tsv(staging / 'BRIDGEA22_ALIAS_MAP.tsv', pl2b._alias_fields(pl2b.read_tsv(safe_file(str(PL2A / 'BRIDGEPL2A_TRUTH_PAIR_REGISTRY.tsv')))), aliases)
    pl2b.write_tsv(staging / 'BRIDGEA22_REACTION_REGISTRY.tsv', list(reactions[0]), reactions)
    pl2b.write_tsv(staging / 'BRIDGEA22_CONTEXT_SUMMARY.tsv', pl2b._summary_fields(), summaries)
    pl2b.atomic_json(staging / 'BRIDGEA22_CASE_OUTCOMES.parts.json', parts)
    source_audit = {'schema': 'bridge.a22.source_audit.v1', 'predecessor_manifest_sha256': fingerprint['predecessor_manifest_sha256'], 'a21_feature_logical_content_sha256': m21['feature_logical_content_sha256'], 'pl2a_artifact_sha256': mp['artifact_sha256'], 'a20_artifact_sha256': m20['artifact_sha256'], 'a21_artifact_sha256': m21['artifact_sha256'], 'implementation_sha256': codes, 'a2_g_used': False, 'reaction_b_access_after_support_gate': True}
    pl2b.atomic_json(staging / 'BRIDGEA22_SOURCE_AUDIT.json', source_audit)
    contract = a22.analysis_contract()
    contract.update({'partition_count': PART_COUNT, 'reaction_block_size': BLOCK_SIZE, 'primary_pair_weighting': 'DISTINCT_TRUTH_PAIR_ONCE', 'context_summary_fraction_denominator': 'NON_TIE_CASES'})
    pl2b.atomic_json(staging / 'BRIDGEA22_ANALYSIS_CONTRACT.json', contract)
    overall = next(r for r in summaries if r['summary_scope'] == 'OVERALL' and r['lambda_total'] == scalar.PRIMARY_LAMBDA and r['reliability_q'] == 0.5)
    qc = {'schema': 'bridge.a22.qc.v1', 'status': 'PASS', 'counts': {'truth_selection_aliases': 480, 'truth_pair_aliases': 1200, 'all_distinct_truth_pairs': 885, 'matched_primary_truth_pairs': 836, 'matched_registry_rows': 836, 'nonprimary_truth_pairs': 49, 'reactions': 4179, 'case_rows': parts['total_data_rows'], 'truth_tie_cases': overall['truth_tie_count'], 'non_tie_cases': overall['non_tie_count'], 'context_summary_rows': len(summaries)}, 'support': support_qc, 'lambda_zero_maximum_deviation': max(r['lambda_zero_maximum_deviation'] for r in records), 'scalar_batch_maximum_deviation': max(r['scalar_batch_maximum_deviation'] for r in records), 'scalar_batch_break_even_gain_scale_maximum_deviation': max(r['scalar_batch_break_even_gain_scale_maximum_deviation'] for r in records), 'degenerate_break_even_sentinel_mismatches': sum(r.get('degenerate_break_even_sentinel_mismatches', 0) for r in records), 'case_logical_content_sha256': parts['logical_content_sha256'], 'a2_g_used': False}
    pl2b.atomic_json(staging / 'BRIDGEA22_QC.json', qc)
    for rec in records:
        shutil.copy2(work / rec['filename'], staging / rec['filename'])
    names = sorted(p.name for p in staging.iterdir())
    hashes = {name: pl2b.sha256_file(staging / name) for name in names}
    manifest = {'schema': 'bridge.a22.sign_only_utility.manifest.v1', 'status': 'BRIDGE_A22_SIGN_ONLY_UTILITY_FROZEN', 'fingerprint': fingerprint, 'counts': qc['counts'], 'artifact_sha256': hashes, 'case_logical_content_sha256': parts['logical_content_sha256'], 'case_parts_manifest_sha256': hashes['BRIDGEA22_CASE_OUTCOMES.parts.json'], 'a2_g_used': False, 'primary_benchmark_unit': 'distinct_truth_pair_once', 'primary_lambda': scalar.PRIMARY_LAMBDA, 'a1_vs_a2_synthesis_computed': False}
    pl2b.atomic_json(staging / 'BRIDGEA22_MANIFEST.json', manifest)
    for p in staging.iterdir():
        os.replace(p, out / p.name)
    staging.rmdir()
    validate_final(out, fingerprint)
    finalize_gain_storage(work)
    return {'status': manifest['status'], 'manifest': manifest}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.parse_args()
    result = produce()
    print(json.dumps({'status': result['status'], 'counts': result['manifest']['counts'], 'case_logical_content_sha256': result['manifest']['case_logical_content_sha256']}, sort_keys=True))


if __name__ == '__main__':
    main()
