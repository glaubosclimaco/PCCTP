"""Bind existing TRAIN starts to completed audits; never run an optimizer."""
import hashlib
import json
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from experiments.lkh.raw_route_audit import load, check


def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def load_starts():
    source_audit = HERE/'ruin_total60_absolute_audit.json'
    audit = json.loads(source_audit.read_text())
    assert audit['status'] == 'PASS' and audit['batch_complete'] and audit['runs'] == 54
    rows = {(r['instance'], r['seed'], r['mode']): r for r in audit['rows']}
    assert len(rows) == 54
    cases_file = ROOT/'experiments/original_mdm/hard_training_baseline/summary.json'
    roles_file = ROOT/'experiments/ml_lns/decoder_initialized_states/protocol.json'
    profile_file = HERE/'hard_train_instance_profile.json'
    cases = json.loads(cases_file.read_text())['rows']
    roles = json.loads(roles_file.read_text())['roles']
    profiles = {r['instance']: r for r in json.loads(profile_file.read_text())['rows']}
    assert len(cases) == len(profiles) == 9
    assert {c['instance'] for c in cases} == set(profiles)
    starts = []
    for case in cases:
        name, geometry = case['instance'], case['geometry']
        assert roles[geometry] == 'train' and profiles[name]['geometry'] == geometry
        data = load(name)
        assert data['input_hash'] == profiles[name]['input_sha256']
        for seed in (42, 43, 44):
            source = HERE/'ruin_total60_absolute'/f'{name}.seed{seed}.hybrid.json'
            hashed = digest(source)
            assert rows[(name, seed, 'hybrid')]['result_sha256'] == hashed
            record = json.loads(source.read_text())
            assert (record['instance'], record['seed'], record['mode']) == (name, seed, 'hybrid')
            check(data, record['pre_route'], record['pre_cost'])
            starts.append(dict(instance=name, geometry=geometry, source_seed=seed,
                               route=record['pre_route'], cost=record['pre_cost'],
                               source_file=str(source.relative_to(ROOT)), source_sha256=hashed,
                               input_sha256=data['input_hash']))
    assert len(starts) == 27
    random.Random(20261020).shuffle(starts)
    files = [Path(__file__), HERE/'raw_route_audit.py', source_audit, cases_file, roles_file, profile_file]
    return dict(scope='Prospective conditional TRAIN starts only; generation excluded; no optimizer.',
                starts=starts, scheduled_pairs=27, scheduled_phases=54, budget_seconds=30.,
                repair_seed_rule='source_seed+49',
                provenance_sha256={str(p.relative_to(ROOT)): digest(p) for p in files})


if __name__ == '__main__':
    manifest = load_starts()
    (HERE/'cached30_starts.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(json.dumps(dict(status='PASS', starts=len(manifest['starts']),
                          geometries=len({s['geometry'] for s in manifest['starts']}),
                          native_optimizer_calls=0)))
