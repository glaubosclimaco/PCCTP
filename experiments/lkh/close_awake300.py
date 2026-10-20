"""Freeze a completed audited batch and descriptive geometry sensitivity.

No optimizer, observer refresh, reserved data, new gate or inferential test.
Run only after the existing queue has finished the complete-batch observers.
"""
import hashlib
import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from experiments.lkh.raw_route_audit import load, check


def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def close():
    batch = HERE/'ruin_total300_awake_train'
    assert (batch/'completion.json').exists(), 'Wait for all81 runs; no partial closing analysis'
    names = ['ruin_total300_awake_train_audit.json', 'ruin_total300_awake_train_summary.json',
             'ruin_total300_awake_environment_audit.json', 'ruin_total300_awake_train_anytime.json',
             'ruin_total300_awake_train_mining_exposure_bounds.json']
    paths = [HERE/name for name in names]
    audit, summary, env, curve, mining = [json.loads(p.read_text()) for p in paths]
    assert audit['status'] == 'PASS' and audit['batch_complete'] and audit['runs'] == 81
    assert summary['status'] == curve['status'] == mining['status'] == 'COMPLETE'
    assert summary['runs'] == curve['runs'] == mining['runs'] == 81
    assert summary['completed_geometries'] == 9 and not summary['feasibility_failures']
    assert env['batch_complete'] and env['eligible_complete_environment'] and env['validated_blocks'] == 27
    assert curve['audit_sha256'] == digest(paths[0])
    protocol = json.loads((batch/'protocol.json').read_text())
    cases = {c['instance'] for c in protocol['cases']}
    records = {(r['instance'], r['seed'], r['mode']): r for r in audit['rows']}
    assert len(records) == len(audit['rows']) == 81
    assert set(records) == {(n, s, m) for n in cases for s in (70, 71, 72)
                            for m in ('mdm', 'hybrid10', 'hybrid60')}
    assert {(r['instance'],r['seed']) for r in env['rows']} == {(n,s) for n in cases for s in (70,71,72)}
    for row in env['rows']:
        path = batch/f"environment.{row['instance']}.seed{row['seed']}.json"
        assert digest(path) == row['environment_sha256']
    result_hashes, input_hashes, geometries = {}, {}, []
    for name in sorted(cases):
        data = load(name); input_hashes[name] = data['input_hash']
        costs = {m: [] for m in ('mdm', 'hybrid10', 'hybrid60')}
        for seed in (70, 71, 72):
            for mode in costs:
                row = records[name, seed, mode]
                path = batch/f'{name}.seed{seed}.{mode}.json'
                assert row['feasible'] and digest(path) == row['result_sha256']
                record = json.loads(path.read_text())
                assert record['cost'] == row['cost']
                check(data, record['route'], record['cost'])
                result_hashes[path.name] = digest(path)
                costs[mode].append(Fraction(str(record['cost'])))
        means = {m: sum(values)/3 for m, values in costs.items()}
        gains = {m: 100*(means['mdm']-means[m])/means['mdm'] for m in ('hybrid10', 'hybrid60')}
        geometries.append(dict(instance=name, means={m: float(v) for m,v in means.items()},
            exact_gains={m: str(v) for m,v in gains.items()},
            paired_seed_counts={m: dict(better=sum(x<y for x,y in zip(costs[m],costs['mdm'])),
                equal=sum(x==y for x,y in zip(costs[m],costs['mdm'])),
                worse=sum(x>y for x,y in zip(costs[m],costs['mdm']))) for m in gains}))
    comparisons = {}
    for mode in ('hybrid10','hybrid60'):
        gains = [Fraction(g['exact_gains'][mode]) for g in geometries]
        mean = sum(gains)/9
        wins, losses = sum(g>Fraction(1,10**9) for g in gains),sum(g< -Fraction(1,10**9) for g in gains)
        gate = mean >= Fraction(1,5) and wins >= 6 and losses <= 1
        old = summary['comparisons'][mode]
        assert abs(float(mean)-old['mean_gain_percent']) < 1e-9
        assert (wins,losses,gate) == (old['wins'],old['losses'],old['gate_passed'])
        leave_one_out = {g['instance']: str((sum(gains)-gain)/8) for g,gain in zip(geometries,gains)}
        comparisons[mode] = dict(exact_mean_gain_percent=str(mean),mean_gain_percent=float(mean),
            wins=wins,losses=losses,ties=9-wins-losses,gate_passed=gate,
            descriptive_leave_one_geometry_out=leave_one_out,
            descriptive_range=[float(min(map(Fraction,leave_one_out.values()))),
                               float(max(map(Fraction,leave_one_out.values())))])
    report = dict(status='COMPLETE_DESCRIPTIVE_TRAIN', runs=81, geometries=geometries, comparisons=comparisons,
        scope='Selected nine TRAIN geometries; exploratory complete300s screen only. '
              'Leave-one-out is sensitivity, not exclusion, inference, a new gate or independent confirmation. '
              'Not superiority over completed MDM70.',
        input_sha256=input_hashes,result_sha256=result_hashes,
        upstream_sha256={str(p.relative_to(ROOT)):digest(p) for p in paths},
        source_sha256={str(Path(__file__).relative_to(ROOT)):digest(Path(__file__)),
                       'experiments/lkh/raw_route_audit.py':digest(HERE/'raw_route_audit.py')})
    snapshot = HERE/'awake300_completed_snapshots'/digest(paths[0])
    snapshot.mkdir(parents=True,exist_ok=True)
    for source in paths+[batch/'protocol.json',batch/'completion.json',batch/'environment_audit.json']:
        # Distinct parent prefix for the batch's own generic filenames.
        target = snapshot/(source.parent.name+'__'+source.name)
        content = source.read_bytes()
        if target.exists(): assert target.read_bytes() == content
        else:
            with target.open('xb') as handle: handle.write(content)
    report['snapshot_directory'] = str(snapshot.relative_to(ROOT))
    target = HERE/'awake300_closing_analysis.json'
    content = json.dumps(report,indent=2)+'\n'
    if target.exists(): assert target.read_text() == content
    else:
        with target.open('x') as handle: handle.write(content)
    print(json.dumps(dict(status=report['status'],comparisons=comparisons)))


if __name__ == '__main__': close()
