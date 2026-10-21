"""Supplementary seed dispersion; never changes the prespecified primary gate."""
import hashlib
import json
import math
import statistics
from pathlib import Path

HERE = Path(__file__).resolve().parent


def describe(values):
    assert len(values) == 10
    assert all(type(x) in (int, float) and math.isfinite(x) for x in values)
    return dict(n=10, mean=statistics.mean(values), median=statistics.median(values),
                minimum=min(values), maximum=max(values), sample_sd=statistics.stdev(values))


def report(audit, protocol):
    assert audit['status'] == 'PASS'
    assert protocol['seeds'] == list(range(81, 91))
    expected = {c['instance'] for c in protocol['cases']}
    assert len(expected) == len(protocol['cases']) == 9
    keys = [(r['instance'], r['seed']) for r in audit['rows']]
    assert len(keys) == len(set(keys))
    assert all(name in expected and seed in range(81, 91) for name, seed in keys)
    geometries = []
    for case in protocol['cases']:
        rows = sorted((r for r in audit['rows'] if r['instance'] == case['instance']),
                      key=lambda r: r['seed'])
        complete = (len(rows) == 10 and all(r['status'] == 'PAIR_COMPLETED' for r in rows))
        item = dict(instance=case['instance'], geometry=case['geometry'], observed_pairs=len(rows),
                    complete=complete)
        if complete:
            reference = [r['reference_cost'] for r in rows]
            assert all(x > 0 for x in reference)
            item['reference_cost'] = describe(reference)
            item['reference_defined_seconds'] = describe([r['duration'] for r in rows])
            item['candidates'] = {}
            for mode in ('hybrid10', 'hybrid60'):
                costs = [r['candidate_costs'][mode] for r in rows]
                assert all(x > 0 for x in costs)
                deltas = [a-b for a, b in zip(reference, costs)]
                item['candidates'][mode] = dict(cost=describe(costs), paired_cost_reduction=describe(deltas),
                    wins=sum(x > 0 for x in deltas), ties=sum(x == 0 for x in deltas),
                    losses=sum(x < 0 for x in deltas))
            item['seed_pairs'] = [dict(seed=r['seed'], reference_cost=r['reference_cost'],
                duration=r['duration'], candidate_costs=r['candidate_costs']) for r in rows]
        geometries.append(item)
    return dict(status='DESCRIPTIVE_ONLY', planned_geometries=9, planned_pairs=90,
        completed_geometries=sum(g['complete'] for g in geometries), geometries=geometries,
        scope='Seed variability conditional on each selected TRAIN geometry. Sample SD uses n-1; no normality assumption, CI, p-value, cross-geometry inference or gate change. Incomplete geometries retain counts without cost summaries; paired reductions are raw costs, not the primary geometry percentage effect.')


if __name__ == '__main__':
    a, p = HERE/'paired70_train_audit.json', HERE/'paired70_train/protocol.json'
    if not a.exists() or not p.exists():
        print('NO_RESULTS'); raise SystemExit(0)
    result = report(json.loads(a.read_text()), json.loads(p.read_text()))
    result['source_sha256'] = {str(x): hashlib.sha256(x.read_bytes()).hexdigest() for x in (a, p, Path(__file__))}
    (HERE/'paired70_dispersion.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'geometries'}))
