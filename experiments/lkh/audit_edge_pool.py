"""Audit offline edge-pool outputs independently of the model reader."""
import hashlib
import json
from pathlib import Path
from raw_route_audit import load, check

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent


def edges(route):
    return {tuple(sorted((route[i - 1], v))) for i, v in enumerate(route)}


def audit(row):
    d = load(row['instance'])
    routes = row['input_routes']
    assert routes
    costs = []
    for route in routes:
        cost = sum(d['matrix'][route[i - 1] * d['n'] + v]
                   for i, v in enumerate(route))
        costs.append(check(d, route, cost))
    pool = set().union(*(edges(route) for route in routes))
    recorded = [tuple(edge) for edge in row['pool_edges']]
    assert len(recorded) == len(set(recorded))
    assert set(recorded) == pool
    assert row['initial_cost'] == min(costs)
    for key in ('result', 'raw_result'):
        r = row[key]
        check(d, r['route'], r['cost'])
        assert edges(r['route']) <= pool
        assert r['cost'] <= min(costs)
        assert r['proof_scope'] == 'restricted edge pool'
    discarded = row['native_return_seconds'] > 60 and row['raw_result']['cost'] < min(costs)
    assert row['late_improvement_discarded'] == discarded
    assert row['native_elapsed'] >= row['native_return_seconds'] >= row['initial_ready_native_seconds']
    if row['late_improvement_discarded']:
        assert discarded and row['result']['cost'] == min(costs)
    else:
        assert row['result'] == row['raw_result']
    assert row['feasible_before_deadline'] == (row['initial_ready_native_seconds'] <= 60)
    for relative, digest in row['input_hashes'].items():
        path = ROOT / relative
        assert path.is_relative_to(ROOT)
        assert hashlib.sha256(path.read_bytes()).hexdigest() == digest
    return len(routes) + 2


if __name__ == '__main__':
    folder = HERE / 'edge_pool_study'
    files = sorted(p for p in folder.glob('*.json') if p.name != 'completion.json')
    checked = sum(audit(json.loads(p.read_text())) for p in files)
    complete = (folder / 'completion.json').exists()
    if complete:
        assert len(files) == 9
    report = dict(status='PASS' if files else 'NO_RESULTS', cases=len(files),
                  routes_checked=checked, batch_complete=complete,
                  scope='raw scalar feasibility, costs, exact pool union, input hashes; no global optimality or deadline certification',
                  result_hashes={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
    (HERE / 'edge_pool_audit.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report))
