"""Independent raw scalar checks of expansion arms and frozen decision rule."""
import json, random, hashlib
from pathlib import Path
from raw_route_audit import load, check
from audit_edge_pool import edges
HERE = Path(__file__).resolve().parent

if __name__ == '__main__':
    folder = HERE / 'edge_expansion'; rows = []; hashes = {}; checked = 0
    for source in sorted((HERE / 'edge_pool_study').glob('g_*.json')):
        old = json.loads(source.read_text()); d = load(old['instance']); n = d['n']
        visits = [v for v in range(n) if d['types'][v] != 2]
        base = set().union(*(edges(r) for r in old['input_routes']))
        nearest = set()
        for a in visits:
            ranked = sorted((b for b in visits if b != a), key=lambda b: (d['matrix'][a*n+b], b))
            nearest.update(tuple(sorted((a,b))) for b in ranked[:5])
        count = len(nearest - base)
        absent = [(a,b) for i,a in enumerate(visits) for b in visits[i+1:] if (a,b) not in base]
        pools = dict(nearest=base | nearest, random=base | set(random.Random(42).sample(absent,count)))
        pair = {}
        for arm in ('nearest','random'):
            p = folder / (source.stem + '.' + arm + '.json'); r = json.loads(p.read_text())
            assert r['source_sha256'] == hashlib.sha256(source.read_bytes()).hexdigest()
            assert r['arm'] == arm and r['instance'] == old['instance']
            assert r['initial_cost'] == old['initial_cost']
            recorded = [tuple(e) for e in r['pool_edges']]
            assert len(recorded) == len(set(recorded)) and set(recorded) == pools[arm]
            assert r['added_edges'] == count
            for key in ('result','raw_result'):
                result = r[key]; check(d,result['route'],result['cost']); checked += 1
                assert edges(result['route']) <= pools[arm]
                assert result['cost'] <= r['initial_cost']
                assert result['proof_scope'] == 'restricted edge pool'
            late = r['native_return_seconds'] > 60 and r['raw_result']['cost'] < r['initial_cost']
            assert r['late_improvement_discarded'] == late
            assert r['native_return_seconds'] >= r['initial_ready_native_seconds']
            if late: assert r['result']['cost'] == r['initial_cost']
            else: assert r['result'] == r['raw_result']
            pair[arm] = r['result']['cost']; hashes[p.name] = hashlib.sha256(p.read_bytes()).hexdigest()
        rows.append(dict(instance=old['instance'],initial_cost=old['initial_cost'],**pair))
    gains = [100*(r['initial_cost']-r['nearest'])/r['initial_cost'] for r in rows]
    wins = sum(g>0 for g in gains); control_wins = sum(r['nearest']<r['random'] for r in rows)
    mean = sum(gains)/len(gains)
    report = dict(status='PASS',cases=len(rows),routes_checked=checked,
                  scope='raw costs, feasibility, reconstructed candidate sets and discard rule; no global bound or runtime certification',
                  mean_gain_percent=mean,wins_over_input=wins,wins_over_random=control_wins,
                  gate_passed=wins>=3 and mean>=.2 and control_wins>=5,rows=rows,result_hashes=hashes)
    (HERE/'edge_expansion_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('rows','result_hashes')}))
