"""Raw scalar audit; verify completeness before interpreting full-model bounds."""
import json, hashlib
from pathlib import Path
from raw_route_audit import load, check
from audit_edge_pool import edges
HERE = Path(__file__).resolve().parent

if __name__ == '__main__':
    folder = HERE/'full_graph60'; rows=[]; hashes={}; checked=0
    for p in sorted(folder.glob('g_*.json')):
        r=json.loads(p.read_text()); d=load(r['instance'])
        source=HERE/'edge_pool_study'/(r['instance']+'.json')
        old=json.loads(source.read_text())
        assert r['source_sha256']==hashlib.sha256(source.read_bytes()).hexdigest()
        assert r['initial_cost']==old['initial_cost'] and r['arm']=='full'
        visits=[v for v,t in enumerate(d['types']) if t!=2]
        full={(a,b) for i,a in enumerate(visits) for b in visits[i+1:]}
        assert set(map(tuple,r['pool_edges']))==full and len(r['pool_edges'])==len(full)
        for key in ('result','raw_result'):
            result=r[key];check(d,result['route'],result['cost']);checked+=1
            assert edges(result['route'])<=full and result['cost']<=r['initial_cost']
        late=r['native_return_seconds']>60 and r['raw_result']['cost']<r['initial_cost']
        assert r['late_improvement_discarded']==late
        if late:assert r['result']['cost']==r['initial_cost']
        else:assert r['result']==r['raw_result']
        rows.append(dict(instance=r['instance'],initial_cost=r['initial_cost'],cost=r['result']['cost'],
                         gain_percent=100*(r['initial_cost']-r['result']['cost'])/r['initial_cost'],
                         native_seconds=r['native_return_seconds'],complete_graph_verified=True))
        hashes[p.name]=hashlib.sha256(p.read_bytes()).hexdigest()
    complete=(folder/'completion.json').exists()
    if complete:assert len(rows)==9
    mean=sum(r['gain_percent'] for r in rows)/len(rows) if rows else None
    wins=sum(r['gain_percent']>0 for r in rows)
    report=dict(status='PASS' if rows else 'NO_RESULTS',cases=len(rows),routes_checked=checked,
                batch_complete=complete,mean_gain_percent=mean,wins=wins,
                gate_passed=complete and wins>=3 and mean>=.2,
                scope='raw scalar routes and complete graph; no independent solver-bound or elapsed-clock certification',
                rows=rows,result_hashes=hashes)
    (HERE/'full_graph60_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('rows','result_hashes')}))
