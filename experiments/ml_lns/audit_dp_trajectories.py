"""Independent final-route, monotonicity, provenance and footprint audit."""
import hashlib
import json
from pathlib import Path
from pilot import read_instance,validate

HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[1]

def main():
    checked=[]
    for path in sorted((HERE/'dp_trajectory_development').glob('g_*.json')):
        row=json.loads(path.read_text())
        protocol=json.loads((path.parent/f"protocol.{row['protocol']}.json").read_text())
        case=next(c for c in protocol['evaluation_cases'] if c['instance']==row['instance'] and c['seed']==row['seed'])
        source=HERE/'decoder_initialized_states'/f"{row['instance']}.seed{row['seed']}.json"
        assert hashlib.sha256(source.read_bytes()).hexdigest()==case['baseline_sha256']
        baseline=json.loads(source.read_text());assert baseline['role']=='development'
        d=read_instance(ROOT/'data'/row['instance']);mapping={int(v):i for i,v in enumerate(d['visits'])}
        assert abs(validate(d,[mapping[v] for v in row['route']])-row['cost'])<1e-6
        initial=next(s for s in baseline['snapshots'] if s['iteration']==1)
        assert row['initial_cost']==initial['cost']
        value=row['initial_cost'];elapsed=0.
        for r in row['trace']:
            assert r['cost']<=value and r['elapsed']>=elapsed
            assert set(r['mandatory_released'])<=d['mandatory']
            assert not set(r['optional_released'])&d['mandatory']
            assert r['components']<=2*(len(r['optional_released'])+len(r['mandatory_released']))
            value=r['cost'];elapsed=r['elapsed']
        assert row['cost']==value and row['queries']==len(row['trace'])
        assert row['total_seconds']>=elapsed
        assert row['size_limits']==sum(r['reason']=='component_limit' for r in row['trace'])
        checked.append(dict(file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    result=dict(status='PASS',trajectories=len(checked),files=checked,
                scope='Completed files only; not certification of timing fairness or novelty')
    (HERE/'dp_trajectory_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='files'}))

if __name__=='__main__':main()
