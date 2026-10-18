"""Audit completed scheduled pilot files; no timing superiority inference."""
import hashlib
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from experiments.ml_lns.pilot import read_instance,validate

def main():
    out=HERE/'scheduled_development';rows=[]
    for p in sorted(out.glob('g_*.json')):
        if p.name.endswith('.adapter_check.json'):continue
        r=json.loads(p.read_text())
        if 'operator_policy' not in r:continue
        fingerprint=p.stem.split('.')[-1]
        protocol=json.loads((out/f'protocol.{fingerprint}.json').read_text())
        case=next(c for c in protocol['cases'] if c['instance']==r['instance'] and c['seed']==r['seed'])
        assert case['role']=='development'
        baseline=ROOT/'experiments/ml_lns/decoder_initialized_states'/f"{r['instance']}.seed{r['seed']}.json"
        assert hashlib.sha256(baseline.read_bytes()).hexdigest()==case['source_sha256']
        assert r['method']=='scheduled_MDM_'+r['operator_policy']
        d=read_instance(ROOT/'data'/r['instance']);mapping={int(v):i for i,v in enumerate(d['visits'])}
        for state in [r,*r['snapshots']]:assert abs(validate(d,[mapping[v] for v in state['route']])-state['cost'])<1e-6
        costs=[s['cost'] for s in r['snapshots']];assert all(a>=b for a,b in zip(costs,costs[1:]))
        if r['completed']:assert r['iterations']==70 and not r['timed_out']
        rows.append(dict(file=p.name,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
                         policy=r['operator_policy'],instance=r['instance'],cost=r['cost'],
                         seconds=r['process_wall_seconds'],completed=r['completed']))
    result=dict(status='PASS',files=len(rows),rows=rows,
                scope='Completed files at audit time; fixed iterations, sparse checkpoints and additional timing wrapper; no speed claim')
    (HERE/'scheduled_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='rows'}))

if __name__=='__main__':main()
