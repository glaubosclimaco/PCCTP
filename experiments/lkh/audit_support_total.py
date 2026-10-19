"""Audit all saved routes, preserving infeasible-before-deadline observations."""
import json,sys,hashlib,statistics
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from experiments.ml_lns.pilot import read_instance,validate
HERE=Path(__file__).resolve().parent;OUT=HERE/'support_total'

def audit():
    assert (OUT/'completion.json').exists(),'Batch still running'
    protocol=json.loads((OUT/'protocol.json').read_text());modes=protocol['modes'];rows=[];hashes={};failures={m:0 for m in modes};times={m:[] for m in modes}
    for c in protocol['cases']:
        name=c['instance'];d=read_instance(ROOT/'data'/name);mapping={int(v):i for i,v in enumerate(d['visits'])};pairs=[]
        for seed in protocol['seeds']:
            pair={}
            for m in modes:
                p=OUT/f'{name}.seed{seed}.{m}.json';np=OUT/f'{name}.seed{seed}.{m}.native.json'
                r=json.loads(p.read_text());native=json.loads(np.read_text())
                for f in (p,np):hashes[f.name]=hashlib.sha256(f.read_bytes()).hexdigest()
                assert (r['instance'],r['seed'],r['mode'])==(name,seed,m)
                assert native['input_sha256']==hashlib.sha256((ROOT/'data'/name).read_bytes()).hexdigest()
                budget=15 if m=='mdm' else 10
                assert native['budget_seconds']==budget and native['elite_count_bounded']
                previous=float('inf');last=-1
                for e in native['anytime']:
                    assert last<=e['seconds']<=budget and e['cost']<previous
                    assert validate(d,[mapping[v] for v in e['route']])==e['cost']
                    previous=e['cost'];last=e['seconds']
                assert bool(native['anytime'])==r['feasible']
                if r['feasible']:
                    assert validate(d,[mapping[v] for v in r['route']])==r['cost']<=native['cost']
                    if r['cost']<native['cost']:assert r['wall_seconds']<=15 and r['post'] and r['post']['cost']==r['cost']
                    if r['post']:
                        post=r['post'];assert post['initial_cost']==native['cost']
                        assert validate(d,[mapping[v] for v in post['route']])==post['cost']
                else:
                    assert r['cost'] is None and r['route']==[];failures[m]+=1
                times[m].append(r['wall_seconds']);pair[m]=dict(feasible=r['feasible'],cost=r['cost'])
            pairs.append(dict(seed=seed,arms=pair))
        comparisons={}
        for control in ('mdm','fixed_lkh','support_2opt'):
            joint=[p for p in pairs if p['arms'][control]['feasible'] and p['arms']['support_lkh']['feasible']]
            comparisons[control]=dict(jointly_feasible_pairs=len(joint),gain_percent=statistics.mean(100*(p['arms'][control]['cost']-p['arms']['support_lkh']['cost'])/p['arms'][control]['cost'] for p in joint) if joint else None)
            if len(joint)==3:
                diff=statistics.mean(p['arms'][control]['cost']-p['arms']['support_lkh']['cost'] for p in joint)
                comparisons[control]['classification']='win' if diff>0 else 'loss' if diff<0 else 'tie'
            else:comparisons[control]['classification']='incomplete_feasibility'
        rows.append(dict(instance=name,pairs=pairs,comparisons=comparisons))
    summary={}
    for control in ('mdm','fixed_lkh','support_2opt'):
        complete=[r['comparisons'][control] for r in rows if r['comparisons'][control]['jointly_feasible_pairs']==3]
        summary[control]=dict(fully_feasible_geometries=len(complete),wins=sum(r['classification']=='win' for r in complete),losses=sum(r['classification']=='loss' for r in complete),mean_gain_percent=statistics.mean(r['gain_percent'] for r in complete) if complete else None)
    gate=not any(failures.values()) and summary['mdm']['mean_gain_percent']>=.5 and summary['mdm']['wins']>=6 and summary['mdm']['losses']<=1 and summary['fixed_lkh']['wins']>=3 and summary['support_2opt']['wins']>=3
    result=dict(status='PASS',scope='exploratory TRAIN; incomplete feasibility is retained, not silently excluded',feasibility_failures=failures,comparisons=summary,gate_passed=gate,rows=rows,result_hashes=hashes,timing={m:dict(mean=statistics.mean(t),maximum=max(t)) for m,t in times.items()})
    (HERE/'support_total_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('rows','result_hashes')}))
if __name__=='__main__':audit()
