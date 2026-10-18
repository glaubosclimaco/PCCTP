"""Compare constructive lifts and recorded MIP queries in identical fixed neighborhoods."""
import hashlib
import json
from pathlib import Path
import time
from pilot import read_instance,validate
from selection_lift import best_lift

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]

def main():
    out=HERE/'lifting_development'; out.mkdir(exist_ok=True)
    states=[]
    for path in sorted((HERE/'timed_oracle_decoder_initialized_states').glob('*.json')):
        source=json.loads(path.read_text())
        if source['role'] not in ('development','train'): raise ValueError('Reserved labels')
        digest=hashlib.sha256(path.read_bytes()).hexdigest()
        target=out/(path.stem+'.'+digest[:10]+'.json')
        if target.exists():
            states.append(json.loads(target.read_text())); continue
        d=read_instance(ROOT/'data'/source['instance'])
        mapping={int(v):i for i,v in enumerate(d['visits'])}
        route=[mapping[v] for v in source['initial_route']]
        initial=validate(d,route); rows=[]
        for r in source['records']:
            u={mapping[v] for v in r['released']}
            began=time.perf_counter(); lifted=best_lift(d,route,u,limit=3,seed=source['seed'])
            seconds=time.perf_counter()-began
            assert lifted['cost']<=initial+1e-7
            if r['proven_optimal']:
                assert lifted['cost']>=r['cost']-1e-6,('Lift below certified neighborhood optimum',path.name,r,lifted)
            rows.append(dict(released=r['released'],mip_cost=r['cost'],mip_optimal=r['proven_optimal'],
                             mip_query_seconds=r.get('query_wall_seconds'),lift_cost=lifted['cost'],
                             lift_seconds=seconds,feasible_assignments=lifted['feasible_assignments'],
                             lift_route=[int(d['visits'][v]) for v in lifted['route']]))
        row=dict(instance=source['instance'],geometry=source['geometry'],phase=source['phase'],seed=source['seed'],
                 initial_cost=initial,source_sha256=digest,assignment_limit=3,records=rows,
                 mip_improving=sum(r['mip_cost']<initial-1e-7 for r in rows),
                 lift_improving=sum(r['lift_cost']<initial-1e-7 for r in rows),
                 lift_strictly_better=sum(r['lift_cost']<r['mip_cost']-1e-7 for r in rows),
                 lift_matches_improving_mip=sum(r['mip_cost']<initial-1e-7 and r['lift_cost']<=r['mip_cost']+1e-7 for r in rows))
        target.write_text(json.dumps(row,indent=2)+'\n'); states.append(row)
        print(json.dumps({k:v for k,v in row.items() if k!='records'}),flush=True)
    summary=dict(scope='Development fixed-neighborhood diagnostic, no novelty/end-to-end speed claim',
                 states=len(states),queries=sum(len(s['records']) for s in states),
                 **{k:sum(s[k] for s in states) for k in ('mip_improving','lift_improving','lift_strictly_better','lift_matches_improving_mip')})
    (HERE/'lifting_development_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary),flush=True)

if __name__=='__main__': main()
