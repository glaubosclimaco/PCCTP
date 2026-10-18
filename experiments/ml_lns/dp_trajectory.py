"""Development-only complete DP search; no transferred ML ranking or novelty claim.

All online construction, choice and repair time is charged. Each policy uses
the same DP cap and per-query limit. Mandatory anchors retain their visits.
"""
import argparse
import hashlib
import json
from pathlib import Path
import random
import time
from dp_repair import DPRepair,library
from pilot import candidates,read_instance,validate

HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[1]
METHODS=['optional_uniform','anchor_uniform','anchor_random','anchor_alns']

def run(d,initial,method,budget,seed):
    began=time.perf_counter(); model=DPRepair(d); route=initial[:]
    value=validate(d,route); first=value; trace=[]; rng=random.Random(seed)
    weights={k:1. for k in ('random','coverage','prize')}; pool_id=0
    mapping={int(v):i for i,v in enumerate(d['visits'])}
    while time.perf_counter()-began<budget:
        pool=candidates(d,route,12 if method=='optional_uniform' else 8,12,seed+10007*pool_id)
        pool_id+=1
        while pool and time.perf_counter()-began<budget:
            if method=='anchor_alns':
                kinds=sorted({k for k,u in pool})
                kind=rng.choices(kinds,weights=[weights[k] for k in kinds])[0]
                index=rng.choice([i for i,(k,u) in enumerate(pool) if k==kind])
            else:index=rng.randrange(len(pool))
            kind,u=pool.pop(index); anchors=[]
            if method!='optional_uniform':
                if method=='anchor_random':anchors=rng.sample(sorted(d['mandatory']),min(4,len(d['mandatory'])))
                else:anchors=sorted(d['mandatory'],key=lambda v:(min(d['costs'][v,w] for w in u),v))[:4]
            remaining=budget-(time.perf_counter()-began)
            if remaining<=0:break
            result=model.solve(route,u|set(anchors),min(.2,remaining),'fragment_clique',seed=seed)
            proposed=[mapping[v] for v in result['route']]
            assert abs(validate(d,proposed)-result['cost'])<1e-6
            improved=result['cost']<value-1e-8
            weights[kind]=max(.05,.9*weights[kind]+.1*(6. if improved else 0.))
            if improved:route=proposed;value=result['cost']
            trace.append(dict(elapsed=time.perf_counter()-began,cost=value,kind=kind,
                              optional_released=sorted(u),mandatory_released=anchors,
                              query_seconds=result['total_seconds'],components=result['contracted_nodes'],
                              reason=result['reason'],proven_optimal=result['proven_optimal']))
            if improved:break
        if not pool and not trace:break
    return dict(method=method,budget_seconds=budget,total_seconds=time.perf_counter()-began,
                initial_cost=first,cost=validate(d,route),queries=len(trace),trace=trace,
                size_limits=sum(r['reason']=='component_limit' for r in trace),
                route=[int(d['visits'][v]) for v in route])

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--budget',type=float,default=5.)
    p.add_argument('--limit',type=int,default=6)
    args=p.parse_args(); cases=[]; seen=set()
    for path in sorted((HERE/'decoder_initialized_states').glob('g_*.json')):
        b=json.loads(path.read_text())
        if b['role']!='development' or b['geometry'] in seen:continue
        s=next((s for s in b['snapshots'] if s['iteration']==1),None)
        if s is None:continue
        cases.append((path,b,s));seen.add(b['geometry'])
        if len(cases)>=args.limit:break
    protocol=dict(scope='Development control pilot; historical sparse timings, no final speed claim',
                  methods=METHODS,budget=args.budget,query_limit=.2,component_limit=18,
                  optional_k=12,anchored_optional_k=8,mandatory_k=4,
                  evaluation_cases=[dict(instance=b['instance'],geometry=b['geometry'],seed=b['seed'],
                                         baseline_sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p,b,s in cases],
                  code_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in
                               ['dp_trajectory.py','dp_repair.py','path_dp.cpp','pilot.py']})
    fingerprint=hashlib.sha256(json.dumps(protocol,sort_keys=True).encode()).hexdigest()[:12]
    out=HERE/'dp_trajectory_development';out.mkdir(exist_ok=True)
    (out/f'protocol.{fingerprint}.json').write_text(json.dumps(protocol,indent=2)+'\n')
    library() # One-off native compilation excluded; online model loading is inside run.
    for path,b,s in cases:
        d=read_instance(ROOT/'data'/b['instance']);mapping={int(v):i for i,v in enumerate(d['visits'])}
        initial=[mapping[v] for v in s['route']];assert validate(d,initial)==s['cost']
        order=METHODS[:];random.Random(b['instance']+str(b['seed'])).shuffle(order)
        for method in order:
            target=out/f'{path.stem}.{method}.{fingerprint}.json'
            if target.exists():continue
            row=run(d,initial,method,args.budget,b['seed'])
            endpoint=s['seconds']+row['total_seconds']
            prior=[s for s in b['snapshots'] if s['seconds']<=endpoint]
            row.update(instance=b['instance'],geometry=b['geometry'],seed=b['seed'],protocol=fingerprint,
                       initial_MDM_seconds=s['seconds'],recorded_MDM_cost_before_end=min(s['cost'] for s in prior),
                       baseline_final_cost=b['cost'],baseline_completed=b.get('completed',True))
            temporary=target.with_suffix('.partial');temporary.write_text(json.dumps(row,indent=2)+'\n');temporary.replace(target)
            print(json.dumps({k:v for k,v in row.items() if k not in ('trace','route')}),flush=True)

if __name__=='__main__':main()
