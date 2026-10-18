"""Development-only end-to-end pilot from independently validated original MDM states.

Includes all online feature, neighborhood and model construction time. Training
is offline and uses train geometries only. This is not the final benchmark.
"""
import argparse
import hashlib
import json
from pathlib import Path
import random
import time
import numpy as np
from pilot import Repair,candidates,features,free_edge_set,read_instance,validate
from context_study import route_features
from exchange_signature import signature

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
GENERIC=[0,6,7,8,9,10,11]

def train_ranker(states,indices):
    x=np.array([r['features'] for s in states for r in s['records']])[:,indices]
    mean=x.mean(axis=0); scale=x.std(axis=0); scale[scale<1e-10]=1.
    xx=[]; yy=[]
    for state in states:
        records=state['records']
        for i,r in enumerate(records):
            for other in records[:i]:
                delta=r['normalized_gain']-other['normalized_gain']
                if abs(delta)>1e-12:
                    xx.append((np.array(r['features'])[indices]-np.array(other['features'])[indices])/scale)
                    yy.append(delta)
    if xx:
        x=np.array(xx); y=np.array(yy)
        weights=np.linalg.solve(x.T@x+np.eye(len(indices)),x.T@y)
    else: weights=np.zeros(len(indices))
    return dict(indices=indices,mean=mean.tolist(),scale=scale.tolist(),weights=weights.tolist())

def run(d,initial,method,ranker,budget,seed,k=8,count=24):
    began=time.perf_counter(); route=initial[:]; initial_cost=validate(d,route)
    model=Repair(d); trace=[]; rng=random.Random(seed); pools=0; feature_seconds=0.
    while time.perf_counter()-began<budget:
        pool=candidates(d,route,k,count,seed+pools*10007); pools+=1
        if not pool: break
        scored=[]; feature_began=time.perf_counter()
        for kind,u in pool:
            if method=='uniform': score=rng.random()
            else:
                free=free_edge_set(d,route,u,'route_context',seed=seed)
                f=features(d,route,u)+route_features(d,route,u,free)
                if method in ('joint_ridge','joint_savings'): f+=signature(d,route,u)
                if method=='joint_savings': score=f[20]
                else:
                    indices=ranker['indices']
                    score=float(((np.array(f)[indices]-ranker['mean'])/ranker['scale'])@ranker['weights'])
            scored.append((score,kind,u))
        feature_seconds+=time.perf_counter()-feature_began
        for _,kind,u in sorted(scored,key=lambda row:row[0],reverse=True):
            remaining=budget-(time.perf_counter()-began)
            if remaining<.01: break
            result=model.solve(route,u,min(.5,remaining),operator='route_context',seed=seed)
            trace.append(dict(elapsed=time.perf_counter()-began,kind=kind,cost=result['cost'],
                              status=result['status'],proven_optimal=result['proven_optimal']))
            mapping={int(v):i for i,v in enumerate(d['visits'])}
            next_route=[mapping[v] for v in result['route']]
            if result['cost']<validate(d,route)-1e-8:
                route=next_route
                break  # Scores refer to the old incumbent; recompute all of them.
    final_cost=validate(d,route)
    return dict(method=method,budget_seconds=budget,total_seconds=time.perf_counter()-began,
                feature_seconds=feature_seconds,initial_cost=initial_cost,cost=final_cost,
                normalized_gain=(initial_cost-final_cost)/max(1.,initial_cost),
                queries=len(trace),route=[int(d['visits'][v]) for v in route],trace=trace)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--budget',type=float,default=10.)
    parser.add_argument('--phase',type=int,default=70)
    parser.add_argument('--limit',type=int,default=3)
    args=parser.parse_args()
    train=[]; manifest={}
    for path in sorted((HERE/'strong_oracle').glob('*.json')):
        row=json.loads(path.read_text())
        if row['role']=='train':
            train.append(row); manifest[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
    if not train: raise ValueError('No independent training states')
    models=dict(generic_ridge=train_ranker(train,GENERIC),
                joint_ridge=train_ranker(train,GENERIC+list(range(12,22))))
    output=HERE/'iterative_development'; output.mkdir(exist_ok=True)
    freeze=dict(scope='Development only; no validation/test labels; offline fit fixed before trajectory',
                source_hashes=manifest,models=models,phase=args.phase,budget=args.budget)
    fingerprint=hashlib.sha256(json.dumps(freeze,sort_keys=True).encode()).hexdigest()[:12]
    (output/f'models.{fingerprint}.json').write_text(json.dumps(freeze,indent=2)+'\n')
    used=0
    for path in sorted((HERE/'original_states').glob('g_*.json')):
        row=json.loads(path.read_text())
        if row.get('role')!='development': continue
        snapshot=next((s for s in row['snapshots'] if s['iteration']==args.phase),None)
        if snapshot is None: continue
        if args.limit and used>=args.limit: break
        used+=1
        d=read_instance(ROOT/'data'/row['instance'])
        mapping={int(v):i for i,v in enumerate(d['visits'])}
        initial=[mapping[v] for v in snapshot['route']]
        assert abs(validate(d,initial)-snapshot['cost'])<1e-6
        for method in ('uniform','joint_savings','generic_ridge','joint_ridge'):
            target=output/f'{path.stem}.phase{args.phase}.{method}.{fingerprint}.json'
            if target.exists(): continue
            result=run(d,initial,method,models.get(method),args.budget,row['seed'])
            result.update(instance=d['name'],geometry=row['geometry'],seed=row['seed'],
                          initial_phase=args.phase,model_fingerprint=fingerprint,role='development')
            target.write_text(json.dumps(result,indent=2)+'\n')
            print(json.dumps({k:v for k,v in result.items() if k not in ('trace','route')}),flush=True)

if __name__=='__main__': main()
