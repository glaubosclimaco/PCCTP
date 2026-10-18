"""Complete development trajectories with matched extra wall budgets.

Learning is trained only on separate train geometries. Cheap-repair targets
are used, rather than silently treating bounded MIP rewards as identical.
ALNS is a simple transparent control, not a reproduction of a published solver.
"""
import argparse
import hashlib
import json
from pathlib import Path
import random
import time
import numpy as np
from pilot import candidates,features,free_edge_set,read_instance,validate
from context_study import route_features
from exchange_signature import signature
from selection_lift import best_lift
from fragment_repair import FragmentRepair
from nonlinear_ranker import fit,predict

HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[1]
GENERIC=[0,6,7,8,9,10,11]
METHODS=['uniform','joint_savings','alns','learned_generic','learned_joint',
         'hybrid_uniform','hybrid_savings','hybrid_alns','hybrid_generic','hybrid_joint']

def training_data():
    states=[]; hashes={}
    for path in sorted((HERE/'lifting_development').glob('*.json')):
        cheap=json.loads(path.read_text())
        source_path=HERE/'timed_oracle_decoder_initialized_states'/(
            cheap['instance']+f".seed{cheap.get('seed',42)}.phase{cheap['phase']}.k8.json")
        source=json.loads(source_path.read_text())
        if source['role']!='train': continue
        if hashlib.sha256(source_path.read_bytes()).hexdigest()!=cheap['source_sha256']:
            raise ValueError('Changed training provenance')
        for record,repair in zip(source['records'],cheap['records']):
            assert record['released']==repair['released']
            record['normalized_gain']=(source['initial_cost']-repair['lift_cost'])/max(1.,source['initial_cost'])
        states.append(source); hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
    if not states: raise ValueError('No independent train geometry labels')
    return states,hashes

def run(d,initial,method,models,seconds,seed):
    began=time.perf_counter(); route=initial[:]; value=validate(d,route)
    first=value; rng=random.Random(seed); pool_id=0; trace=[]; weights={k:1. for k in ('random','coverage','prize')}
    policy={'hybrid_uniform':'uniform','hybrid_savings':'joint_savings','hybrid_alns':'alns',
            'hybrid_generic':'learned_generic','hybrid_joint':'learned_joint'}.get(method,method)
    mip=FragmentRepair(d,backend='highspy',warm_start=True,connectivity='cuts') if method.startswith('hybrid_') else None
    while time.perf_counter()-began<seconds:
        pool=candidates(d,route,8,24,seed+10007*pool_id); pool_id+=1
        if not pool: break
        records=[]
        for kind,u in pool:
            record=dict(kind=kind,u=u)
            if policy not in ('uniform','alns'):
                free=free_edge_set(d,route,u,'route_context',seed=seed)
                f=features(d,route,u)+route_features(d,route,u,free)
                if policy!='learned_generic': f+=signature(d,route,u)
                record['features']=f
            records.append(record)
        if policy=='uniform': rng.shuffle(records)
        elif policy=='joint_savings': records.sort(key=lambda r:r['features'][20],reverse=True)
        elif policy not in ('alns',):
            model=models['generic' if policy=='learned_generic' else 'joint']
            scores=predict(model,records)
            records=[records[int(i)] for i in np.argsort(-scores,kind='stable')]
        while records and time.perf_counter()-began<seconds:
            if policy=='alns':
                kinds=sorted({r['kind'] for r in records})
                arm=rng.choices(kinds,weights=[weights[k] for k in kinds],k=1)[0]
                index=rng.choice([i for i,r in enumerate(records) if r['kind']==arm])
                record=records.pop(index)
            else: record=records.pop(0)
            query_began=time.perf_counter()
            repair=best_lift(d,route,record['u'],limit=3,seed=seed)
            next_route=repair['route']; next_cost=repair['cost']; mip_used=False
            remaining=seconds-(time.perf_counter()-began)
            if mip is not None and len(trace)%8==7 and remaining>.02:
                solved=mip.solve(route,record['u'],min(.2,remaining),'route_context',seed=seed,start_route=next_route)
                mapping={int(v):i for i,v in enumerate(d['visits'])}
                next_route=[mapping[v] for v in solved['route']]; next_cost=solved['cost']; mip_used=True
            improved=next_cost<value-1e-8
            weights[record['kind']]=.9*weights[record['kind']]+.1*(6. if improved else 1.)
            trace.append(dict(elapsed=time.perf_counter()-began,cost=next_cost,kind=record['kind'],
                              query_seconds=time.perf_counter()-query_began,mip_used=mip_used))
            if improved:
                route=next_route; value=next_cost; break
    assert abs(validate(d,route)-value)<1e-6
    return dict(method=method,budget_seconds=seconds,total_seconds=time.perf_counter()-began,
                initial_cost=first,cost=value,gain=(first-value)/max(1.,first),queries=len(trace),
                trace=trace,route=[int(d['visits'][v]) for v in route])

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--budget',type=float,default=5.)
    p.add_argument('--phase',type=int,default=1)
    p.add_argument('--limit',type=int,default=3)
    p.add_argument('--one-per-geometry',action='store_true')
    p.add_argument('--methods',nargs='+',choices=METHODS,default=METHODS)
    args=p.parse_args(); train,hashes=training_data()
    models=dict(generic=fit(train,GENERIC,'linear'),joint=fit(train,GENERIC+list(range(12,22)),'linear'))
    protocol=dict(scope='Development only; all online feature, choice and repair time included; not final benchmark',
                  code_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
                               ['hybrid_lns.py','selection_lift.py','fragment_repair.py','nonlinear_ranker.py',
                                'highs_backend.py','pilot.py','context_study.py','exchange_signature.py']},
                  training_hashes=hashes,training_geometries=sorted({s['geometry'] for s in train}),
                  models={name:{k:v.tolist() if isinstance(v,np.ndarray) else v for k,v in model.items()} for name,model in models.items()},
                  budget=args.budget,phase=args.phase,methods=args.methods,one_per_geometry=args.one_per_geometry)
    fingerprint=hashlib.sha256(json.dumps(protocol,sort_keys=True).encode()).hexdigest()[:12]
    out=HERE/'hybrid_development'; out.mkdir(exist_ok=True)
    (out/f'protocol.{fingerprint}.json').write_text(json.dumps(protocol,indent=2)+'\n')
    used=0; seen=set()
    for path in sorted((HERE/'decoder_initialized_states').glob('g_*.json')):
        baseline=json.loads(path.read_text())
        if baseline['role']!='development': continue
        if args.one_per_geometry and baseline['geometry'] in seen: continue
        if args.limit and used>=args.limit: break
        snapshot=next((s for s in baseline['snapshots'] if s['iteration']==args.phase),None)
        if snapshot is None: continue
        seen.add(baseline['geometry'])
        used+=1
        d=read_instance(ROOT/'data'/baseline['instance']); mapping={int(v):i for i,v in enumerate(d['visits'])}
        route=[mapping[v] for v in snapshot['route']]
        assert abs(validate(d,route)-snapshot['cost'])<1e-6
        order=args.methods[:]
        ordering_seed=int(hashlib.sha256((baseline['instance']+str(baseline['seed'])).encode()).hexdigest()[:12],16)
        random.Random(ordering_seed).shuffle(order)
        for method in order:
            target=out/(path.stem+f'.phase{args.phase}.{method}.{fingerprint}.json')
            if target.exists(): continue
            result=run(d,route,method,models,args.budget,baseline['seed'])
            endpoint=snapshot['seconds']+result['total_seconds']
            earlier=[s for s in baseline['snapshots'] if s['seconds']<=endpoint]
            result.update(instance=d['name'],geometry=baseline['geometry'],seed=baseline['seed'],
                          initial_phase=args.phase,model_fingerprint=fingerprint,initial_MDM_seconds=snapshot['seconds'],
                          recorded_MDM_cost_before_end=min(s['cost'] for s in earlier),
                          baseline_final_cost=baseline['cost'],baseline_completed=baseline.get('completed',True),
                          note='Sparse historical MDM checkpoints; no confirmed same-machine speed claim')
            target.write_text(json.dumps(result,indent=2)+'\n')
            print(json.dumps({k:v for k,v in result.items() if k not in ('trace','route')}),flush=True)

if __name__=='__main__': main()
