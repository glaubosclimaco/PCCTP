"""Collect fixed-state repair labels from the supplied original MDM snapshots.

Never train on test geometry labels. Cases are resumable and have unique keys.
Current development budget k=8,24 candidates,.5s was frozen before new labels.
"""
import argparse
import json
from pathlib import Path
import time
from pilot import FEATURES, Repair, candidates, features, free_edge_set, read_instance, validate
from context_study import ROUTE_FEATURES, route_features
from exchange_signature import NAMES, signature

ROOT=Path(__file__).resolve().parents[2]
BASE=Path(__file__).with_name('original_states')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--roles',nargs='+',default=['development','train'])
    parser.add_argument('--phases',nargs='+',type=int,default=[1,10,70])
    parser.add_argument('--seeds',nargs='+',type=int,default=[42])
    parser.add_argument('--k',type=int,default=8)
    parser.add_argument('--candidates',type=int,default=24)
    parser.add_argument('--seconds',type=float,default=.5)
    parser.add_argument('--state-dir',type=Path,default=BASE)
    parser.add_argument('--output-dir',type=Path)
    args=parser.parse_args()
    if args.output_dir is None:
        args.output_dir=Path(__file__).with_name('strong_oracle' if args.state_dir.resolve()==BASE.resolve()
                                               else 'strong_oracle_'+args.state_dir.name)
    plan=json.loads((args.state_dir/'protocol.json').read_text())
    args.output_dir.mkdir(exist_ok=True)
    failures=[]
    for path in sorted(args.state_dir.glob('g_*.json')):
        try:
            baseline=json.loads(path.read_text())
        except json.JSONDecodeError:
            print(json.dumps(dict(skipped_incomplete_export=path.name)),flush=True)
            continue
        if not baseline.get('independently_validated') or 'geometry' not in baseline:
            print(json.dumps(dict(skipped_unvalidated_export=path.name)),flush=True)
            continue
        geometry=baseline['geometry']
        if plan['roles'][geometry] not in args.roles: continue
        if baseline['seed'] not in args.seeds: continue
        d=read_instance(ROOT/'data'/baseline['instance'])
        mapping={int(v):i for i,v in enumerate(d['visits'])}
        seen=set(); model=None
        for snapshot in baseline['snapshots']:
            phase=snapshot['iteration']
            if phase not in args.phases: continue
            route=[mapping[v] for v in snapshot['route']]
            # Do not count identical best-so-far states at successive phases twice.
            statekey=tuple(sorted(tuple(sorted((route[i-1],v))) for i,v in enumerate(route)))
            if statekey in seen: continue
            seen.add(statekey)
            target=args.output_dir/(path.stem+f'.phase{phase}.k{args.k}.json')
            if target.exists(): continue
            began=time.perf_counter(); initial=validate(d,route)
            if abs(initial-snapshot['cost'])>1e-6: raise ValueError('Invalid baseline snapshot')
            if model is None: model=Repair(d)
            row=dict(instance=d['name'],geometry=geometry,role=plan['roles'][geometry],
                     seed=baseline['seed'],phase=phase,initial_cost=initial,initial_route=snapshot['route'],
                     baseline_method=baseline['method'],k=args.k,seconds_per_candidate=args.seconds,
                     feature_names=FEATURES+ROUTE_FEATURES+NAMES,generic_feature_indices=[0,6,7,8,9,10,11],
                     records=[])
            for kind,released in candidates(d,route,args.k,args.candidates,baseline['seed']):
                free=free_edge_set(d,route,released,'route_context',seed=baseline['seed'])
                f=features(d,route,released)+route_features(d,route,released,free)+signature(d,route,released)
                result=model.solve(route,released,args.seconds,operator='route_context',seed=baseline['seed'])
                row['records'].append(dict(kind=kind,released=[int(d['visits'][v]) for v in sorted(released)],
                    features=f,normalized_gain=(initial-result['cost'])/max(1.,initial),**result))
                row['elapsed_seconds']=time.perf_counter()-began
                # .partial is never treated as a complete case by evaluators.
                target.with_suffix('.partial').write_text(json.dumps(row,indent=2)+'\n')
            target.write_text(json.dumps(row,indent=2)+'\n')
            target.with_suffix('.partial').unlink(missing_ok=True)
            print(json.dumps(dict(instance=d['name'],role=row['role'],phase=phase,queries=len(row['records']),
                improvements=sum(r['normalized_gain']>1e-9 for r in row['records']),
                best_gain=max((r['normalized_gain'] for r in row['records']),default=0.),seconds=row['elapsed_seconds'])),flush=True)


if __name__=='__main__': main()
