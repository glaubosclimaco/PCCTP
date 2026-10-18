"""Development diagnostic of neighborhood radius after original MDM search.

Uses exact path contraction. No training or confirmatory inference. Radius
changes are logged explicitly rather than silently tuning a held-out test.
"""
import argparse
import json
from pathlib import Path
import time
from pilot import candidates,read_instance,validate
from fragment_repair import FragmentRepair

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--geometries',nargs='+',default=['g_gr137','g_pr144','g_u159'])
    p.add_argument('--radii',nargs='+',type=int,default=[8,16,24])
    p.add_argument('--count',type=int,default=12)
    p.add_argument('--seconds',type=float,default=.5)
    args=p.parse_args()
    output=HERE/'radius_development'; output.mkdir(exist_ok=True)
    for path in sorted((HERE/'original_states').glob('g_*.json')):
        baseline=json.loads(path.read_text())
        if baseline['geometry'] not in args.geometries: continue
        if baseline['role'] not in ('development','train'): raise ValueError('Reserved geometry')
        snapshot=next((s for s in baseline['snapshots'] if s['iteration']==70),None)
        if not snapshot: continue
        d=read_instance(ROOT/'data'/baseline['instance'])
        mapping={int(v):i for i,v in enumerate(d['visits'])}
        route=[mapping[v] for v in snapshot['route']]; initial=validate(d,route)
        model=FragmentRepair(d)
        for k in args.radii:
            target=output/f'{path.stem}.k{k}.count{args.count}.s{args.seconds}.json'
            if target.exists(): continue
            began=time.perf_counter(); records=[]
            for kind,u in candidates(d,route,k,args.count,baseline['seed']):
                result=model.solve(route,u,args.seconds,operator='route_context',seed=baseline['seed'])
                records.append(dict(kind=kind,released=sorted(u),**result))
            row=dict(instance=d['name'],geometry=baseline['geometry'],role=baseline['role'],
                     phase=70,k=k,initial_cost=initial,seconds_per_candidate=args.seconds,
                     elapsed=time.perf_counter()-began,records=records,
                     best_gain=max((initial-r['cost'])/initial for r in records),
                     improving=sum(r['cost']<initial-1e-8 for r in records),
                     optimal=sum(r['proven_optimal'] for r in records),queries=len(records))
            target.write_text(json.dumps(row,indent=2)+'\n')
            print(json.dumps({k:v for k,v in row.items() if k!='records'}),flush=True)

if __name__=='__main__': main()
