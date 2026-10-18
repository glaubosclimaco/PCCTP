"""TRAIN fixed-state complementarity diagnostic, not a full algorithm comparison."""
import sys,json,time,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from experiments.ml_lns.pilot import read_instance,validate,route_cost
from experiments.lkh.screen import polish
HERE=Path(__file__).resolve().parent

def two_opt(d,route,end):
    route=route[:];best=route_cost(d,route)
    while time.perf_counter()<end:
        changed=False
        for i in range(len(route)-1):
            for j in range(i+2,len(route)):
                if time.perf_counter()>=end:return route
                candidate=route[:i+1]+list(reversed(route[i+1:j+1]))+route[j+1:]
                value=route_cost(d,candidate)
                if value<best:
                    route,best,changed=candidate,value,True;break
            if changed:break
        if not changed:break
    return route

def run(d,initial,mode,budget=5.):
    began=time.perf_counter();end=began+budget;route=initial[:];best=validate(d,route);trace=[];evaluated=0
    mapping={int(v):i for i,v in enumerate(d['visits'])}
    def order(candidate):
        remaining=end-time.perf_counter()
        if remaining<=0:return candidate
        if mode=='support_2opt':return two_opt(d,candidate,end)
        result=polish(d,candidate,seconds=min(.2,remaining))
        return [mapping[v] for v in result['route']]
    candidate=order(route)
    value=validate(d,candidate)
    if time.perf_counter()<=end and value<best:route,best=candidate,value
    trace.append(dict(seconds=time.perf_counter()-began,cost=best,route=[int(d['visits'][v]) for v in route]))
    while mode!='fixed_lkh' and time.perf_counter()<end:
        pool=[];selected=set(route);outside=sorted(set(range(len(d['visits'])))-selected)
        optional=sorted(selected-d['mandatory'])
        def offer(candidate):
            try:cost=validate(d,candidate)
            except ValueError:return
            pool.append((cost,candidate))
        for v in optional:
            if time.perf_counter()>=end:break
            offer([u for u in route if u!=v])
            for w in outside:
                if time.perf_counter()>=end:break
                offer([w if u==v else u for u in route])
        for w in outside:
            if time.perf_counter()>=end:break
            pos=min(range(len(route)),key=lambda i:d['costs'][route[i-1],w]+d['costs'][w,route[i]]-d['costs'][route[i-1],route[i]])
            offer(route[:pos]+[w]+route[pos:])
        improved=False
        for _,candidate in sorted(pool,key=lambda p:(p[0],p[1]))[:8]:
            if time.perf_counter()>=end:break
            candidate=order(candidate);evaluated+=1;value=validate(d,candidate)
            if time.perf_counter()<=end and value<best:
                route,best,improved=candidate,value,True
                trace.append(dict(seconds=time.perf_counter()-began,cost=best,route=[int(d['visits'][v]) for v in route]))
                break
        if not improved:break
    return dict(mode=mode,cost=best,route=[int(d['visits'][v]) for v in route],initial_cost=validate(d,initial),wall_seconds=time.perf_counter()-began,trace=trace,evaluated=evaluated)

if __name__=='__main__':
    out=HERE/'support_diagnostic';out.mkdir(exist_ok=True)
    cases=json.loads((ROOT/'experiments/original_mdm/hard_training_baseline/summary.json').read_text())['rows']
    split=json.loads((ROOT/'experiments/ml_lns/decoder_initialized_states/protocol.json').read_text())
    assert all(split['roles'][c['geometry']]=='train' for c in cases)
    source_hash=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    rows=[]
    for case in cases:
        name=case['instance'];d=read_instance(ROOT/'data'/name);mapping={int(v):i for i,v in enumerate(d['visits'])}
        initial=json.loads((HERE/'total_time'/f'{name}.seed42.mdm.json').read_text())
        route=[mapping[v] for v in initial['route']]
        for mode in ('fixed_lkh','support_2opt','support_lkh'):
            target=out/f'{name}.{mode}.json'
            if target.exists():r=json.loads(target.read_text());assert r['source_hash']==source_hash
            else:
                r=run(d,route,mode);r.update(instance=name,source_hash=source_hash)
                target.write_text(json.dumps(r,indent=2)+'\n')
            rows.append(r);print(json.dumps({k:v for k,v in r.items() if k not in ('route','trace')}),flush=True)
    (out/'completion.json').write_text(json.dumps(dict(runs=len(rows),source_hash=source_hash))+'\n')
