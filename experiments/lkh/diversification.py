"""Exploratory TRAIN diversification with preserved best feasible incumbent."""
import json,time,random,hashlib,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from experiments.lkh import support_search
from experiments.lkh.compiled_two_opt import two_opt
from experiments.lkh.feasible_fallback import construct
HERE=Path(__file__).resolve().parent

def run(name,mode):
    began=time.perf_counter();end=began+15;d=support_search.read_instance(ROOT/'data'/name);f=construct(d);best=f['route'];cost=support_search.validate(d,best)
    rng=random.Random(42);trace=[];iterations=0
    if time.perf_counter()<=end:trace.append(dict(seconds=time.perf_counter()-began,cost=cost,route=[int(d['visits'][v]) for v in best]))
    while time.perf_counter()<end:
        candidate=best[:]
        if iterations:
            i,j=sorted(rng.sample(range(len(candidate)),2));candidate[i:j+1]=reversed(candidate[i:j+1])
            if mode!='fixed_lkh':
                for _ in range(2):
                    for attempt in range(20):
                        if time.perf_counter()>=end:break
                        optional=sorted(set(candidate)-d['mandatory']);outside=sorted(set(range(len(d['visits'])))-set(candidate))
                        if not optional or not outside:break
                        old=rng.choice(optional);new=rng.choice(outside);trial=[new if v==old else v for v in candidate]
                        try:support_search.validate(d,trial)
                        except ValueError:continue
                        candidate=trial;break
        remaining=end-time.perf_counter()
        if remaining<=0:break
        r=support_search.run(d,candidate,mode,min(1.,remaining));iterations+=1
        mapping={int(v):i for i,v in enumerate(d['visits'])};route=[mapping[v] for v in r['route']]
        value=support_search.validate(d,route)
        if time.perf_counter()<=end and value<cost:
            best,cost=route,value;trace.append(dict(seconds=time.perf_counter()-began,cost=cost,route=r['route']))
    return dict(instance=name,mode=mode,cost=cost,route=[int(d['visits'][v]) for v in best],trace=trace,iterations=iterations,wall_seconds=time.perf_counter()-began)

if __name__=='__main__':
    support_search.two_opt=two_opt;out=HERE/'diversification';out.mkdir(exist_ok=True)
    cases=json.loads((ROOT/'experiments/original_mdm/hard_training_baseline/summary.json').read_text())['rows']
    split=json.loads((ROOT/'experiments/ml_lns/decoder_initialized_states/protocol.json').read_text())
    assert len(cases)==9 and all(split['roles'][c['geometry']]=='train' for c in cases)
    protocol=dict(scope='TRAIN exploratory, not confirmatory',seed=42,budget=15,source_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in ['diversification.py','support_search.py','compiled_two_opt.py','feasible_fallback.py']})
    target=out/'protocol.json'
    if target.exists():assert json.loads(target.read_text())==protocol
    else:target.write_text(json.dumps(protocol,indent=2)+'\n')
    for c in cases:
        for mode in ('fixed_lkh','support_2opt','support_lkh'):
            target=out/f"{c['instance']}.{mode}.json"
            if target.exists():r=json.loads(target.read_text())
            else:r=run(c['instance'],mode);target.write_text(json.dumps(r,indent=2)+'\n')
            print(json.dumps({k:v for k,v in r.items() if k not in ('route','trace')}),flush=True)
    (out/'completion.json').write_text(json.dumps(dict(runs=27))+'\n')
