"""TRAIN equal-time native MDM, LKH, and random contracted exact repair."""
import json,time,sys,random,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'experiments/ml_lns'));sys.path.insert(0,str(ROOT))
from experiments.lkh import support_common_fallback as common
from experiments.lkh.screen import polish
from pilot import candidates,validate
from fragment_repair import FragmentRepair
HERE=Path(__file__).resolve().parent
old_search=common.search

def search(d,initial,mode,budget):
    if mode!='exact_random':return old_search(d,initial,mode,budget)
    began=time.perf_counter();end=began+budget;best=initial[:];cost=validate(d,best);mapping={int(v):i for i,v in enumerate(d['visits'])};trace=[];records=[]
    remaining=end-time.perf_counter()
    if remaining>0:
        r=polish(d,best,seconds=min(.2,remaining));candidate=[mapping[v] for v in r['route']]
        if time.perf_counter()<=end and r['cost']<cost:best,cost=candidate,r['cost']
    trace.append(dict(seconds=time.perf_counter()-began,cost=cost,route=[int(d['visits'][v]) for v in best]))
    model=FragmentRepair(d,backend='highspy',warm_start=True,connectivity='cuts');t=0
    while time.perf_counter()<end:
        pool=candidates(d,best,8,1,42+10007*t);t+=1
        if not pool:break
        remaining=end-time.perf_counter()
        if remaining<=0:break
        _,released=pool[0];r=model.solve(best,released,min(3.,remaining),operator='route_context',seed=42)
        elapsed=time.perf_counter()-began;records.append(dict(seconds=elapsed,released=[int(d['visits'][v]) for v in sorted(released)],result=r))
        if elapsed<=budget and r['cost']<cost:
            best=[mapping[v] for v in r['route']];cost=validate(d,best);trace.append(dict(seconds=elapsed,cost=cost,route=r['route']))
    return dict(mode=mode,initial_cost=validate(d,initial),cost=cost,route=[int(d['visits'][v]) for v in best],trace=trace,records=records,wall_seconds=time.perf_counter()-began)

if __name__=='__main__':
    out=HERE/'exact_random_total';out.mkdir(exist_ok=True);common.OUT=out;common.search=search
    cases=json.loads((ROOT/'experiments/original_mdm/hard_training_baseline/summary.json').read_text())['rows'];split=json.loads((ROOT/'experiments/ml_lns/decoder_initialized_states/protocol.json').read_text())
    assert len(cases)==9 and all(split['roles'][c['geometry']]=='train' for c in cases)
    modes=['mdm','fixed_lkh','exact_random'];protocol=dict(scope='exploratory TRAIN',cases=cases,seeds=[42,43,44],modes=modes,source_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in ['exact_random_total.py','support_common_fallback.py']})
    p=out/'protocol.json'
    if p.exists():assert json.loads(p.read_text())==protocol
    else:p.write_text(json.dumps(protocol,indent=2)+'\n')
    for c in cases:
        for seed in (42,43,44):
            order=modes[:];random.Random(c['instance']+str(seed)).shuffle(order)
            rows={m:common.run(c['instance'],seed,m) for m in order}
            print(json.dumps(dict(instance=c['instance'],seed=seed,costs={m:r['cost'] for m,r in rows.items()})),flush=True)
    (out/'completion.json').write_text(json.dumps(dict(runs=81))+'\n')
