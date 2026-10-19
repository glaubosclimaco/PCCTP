"""TRAIN viability screen from cheap feasible construction, without native MDM."""
import json,time,sys,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from experiments.lkh.feasible_fallback import construct
from experiments.lkh.support_search import run
from experiments.ml_lns.pilot import read_instance,validate
HERE=Path(__file__).resolve().parent
if __name__=='__main__':
    out=HERE/'direct_support';out.mkdir(exist_ok=True)
    cases=json.loads((ROOT/'experiments/original_mdm/hard_training_baseline/summary.json').read_text())['rows']
    split=json.loads((ROOT/'experiments/ml_lns/decoder_initialized_states/protocol.json').read_text())
    assert len(cases)==9 and all(split['roles'][c['geometry']]=='train' for c in cases)
    protocol=dict(scope='TRAIN viability only; cached MDM comparison not confirmatory',budget=15,source_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in ['direct_support.py','support_search.py','feasible_fallback.py']})
    target=out/'protocol.json'
    if target.exists():assert json.loads(target.read_text())==protocol
    else:target.write_text(json.dumps(protocol,indent=2)+'\n')
    for c in cases:
        name=c['instance'];target=out/f'{name}.json'
        if target.exists():r=json.loads(target.read_text())
        else:
            began=time.perf_counter();d=read_instance(ROOT/'data'/name);f=construct(d);elapsed=time.perf_counter()-began
            route=[int(d['visits'][v]) for v in f['route']];trace=[]
            if elapsed<=15:trace.append(dict(seconds=elapsed,cost=f['cost'],route=route))
            remaining=15-(time.perf_counter()-began);post=None
            if remaining>0:
                post_start=time.perf_counter()-began;post=run(d,f['route'],'support_lkh',remaining)
                trace.extend(dict(e,seconds=e['seconds']+post_start) for e in post['trace'] if e['seconds']+post_start<=15)
            best=min(trace,key=lambda e:e['cost']) if trace else None
            r=dict(instance=name,cost=best['cost'] if best else None,route=best['route'] if best else [],trace=trace,fallback_cost=f['cost'],post=post,wall_seconds=time.perf_counter()-began)
            target.write_text(json.dumps(r,indent=2)+'\n')
        print(json.dumps({k:r[k] for k in ['instance','cost','fallback_cost','wall_seconds']}),flush=True)
    (out/'completion.json').write_text(json.dumps(dict(cases=9))+'\n')
