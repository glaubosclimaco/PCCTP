"""TRAIN-only four-arm equal-total-time evaluation."""
import sys,json,time,random,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from experiments.lkh.support_search import run as search
from experiments.original_mdm.timed_experiment import run_one
from experiments.ml_lns.pilot import read_instance,validate
HERE=Path(__file__).resolve().parent;OUT=HERE/'support_total'

def run(name,seed,mode):
    target=OUT/f'{name}.seed{seed}.{mode}.json'
    if target.exists():return json.loads(target.read_text())
    began=time.clock_gettime(time.CLOCK_MONOTONIC)
    native=run_one(name,seed,'original',15. if mode=='mdm' else 10.,OUT/f'{name}.seed{seed}.{mode}.native.json',True,parent_started=began)
    r=dict(instance=name,seed=seed,mode=mode,cost=native['cost'],route=native['route'],feasible=native['feasible_before_deadline'],post=None)
    if mode!='mdm' and r['feasible']:
        d=read_instance(ROOT/'data'/name);mapping={int(v):i for i,v in enumerate(d['visits'])}
        remaining=15.-(time.clock_gettime(time.CLOCK_MONOTONIC)-began)
        if remaining>0:
            post=search(d,[mapping[v] for v in native['route']],mode,remaining)
            r['post']=post
            if time.clock_gettime(time.CLOCK_MONOTONIC)-began<=15. and post['cost']<r['cost']:
                r['cost']=post['cost'];r['route']=post['route']
    r['wall_seconds']=time.clock_gettime(time.CLOCK_MONOTONIC)-began
    r['overrun_seconds']=max(0,r['wall_seconds']-15.)
    if r['feasible']:
        d=read_instance(ROOT/'data'/name);mapping={int(v):i for i,v in enumerate(d['visits'])}
        assert validate(d,[mapping[v] for v in r['route']])==r['cost']
    target.write_text(json.dumps(r,indent=2)+'\n');return r

if __name__=='__main__':
    OUT.mkdir(exist_ok=True)
    cases=json.loads((ROOT/'experiments/original_mdm/hard_training_baseline/summary.json').read_text())['rows']
    split=json.loads((ROOT/'experiments/ml_lns/decoder_initialized_states/protocol.json').read_text())
    assert len(cases)==9 and all(split['roles'][c['geometry']]=='train' for c in cases)
    modes=['mdm','fixed_lkh','support_2opt','support_lkh']
    protocol=dict(cases=cases,seeds=[42,43,44],modes=modes,total_budget=15,native_budget=10,source_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in ['support_total.py','support_search.py']})
    p=OUT/'protocol.json'
    if p.exists():assert json.loads(p.read_text())==protocol
    else:p.write_text(json.dumps(protocol,indent=2)+'\n')
    for c in cases:
        for seed in (42,43,44):
            order=modes[:];random.Random(c['instance']+str(seed)).shuffle(order)
            rows={m:run(c['instance'],seed,m) for m in order}
            print(json.dumps(dict(instance=c['instance'],seed=seed,costs={m:r['cost'] for m,r in rows.items()})),flush=True)
    (OUT/'completion.json').write_text(json.dumps(dict(runs=108))+'\n')
