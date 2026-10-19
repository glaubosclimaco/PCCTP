"""Fixed-state control study; preserves the original Python control results."""
import json,sys,hashlib,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from experiments.lkh import support_search
from experiments.lkh.compiled_two_opt import two_opt,BIN
HERE=Path(__file__).resolve().parent
if __name__=='__main__':
    verification=json.loads((HERE/'two_opt_verification.json').read_text());assert verification['status']=='PASS'
    support_search.two_opt=two_opt
    out=HERE/'strong_state_support';out.mkdir(exist_ok=True)
    cases=json.loads((ROOT/'experiments/original_mdm/hard_training_baseline/summary.json').read_text())['rows']
    split=json.loads((ROOT/'experiments/ml_lns/decoder_initialized_states/protocol.json').read_text())
    assert len(cases)==9 and all(split['roles'][c['geometry']]=='train' for c in cases)
    protocol=dict(scope='fixed-state TRAIN; thirty additional seconds from fixed-work states; not end-to-end',source_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in ['strong_state_support.py','compiled_two_opt.py','support_search.py','two_opt.cpp']},binary_sha256=hashlib.sha256(Path(BIN).read_bytes()).hexdigest())
    target=out/'protocol.json'
    if target.exists():assert json.loads(target.read_text())==protocol
    else:target.write_text(json.dumps(protocol,indent=2)+'\n')
    for index,c in enumerate(cases):
        name=c['instance'];d=support_search.read_instance(ROOT/'data'/name);mapping={int(v):i for i,v in enumerate(d['visits'])}
        initial=json.loads((HERE/'fixed_work'/f'{name}.json').read_text());route=[mapping[v] for v in initial['lkh']['route']]
        modes=('fixed_lkh','support_2opt','support_lkh');modes=modes[index%3:]+modes[:index%3]
        for mode in modes:
            target=out/f'{name}.{mode}.json'
            if target.exists():r=json.loads(target.read_text())
            else:
                native_start=time.clock_gettime(time.CLOCK_MONOTONIC);cpu=time.process_time();r=support_search.run(d,route,mode,budget=30.);elapsed=time.clock_gettime(time.CLOCK_MONOTONIC)-native_start;raw=dict(r);late=elapsed>30 and r['cost']<support_search.validate(d,route)
                if late:r=dict(r,cost=support_search.validate(d,route),route=[int(d['visits'][v]) for v in route])
                r.update(instance=name,raw_result=raw,late_improvement_discarded=late,native_elapsed=elapsed,cpu_seconds=time.process_time()-cpu,source_status=initial['status']);target.write_text(json.dumps(r,indent=2)+'\n')
            print(json.dumps(dict(instance=name,mode=mode,cost=r['cost'],wall_seconds=r['wall_seconds'])),flush=True)
    (out/'completion.json').write_text(json.dumps(dict(runs=27))+'\n')
