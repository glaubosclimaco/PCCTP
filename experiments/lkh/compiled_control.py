"""Fixed-state control study; preserves the original Python control results."""
import json,sys,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from experiments.lkh import support_search
from experiments.lkh.compiled_two_opt import two_opt,BIN
HERE=Path(__file__).resolve().parent
if __name__=='__main__':
    verification=json.loads((HERE/'two_opt_verification.json').read_text());assert verification['status']=='PASS'
    support_search.two_opt=two_opt
    out=HERE/'compiled_control';out.mkdir(exist_ok=True)
    cases=json.loads((ROOT/'experiments/original_mdm/hard_training_baseline/summary.json').read_text())['rows']
    split=json.loads((ROOT/'experiments/ml_lns/decoder_initialized_states/protocol.json').read_text())
    assert len(cases)==9 and all(split['roles'][c['geometry']]=='train' for c in cases)
    protocol=dict(scope='fixed-state TRAIN; five additional seconds; not end-to-end',source_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in ['compiled_control.py','compiled_two_opt.py','support_search.py','two_opt.cpp']},binary_sha256=hashlib.sha256(Path(BIN).read_bytes()).hexdigest())
    target=out/'protocol.json'
    if target.exists():assert json.loads(target.read_text())==protocol
    else:target.write_text(json.dumps(protocol,indent=2)+'\n')
    for c in cases:
        name=c['instance'];d=support_search.read_instance(ROOT/'data'/name);mapping={int(v):i for i,v in enumerate(d['visits'])}
        initial=json.loads((HERE/'total_time'/f'{name}.seed42.mdm.json').read_text());route=[mapping[v] for v in initial['route']]
        for mode in ('fixed_lkh','support_2opt','support_lkh'):
            target=out/f'{name}.{mode}.json'
            if target.exists():r=json.loads(target.read_text())
            else:
                r=support_search.run(d,route,mode);r.update(instance=name);target.write_text(json.dumps(r,indent=2)+'\n')
            print(json.dumps(dict(instance=name,mode=mode,cost=r['cost'],wall_seconds=r['wall_seconds'])),flush=True)
    (out/'completion.json').write_text(json.dumps(dict(runs=27))+'\n')
