"""Exploratory 30-second training persistence; never unlocks test data."""
import sys,json,random,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from experiments.lkh import total_time
if __name__=='__main__':
    out=Path(__file__).with_name('horizon30');out.mkdir(exist_ok=True);total_time.OUT=out
    cases=json.loads((ROOT/'experiments/original_mdm/hard_training_baseline/summary.json').read_text())['rows']
    split=json.loads((ROOT/'experiments/ml_lns/decoder_initialized_states/protocol.json').read_text())
    assert len(cases)==9 and all(split['roles'][c['geometry']]=='train' for c in cases)
    protocol=dict(cases=cases,seeds=[42,43,44],total_budget=30.,polish_budget=1.,scope='exploratory TRAIN only; does not unlock test',source_sha256=hashlib.sha256(Path(total_time.__file__).read_bytes()).hexdigest())
    target=out/'protocol.json'
    if target.exists():assert json.loads(target.read_text())==protocol
    else:target.write_text(json.dumps(protocol,indent=2)+'\n')
    for case in cases:
        for seed in (42,43,44):
            order=['mdm','hybrid'];random.Random(case['instance']+str(seed)).shuffle(order)
            pair={m:total_time.run(case['instance'],seed,m,30.,1.) for m in order}
            print(json.dumps(dict(instance=case['instance'],seed=seed,mdm_cost=pair['mdm']['cost'],hybrid_cost=pair['hybrid']['cost'])),flush=True)
    (out/'completion.json').write_text(json.dumps(dict(cases=9,seeds=[42,43,44]))+'\n')
