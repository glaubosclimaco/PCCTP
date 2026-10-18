"""Frozen seven-geometry validation, gated on completed TRAIN audit only."""
import json,random,hashlib
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from experiments.lkh import total_time
HERE=Path(__file__).resolve().parent


if __name__=='__main__':
    audit=json.loads((HERE/'total_time_audit.json').read_text())
    assert audit['status']=='PASS' and audit['gate_passed']
    split_path=ROOT/'experiments/ml_lns/decoder_initialized_states/protocol.json'
    split=json.loads(split_path.read_text())
    cases=[dict(geometry=g,instance=split['instances'][g]) for g,role in sorted(split['roles'].items()) if role=='validation']
    assert len(cases)==7
    out=HERE/'validation';out.mkdir(exist_ok=True);total_time.OUT=out
    protocol=dict(cases=cases,seeds=[42,43,44],arms=['mdm','hybrid'],allocation=[9.,1.],total_budget=10.,
                  original_split_sha256=hashlib.sha256(split_path.read_bytes()).hexdigest(),
                  train_audit_sha256=hashlib.sha256((HERE/'total_time_audit.json').read_bytes()).hexdigest(),
                  test_opened=False)
    p=out/'protocol.json'
    if p.exists():assert json.loads(p.read_text())==protocol
    else:p.write_text(json.dumps(protocol,indent=2)+'\n')
    for case in cases:
        assert split['roles'][case['geometry']]=='validation'
        for seed in (42,43,44):
            modes=['mdm','hybrid'];random.Random(case['instance']+str(seed)).shuffle(modes)
            pair={m:total_time.run(case['instance'],seed,m) for m in modes}
            print(json.dumps(dict(instance=case['instance'],seed=seed,mdm_cost=pair['mdm']['cost'],hybrid_cost=pair['hybrid']['cost'])),flush=True)
    (out/'completion.json').write_text(json.dumps(dict(cases=7,seeds=[42,43,44],arms=['mdm','hybrid']))+'\n')
