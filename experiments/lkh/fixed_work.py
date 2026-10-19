"""Longer fixed-iteration baseline diagnostic, preserving censored cases."""
import sys,json,time,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from experiments.original_mdm.timed_experiment import run_one
from experiments.ml_lns.pilot import read_instance,validate
from experiments.lkh.screen import polish
HERE=Path(__file__).resolve().parent
if __name__=='__main__':
    out=HERE/'fixed_work';out.mkdir(exist_ok=True)
    cases=json.loads((ROOT/'experiments/original_mdm/hard_training_baseline/summary.json').read_text())['rows'];split=json.loads((ROOT/'experiments/ml_lns/decoder_initialized_states/protocol.json').read_text())
    assert len(cases)==9 and all(split['roles'][c['geometry']]=='train' for c in cases)
    build=ROOT/'experiments/original_mdm/build_timed_elite_bounded'
    protocol=dict(scope='TRAIN fixed70-iteration diagnostic, additional LKH time',seed=42,safety_timeout=1800,binary_sha256=hashlib.sha256((build/'original_mdm').read_bytes()).hexdigest(),build_manifest=json.loads((build/'manifest.json').read_text()),cases=cases)
    p=out/'protocol.json'
    if p.exists():assert json.loads(p.read_text())==protocol
    else:p.write_text(json.dumps(protocol,indent=2)+'\n')
    for c in cases:
        name=c['instance'];p=out/f'{name}.json'
        if p.exists():continue
        began=time.clock_gettime(time.CLOCK_MONOTONIC)
        native=run_one(name,42,'original',0.,out/f'{name}.native.json',True,parent_started=began,hard_timeout_seconds=1800)
        raw=out/f'{name}.native.raw';checkpoint=Path(str(raw)+'.checkpoint');payload=None
        for candidate in (raw,checkpoint):
            if candidate.exists():
                try:payload=json.loads(candidate.read_text());break
                except json.JSONDecodeError:pass
        if payload is None:
            r=dict(instance=name,status='NO_VALID_EXPORT',native=native)
        else:
            d=read_instance(ROOT/'data'/name);mapping={int(v):i for i,v in enumerate(d['visits'])};route=[mapping[v] for v in payload['route']];cost=validate(d,route);assert cost==payload['cost']
            completed=bool(payload['completed']) and not native['hard_timeout'];last=max((s['iteration'] for s in payload['snapshots']),default=0)
            if completed:assert last==70
            lkh=polish(d,route,seconds=2.)
            r=dict(instance=name,status='COMPLETED70' if completed else 'CENSORED',last_iteration=last,native_cost=cost,native_route=payload['route'],native=native,lkh=lkh,published_mdm_best=c['published']['best'])
        p.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items() if k not in ('native','lkh','native_route')}),flush=True)
    (out/'completion.json').write_text(json.dumps(dict(cases=9))+'\n')
