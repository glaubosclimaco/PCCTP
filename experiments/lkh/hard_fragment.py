"""Offline contracted exact-neighborhood headroom on hard TRAIN only."""
import json,sys,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'experiments/ml_lns'));sys.path.insert(0,str(ROOT))
from pilot import read_instance,candidates
from fragment_repair import FragmentRepair
HERE=Path(__file__).resolve().parent
if __name__=='__main__':
    out=HERE/'hard_fragment';out.mkdir(exist_ok=True)
    cases=json.loads((ROOT/'experiments/original_mdm/hard_training_baseline/summary.json').read_text())['rows'];split=json.loads((ROOT/'experiments/ml_lns/decoder_initialized_states/protocol.json').read_text())
    assert len(cases)==9 and all(split['roles'][c['geometry']]=='train' for c in cases)
    protocol=dict(k=8,count=9,seed=42,seconds=3,backend='highspy',warm_start=True,connectivity='cuts',operator='route_context',scope='offline TRAIN only; local proof is not global',source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    p=out/'protocol.json'
    if p.exists():assert json.loads(p.read_text())==protocol
    else:p.write_text(json.dumps(protocol,indent=2)+'\n')
    for c in cases:
        name=c['instance'];d=read_instance(ROOT/'data'/name);mapping={int(v):i for i,v in enumerate(d['visits'])};initial=json.loads((HERE/'compiled_control'/f'{name}.fixed_lkh.json').read_text());route=[mapping[v] for v in initial['route']]
        model=FragmentRepair(d,backend='highspy',warm_start=True,connectivity='cuts')
        for i,(kind,released) in enumerate(candidates(d,route,8,9,42)):
            p=out/f'{name}.neighbor{i}.json'
            if p.exists():continue
            result=model.solve(route,released,3.,operator='route_context',seed=42)
            r=dict(instance=name,index=i,kind=kind,released=[int(d['visits'][v]) for v in sorted(released)],initial_cost=initial['cost'],initial_route=initial['route'],result=result)
            p.write_text(json.dumps(r,indent=2)+'\n')
            print(json.dumps(dict(instance=name,index=i,kind=kind,cost=result['cost'],status=result['status'])),flush=True)
    (out/'completion.json').write_text(json.dumps(dict(cases=9))+'\n')
