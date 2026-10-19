"""Offline TRAIN candidate headroom, before training any new model."""
import json,random,sys,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from experiments.ml_lns.pilot import read_instance,validate
from experiments.lkh.screen import polish
HERE=Path(__file__).resolve().parent
if __name__=='__main__':
    out=HERE/'candidate_oracle';out.mkdir(exist_ok=True)
    cases=json.loads((ROOT/'experiments/original_mdm/hard_training_baseline/summary.json').read_text())['rows'];split=json.loads((ROOT/'experiments/ml_lns/decoder_initialized_states/protocol.json').read_text())
    assert len(cases)==9 and all(split['roles'][c['geometry']]=='train' for c in cases)
    for case in cases:
        name=case['instance'];target=out/f'{name}.json'
        if target.exists():continue
        d=read_instance(ROOT/'data'/name);mapping={int(v):i for i,v in enumerate(d['visits'])};initial_path=HERE/'compiled_control'/f'{name}.fixed_lkh.json';initial=json.loads(initial_path.read_text());route=[mapping[v] for v in initial['route']]
        pool={};selected=set(route);outside=sorted(set(range(len(d['visits'])))-selected)
        def offer(candidate):
            try:cost=validate(d,candidate)
            except ValueError:return
            key=tuple(sorted(candidate))
            if key not in pool or cost<pool[key][0]:pool[key]=(cost,candidate)
        for v in sorted(selected-d['mandatory']):
            offer([u for u in route if u!=v])
            for w in outside:offer([w if u==v else u for u in route])
        for w in outside:
            i=min(range(len(route)),key=lambda i:d['costs'][route[i-1],w]+d['costs'][w,route[i]]-d['costs'][route[i-1],route[i]])
            offer(route[:i]+[w]+route[i:])
        candidates=sorted(pool.values(),key=lambda p:(p[0],p[1]));indices=list(range(min(8,len(candidates))))
        remaining=list(range(8,len(candidates)));random.Random(42).shuffle(remaining);indices+=remaining[:40];records=[]
        for rank in indices:
            cost,candidate=candidates[rank];r=polish(d,candidate,seconds=.2)
            records.append(dict(raw_rank=rank,raw_cost=cost,raw_route=[int(d['visits'][v]) for v in candidate],result=r))
        target.write_text(json.dumps(dict(instance=name,initial_cost=initial['cost'],initial_route=initial['route'],initial_sha256=hashlib.sha256(initial_path.read_bytes()).hexdigest(),pool_size=len(candidates),records=records),indent=2)+'\n')
        print(json.dumps(dict(instance=name,pool_size=len(candidates),evaluated=len(records))),flush=True)
    (out/'completion.json').write_text(json.dumps(dict(cases=9))+'\n')
