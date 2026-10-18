"""Paired fixed 9+1 allocation versus 10-second MDM on hard TRAIN only."""
import json,time,random
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from experiments.original_mdm.timed_experiment import run_one
from experiments.ml_lns.pilot import read_instance,validate
from experiments.lkh.screen import polish
OUT=Path(__file__).with_name('total_time')


def run(name,seed,mode):
    target=OUT/f'{name}.seed{seed}.{mode}.json'
    if target.exists():return json.loads(target.read_text())
    began=time.clock_gettime(time.CLOCK_MONOTONIC)
    native=run_one(name,seed,'original',10. if mode=='mdm' else 9.,
                   OUT/f'{name}.seed{seed}.{mode}.native.json',True,parent_started=began)
    result=dict(instance=name,seed=seed,mode=mode,cost=native['cost'],route=native['route'],
                feasible=native['feasible_before_deadline'],lkh_used=False,native_overrun=native['overrun_seconds'],
                hard_timeout=native['hard_timeout'])
    if mode=='hybrid' and result['feasible']:
        d=read_instance(ROOT/'data'/name);mapping={int(v):i for i,v in enumerate(d['visits'])}
        route=[mapping[v] for v in native['route']]
        remaining=10.-(time.clock_gettime(time.CLOCK_MONOTONIC)-began)
        if remaining>0:
            lkh=polish(d,route,seconds=min(1.,remaining))
            # Recheck the OUTER deadline, including call/return overhead.
            if time.clock_gettime(time.CLOCK_MONOTONIC)-began<=10. and lkh['cost']<result['cost']:
                result['cost']=lkh['cost'];result['route']=lkh['route']
            result['lkh_used']=True;result['lkh']=lkh
    result['wall_seconds']=time.clock_gettime(time.CLOCK_MONOTONIC)-began
    result['overrun']=max(0.,result['wall_seconds']-10.)
    if result['feasible']:
        d=read_instance(ROOT/'data'/name);mapping={int(v):i for i,v in enumerate(d['visits'])}
        assert validate(d,[mapping[v] for v in result['route']])==result['cost']
    target.write_text(json.dumps(result,indent=2)+'\n');return result


if __name__=='__main__':
    OUT.mkdir(exist_ok=True)
    cases=json.loads((ROOT/'experiments/original_mdm/hard_training_baseline/summary.json').read_text())['rows']
    assert len(cases)==9
    for case in cases:
        for seed in (42,43,44):
            order=['mdm','hybrid'];random.Random(case['instance']+str(seed)).shuffle(order)
            pair={mode:run(case['instance'],seed,mode) for mode in order}
            print(json.dumps(dict(instance=case['instance'],seed=seed,
                mdm_cost=pair['mdm']['cost'],hybrid_cost=pair['hybrid']['cost'],
                lkh_used=pair['hybrid']['lkh_used'])),flush=True)
    (OUT/'completion.json').write_text(json.dumps(dict(cases=9,seeds=[42,43,44],arms=['mdm','hybrid']))+'\n')
