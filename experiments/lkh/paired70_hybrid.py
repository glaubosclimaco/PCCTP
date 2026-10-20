"""Candidate with a reference-defined duration; no reference solution is an input."""
import json,math,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from experiments.lkh.feasible_fallback import construct
from experiments.lkh.ruin_repair_absolute import run as search
from experiments.original_mdm.timed_experiment import run_one
from experiments.ml_lns.pilot import read_instance,validate


def select_events(cost,route,events,parent_started,budget):
    assert math.isfinite(parent_started) and math.isfinite(budget) and budget>0
    origin=None
    for event in events:
        timestamp=event['native_timestamp'];value=event['cost']
        assert math.isfinite(timestamp) and timestamp>=parent_started
        assert math.isfinite(value)
        if timestamp<=parent_started+budget and (cost is None or value<cost):
            cost,route,origin=value,event['route'][:],event
    return cost,route,origin


def run(name,seed,mode,budget,out):
    assert mode in ('hybrid10','hybrid60')
    assert math.isfinite(budget) and budget>0
    out=Path(out);target=out/f'{name}.seed{seed}.{mode}.json'
    if target.exists():
        existing=json.loads(target.read_text())
        assert existing['budget_seconds']==budget and existing['mode']==mode
        return existing
    out.mkdir(exist_ok=True)
    began=time.clock_gettime(time.CLOCK_MONOTONIC)
    d=read_instance(ROOT/'data'/name)
    f=construct(d)
    fallback=dict(cost=f['cost'],route=[int(d['visits'][v]) for v in f['route']],
                  seconds=time.clock_gettime(time.CLOCK_MONOTONIC)-began)
    cost,route,origin=(fallback['cost'],fallback['route'],'fallback') if fallback['seconds']<=budget else (None,[],None)
    initial_budget=min(float(mode[6:]),budget)
    native=run_one(name,seed,'original',initial_budget,
        out/f'{name}.seed{seed}.{mode}.native.json',True,
        parent_started=began,hard_timeout_seconds=budget+30.)
    native_events=[dict(native_timestamp=began+e['seconds'],cost=e['cost'],route=e['route'])
                   for e in native['anytime']]
    assert all(e['native_timestamp']<=began+initial_budget+1e-9 for e in native_events)
    cost,route,accepted=select_events(cost,route,native_events,began,budget)
    if accepted is not None:origin='native'
    pre_cost,pre_route=cost,route[:]
    offset=time.clock_gettime(time.CLOCK_MONOTONIC)-began;post=None
    if route and offset<budget:
        post=search(name,'random',budget-offset,initial_override=route,seed=seed)
        cost,route,accepted=select_events(cost,route,post['trace'],began,budget)
        if accepted is not None:origin='post_event'
    elapsed=time.clock_gettime(time.CLOCK_MONOTONIC)-began
    if route:
        mapping={int(v):i for i,v in enumerate(d['visits'])}
        assert validate(d,[mapping[v] for v in route])==cost
    result=dict(instance=name,seed=seed,mode=mode,budget_seconds=budget,
        initialization_budget_seconds=initial_budget,native_started=began,
        fallback=fallback,native=native,pre_cost=pre_cost,pre_route=pre_route,
        post_offset=offset,post=post,cost=cost,route=route,feasible=bool(route),
        origin=origin,wall_seconds=elapsed,
        scope='fresh candidate; only reference duration supplied; reference route/cost unavailable')
    temp=target.with_suffix('.tmp');temp.write_text(json.dumps(result,indent=2)+'\n');temp.replace(target)
    return result
