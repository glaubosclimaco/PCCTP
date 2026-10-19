"""Offline route recombination; refuses to overlap the fixed-work batch."""
import sys,json,time,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'experiments/ml_lns'))
from pilot import read_instance,validate
from highs_warm_start import solve
HERE=Path(__file__).resolve().parent

def edges(route):return {tuple(sorted((route[i-1],v))) for i,v in enumerate(route)}
if __name__=='__main__':
    assert (HERE/'fixed_work/completion.json').exists(),'Wait for all nine fixed-work cases; do not run concurrently'
    verification=json.loads((HERE/'edge_pool_verification.json').read_text());assert verification['status']=='PASS'
    cases=json.loads((ROOT/'experiments/original_mdm/hard_training_baseline/summary.json').read_text())['rows'];split=json.loads((ROOT/'experiments/ml_lns/decoder_initialized_states/protocol.json').read_text());assert all(split['roles'][c['geometry']]=='train' for c in cases)
    out=HERE/'edge_pool_study';out.mkdir(exist_ok=True)
    for c in cases:
        name=c['instance'];target=out/f'{name}.json'
        if target.exists():continue
        began=time.perf_counter();native_began=time.clock_gettime(time.CLOCK_MONOTONIC);cpu=time.process_time()
        d=read_instance(ROOT/'data'/name);mapping={int(v):i for i,v in enumerate(d['visits'])};routes=[];hashes={};missing=[]
        def offer(route):
            route=[mapping[v] for v in route];validate(d,route)
            if route not in routes:routes.append(route)
        path=HERE/'fixed_work'/f'{name}.json';fixed=json.loads(path.read_text());hashes[str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest()
        if fixed['status']!='NO_VALID_EXPORT':offer(fixed['native_route']);offer(fixed['lkh']['route'])
        else:missing.append('fixed_work valid export')
        for mode in ('mdm','fixed_lkh'):
            for seed in (42,43,44):
                path=HERE/'support_common_fallback'/f'{name}.seed{seed}.{mode}.json';r=json.loads(path.read_text());hashes[str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest()
                assert r['feasible'];offer(r['route'])
        assert routes
        pool=set().union(*(edges(route) for route in routes));initial=min(routes,key=lambda route:validate(d,route));initial_cost=validate(d,initial)
        initial_ready_native_seconds=time.clock_gettime(time.CLOCK_MONOTONIC)-native_began
        try:r=solve(d,initial,True,seconds=60.,began=began,allowed_edges=pool)
        except RuntimeError as error:
            if str(error)!='Setup exhausted budget':raise
            r=dict(cost=initial_cost,route=[int(d['visits'][v]) for v in initial],status='Setup exhausted budget',proven_optimal=False,proof_scope='restricted edge pool')
        raw_result=r
        native_return_seconds=time.clock_gettime(time.CLOCK_MONOTONIC)-native_began
        late_improvement_discarded=native_return_seconds>60 and r['cost']<initial_cost
        if late_improvement_discarded:
            r=dict(cost=initial_cost,route=[int(d['visits'][v]) for v in initial],status='Outer native-clock deadline exceeded; improvement discarded',proven_optimal=False,proof_scope='restricted edge pool')
        returned=[mapping[v] for v in r['route']];assert validate(d,returned)==r['cost']<=initial_cost and edges(returned)<=pool
        row=dict(instance=name,scope='offline edge pool; input generation time excluded; bounds local to pool',initial_cost=initial_cost,input_routes=[[int(d['visits'][v]) for v in route] for route in routes],pool_edges=[[int(d['visits'][a]),int(d['visits'][b])] for a,b in sorted(pool)],result=r,raw_result=raw_result,initial_ready_native_seconds=initial_ready_native_seconds,feasible_before_deadline=initial_ready_native_seconds<=60,late_improvement_discarded=late_improvement_discarded,input_hashes=hashes,missing_sources=missing,python_elapsed=time.perf_counter()-began,native_elapsed=time.clock_gettime(time.CLOCK_MONOTONIC)-native_began,cpu_seconds=time.process_time()-cpu)
        target.write_text(json.dumps(row,indent=2)+'\n');print(json.dumps(dict(instance=name,initial_cost=initial_cost,cost=r['cost'],pool_edges=len(pool))),flush=True)
    (out/'completion.json').write_text(json.dumps(dict(cases=9))+'\n')
