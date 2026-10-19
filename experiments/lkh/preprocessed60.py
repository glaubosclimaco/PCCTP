"""Frozen offline nearest-edge expansion with a matched random-edge control."""
import json, time, random, hashlib, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'experiments/ml_lns'))
from pilot import read_instance, validate
from highs_warm_start import solve
from edge_pool_study import edges
sys.path.insert(0,str(ROOT/'experiments/preprocessing'))
from structural import reduce_constraints

if __name__ == '__main__':
    assert (HERE / 'edge_pool_study/completion.json').exists()
    assert json.loads((HERE / 'edge_pool_audit.json').read_text())['batch_complete']
    out = HERE / 'preprocessed60'; out.mkdir(exist_ok=True)
    sources = sorted((HERE / 'edge_pool_study').glob('g_*.json'))
    assert len(sources) == 9
    for case, source in enumerate(sources):
        for arm in ('full',):
            target = out / (source.stem + '.' + arm + '.json')
            if target.exists(): continue
            began = time.perf_counter(); native = time.clock_gettime(time.CLOCK_MONOTONIC)
            cpu = time.process_time(); old = json.loads(source.read_text())
            d = read_instance(ROOT / 'data' / old['instance'])
            prep_start=time.perf_counter()
            reduced=reduce_constraints(range(len(d['visits'])),d['mandatory'],d['prizes'],d['quota'],d['neighborhoods'])
            d=dict(d,mandatory=set(reduced['mandatory']),neighborhoods=[set(s) for s in reduced['neighborhoods']])
            preprocessing_seconds=time.perf_counter()-prep_start
            mapping = {int(v): i for i, v in enumerate(d['visits'])}
            routes = [[mapping[v] for v in r] for r in old['input_routes']]
            initial = min(routes, key=lambda r: validate(d, r)); cost = validate(d, initial)
            pool = set().union(*(edges(r) for r in routes)); n = len(mapping)
            # The matrix below is indexed by internal visit indices.
            near = set()
            for a in range(n):
                ranked = sorted((b for b in range(n) if b != a), key=lambda b: (d['costs'][a, b], b))
                near.update(tuple(sorted((a, b))) for b in ranked[:5])
            count = len(near - pool)
            if arm == 'full': allowed = {(a,b) for a in range(n) for b in range(a+1,n)}; count=len(allowed)-len(pool)
            elif arm == 'nearest': allowed = pool | near
            else:
                absent = [(a, b) for a in range(n) for b in range(a + 1, n) if (a, b) not in pool]
                allowed = pool | set(random.Random(42).sample(absent, count))
            assert len(allowed) == len(pool) + count and edges(initial) <= allowed
            ready = time.clock_gettime(time.CLOCK_MONOTONIC) - native
            try: result = solve(d, initial, True, seconds=60., began=began, allowed_edges=allowed)
            except RuntimeError as error:
                if str(error) != 'Setup exhausted budget': raise
                result = dict(route=[int(d['visits'][v]) for v in initial], cost=cost,
                              status=str(error), proven_optimal=False, proof_scope='restricted edge pool')
            returned = time.clock_gettime(time.CLOCK_MONOTONIC) - native
            raw = result; discarded = returned > 60 and result['cost'] < cost
            if discarded:
                result = dict(route=[int(d['visits'][v]) for v in initial], cost=cost,
                              status='Late improvement discarded', proven_optimal=False, proof_scope='restricted edge pool')
            r = [mapping[v] for v in result['route']]
            assert validate(d, r) == result['cost'] <= cost and edges(r) <= allowed
            row = dict(instance=old['instance'], arm=arm, initial_cost=cost, result=result,
                       raw_result=raw, late_improvement_discarded=discarded,
                       initial_ready_native_seconds=ready, native_return_seconds=returned,
                       cpu_seconds=time.process_time()-cpu, added_edges=count,preprocessing=reduced,preprocessing_seconds=preprocessing_seconds,
                       reduction_source_sha256=hashlib.sha256((ROOT/'experiments/preprocessing/structural.py').read_bytes()).hexdigest(),
                       pool_edges=[[int(d['visits'][a]), int(d['visits'][b])] for a,b in sorted(allowed)],
                       source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                       scope='offline full graph; input generation excluded; complete visit graph')
            target.write_text(json.dumps(row, indent=2)+'\n')
            print(json.dumps(dict(instance=old['instance'], arm=arm, cost=result['cost'])), flush=True)
    (out / 'completion.json').write_text(json.dumps(dict(cases=9, arms=1))+'\n')
