"""Synthetic phase/auditor contract verification; no benchmark or optimizer call."""
import copy
import hashlib
import json
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from experiments.lkh import cached30_phase as module
from experiments.lkh.cached30_phase_audit import audit_phase


def main():
    clock = {'now': 100., 'fail': False, 'runs': 0, 'mode': 'normal'}
    data = dict(n=5, types=[1, 0, 0, 0, 0], prizes=[1]*5, quota=1,
                radius=0, matrix=[int(i != j) for i in range(5) for j in range(5)])
    start = dict(instance='synthetic', source_seed=42, cost=5, route=[0, 1, 2, 3, 4])
    original = module.time, module.core.run, module.core.polish, module.ORIGINAL_POLISH, module.CachedPolisher

    def polish(*args, **kwargs):
        clock['now'] += .03
        return dict(matrix_cache_created=False)

    class FakeCached:
        def __init__(self): self.cache = None
        def __call__(self, *args, **kwargs):
            created = self.cache is None
            if created: self.cache = SimpleNamespace(rows=(('0', '1'), ('1', '0')))
            clock['now'] += .03
            return dict(matrix_cache_created=created, matrix_cache_build_seconds=.01 if created else 0)

    def run(name, arm, budget, initial_override, seed):
        clock['runs'] += 1
        assert initial_override == start['route'] and seed == 91 and budget == 2.
        assert module.core.polish is not polish
        clock['now'] += .02
        native = clock['now']
        if clock['mode'] == 'empty':
            return dict(instance=name, arm=arm, initial_cost=5, cost=5, route=start['route'][:],
                        native_started=native, native_elapsed=0, trace=[], source_sha256=None,
                        attempts=0, feasible_repairs=0, lkh_calls=0, lkh_completed=0,
                        lkh_timeouts=0, lkh_failed=0)
        module.core.polish(None, initial_override, seconds=.2)
        if clock['fail']: raise RuntimeError('synthetic interrupted phase')
        trace = [dict(native_timestamp=100.5, native_seconds=100.5-native,
                      cost=4, route=[0, 1, 2, 3]),
                 dict(native_timestamp=102.01, native_seconds=102.01-native,
                      cost=3, route=[0, 1, 2])]
        if clock['mode'] == 'boundary':
            trace[-1].update(native_timestamp=102., native_seconds=102.-native)
        clock['now'] = 102.1
        return dict(instance=name, arm=arm, initial_cost=5, cost=3, route=[0, 1, 2],
                    native_started=native, native_elapsed=clock['now']-native, trace=trace,
                    source_sha256=None, attempts=1, feasible_repairs=1,
                    lkh_calls=2 if clock['mode'] == 'bad_count' else 1,
                    lkh_completed=1, lkh_timeouts=0, lkh_failed=0)

    checks = []
    try:
        module.time = SimpleNamespace(CLOCK_MONOTONIC=1, clock_gettime=lambda _: clock['now'])
        module.core.run = run
        module.core.polish = module.ORIGINAL_POLISH = polish
        module.CachedPolisher = FakeCached
        records = {}
        with tempfile.TemporaryDirectory(prefix='pcctp-cache-phase-') as folder:
            for arm in ('uncached', 'cached'):
                clock['now'] = 100.
                path = Path(folder)/(arm+'.json')
                record = module.run_phase(start, 91, arm, 2., path)
                assert module.core.polish is polish
                assert json.loads(path.read_text()) == record
                assert not path.with_suffix('.attempt.json').exists()
                assert record['cost'] == 4 and record['raw_phase']['cost'] == 3
                assert record['parent_started']+2 < record['raw_phase']['trace'][-1]['native_timestamp']
                assert record['raw_phase']['trace'][-1]['native_timestamp'] < record['raw_phase']['native_started']+2
                result = audit_phase(record, data, start, arm, 2., module.source_hashes())
                assert result['routes_checked'] == 5 and result['eligible_phase_events'] == 1
                records[arm] = record
                checks.append(arm+': audited parent deadline excludes late improvement retained in raw phase')
            assert records['uncached']['cache_retained_python_bytes'] == 0
            assert records['cached']['cache_retained_python_bytes'] > 0
            assert len(records['cached']['cache_build_events']) == 1
            checks.append('single lazy build measured inside phase; approximate owned memory recorded')

            for mode in ('empty', 'boundary'):
                clock.update(now=100., mode=mode)
                record = module.run_phase(start, 91, 'cached', 2., Path(folder)/(mode+'.json'))
                audit_phase(record, data, start, 'cached', 2., module.source_hashes())
                assert record['cost'] == (5 if mode == 'empty' else 3)
                checks.append(mode+': no-call fallback or event exactly at outer deadline remains eligible')
            clock['mode'] = 'normal'

            mutations = {
                'wrong endpoint': lambda r: r.update(cost=3, route=[0, 1, 2]),
                'wrong start': lambda r: r['start'].update(source_seed=43),
                'wrong source': lambda r: r['source_sha256'].update(fabricated='hash'),
                'wrong budget': lambda r: r.update(budget_seconds=3.),
                'nonfinite wall': lambda r: r.update(wall_seconds=float('nan')),
                'infeasible route': lambda r: r['raw_phase']['trace'][0].update(route=[1, 2, 3, 4]),
                'wrong scalar cost': lambda r: r['raw_phase']['trace'][0].update(cost=3),
                'bad event clock': lambda r: r['raw_phase']['trace'][0].update(native_seconds=0),
                'counter partition': lambda r: r['raw_phase'].update(lkh_completed=0),
                'missing build': lambda r: r.update(cache_build_events=[]),
                'missing footprint': lambda r: r.update(cache_retained_python_bytes=0),
                'unaccounted build time': lambda r: r['cache_build_events'][0].update(build_seconds=1.),
                'nonfinite build clock': lambda r: r['cache_build_events'][0].update(call_started=float('nan')),
                'changed raw final': lambda r: r['raw_phase'].update(cost=4, route=[0, 1, 2, 3]),
            }
            for label, mutate in mutations.items():
                record = copy.deepcopy(records['cached']); mutate(record)
                try: audit_phase(record, data, start, 'cached', 2., module.source_hashes())
                except AssertionError: pass
                else: raise AssertionError('Accepted corruption: '+label)
                checks.append('auditor rejects '+label)

            clock.update(now=100., fail=True)
            target = Path(folder)/'interrupted.json'
            try: module.run_phase(start, 91, 'cached', 2., target)
            except RuntimeError: pass
            else: raise AssertionError('Synthetic failure swallowed')
            assert module.core.polish is polish and not target.exists()
            failure = target.with_suffix('.failed.json')
            assert json.loads(failure.read_text())['error_type'] == 'RuntimeError'
            previous = clock['runs']
            try: module.run_phase(start, 91, 'cached', 2., target)
            except AssertionError: pass
            else: raise AssertionError('Failed phase silently retried')
            assert clock['runs'] == previous
            checks.append('exception restores hook, records failure and blocks silent retry')

            clock.update(now=100., fail=False, mode='bad_count')
            target = Path(folder)/'postsearch-error.json'
            try: module.run_phase(start, 91, 'cached', 2., target)
            except AssertionError: pass
            else: raise AssertionError('Bad call accounting accepted')
            assert module.core.polish is polish and target.with_suffix('.attempt.json').exists()
            previous = clock['runs']
            try: module.run_phase(start, 91, 'cached', 2., target)
            except AssertionError: pass
            else: raise AssertionError('Post-search failure silently retried')
            assert clock['runs'] == previous
            checks.append('post-search assertion leaves attempt marker and blocks silent retry')
    finally:
        module.time, module.core.run, module.core.polish, module.ORIGINAL_POLISH, module.CachedPolisher = original
    sources = list(module.SOURCES)+['experiments/lkh/cached30_phase_audit.py',
                                  'experiments/lkh/raw_route_audit.py',
                                  'experiments/lkh/verify_cached30_phase.py']
    report = dict(status='PASS', scope='Synthetic adapter and independent scalar auditor contracts; '
                  'no native optimizer or benchmark inputs; not evidence of cache search benefit.',
                  checks=checks, source_sha256={p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources})
    (Path(__file__).with_suffix('.json')).write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(status='PASS', checks=len(checks), native_optimizer_calls=0)))


if __name__ == '__main__': main()
