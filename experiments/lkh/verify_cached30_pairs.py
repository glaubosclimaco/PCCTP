"""Stub pair orchestration checks; no benchmark read or optimizer launch."""
import contextlib
import hashlib
import io
import json
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from experiments.lkh import cached30_pairs as module


def main():
    state = dict(calls=[], audited=[], fail=False, sleep=False)
    sources = dict(synthetic='hash')
    starts = [dict(instance='synthetic'+str(i), source_seed=42+i,
                   cost=5, route=[0, 1, 2], input_sha256='synthetic') for i in range(3)]
    original = module.ClockGuard, module.sleep_events_since, module.run_phase, module.audit_phase, module.load

    class Guard:
        def __enter__(self):
            self.started = 100.; self.max_gap = 0.; self.clock_difference = 0.
            return self
        def __exit__(self, *args): return False

    def phase(start, seed, arm, budget, path):
        state['calls'].append((start['instance'], seed, arm))
        if state['fail']: raise RuntimeError('synthetic interruption')
        result = dict(start=start, repair_seed=seed, arm=arm, budget=budget,
                      cost=5 if arm == 'uncached' else 4, source_sha256=sources)
        path.write_text(json.dumps(result)+'\n')
        return result

    def audit(record, data, start, arm, budget, expected):
        assert record['start'] == start and record['arm'] == arm
        assert record['budget'] == budget and record['source_sha256'] == expected
        state['audited'].append((start['instance'], arm))

    checks = []
    def rejected(study, root, protocol='protocol'):
        before = len(state['calls'])
        try: module.execute_pairs(study, root, protocol, sources)
        except AssertionError: pass
        else: raise AssertionError('Expected rejection')
        assert len(state['calls']) == before

    try:
        module.ClockGuard = Guard
        module.sleep_events_since = lambda _: ['synthetic sleep'] if state['sleep'] else []
        module.run_phase, module.audit_phase = phase, audit
        module.load = lambda _: dict(input_hash='synthetic')
        with tempfile.TemporaryDirectory(prefix='pcctp-cache-pairs-') as folder, contextlib.redirect_stdout(io.StringIO()):
            root = Path(folder)/'normal'
            module.execute_pairs(starts, root, 'protocol', sources)
            assert len(state['calls']) == len(state['audited']) == 6
            for start in starts:
                calls = [c for c in state['calls'] if c[0] == start['instance']]
                assert {c[2] for c in calls} == {'uncached', 'cached'}
                assert all(c[1] == start['source_seed']+49 for c in calls)
            checks.append('three pairs: both arms, common start and prescribed repair seeds')
            before = len(state['calls'])
            module.execute_pairs(starts, root, 'protocol', sources)
            assert len(state['calls']) == before and len(state['audited']) == 12
            checks.append('completed pairs reaudited and resumed with zero new calls')
            rejected(starts, root, 'changed protocol')
            checks.append('changed protocol blocks resume before optimizer')
            path = next(root.glob('*.uncached.json')); text = path.read_text()
            path.write_text(text+' ')
            rejected(starts, root)
            path.write_text(text)
            checks.append('changed phase bytes block resume before optimizer')
            env = next(root.glob('environment.*.json')); text = env.read_text()
            env.write_text(text+' ')
            rejected(starts, root)
            env.write_text(text)
            checks.append('changed environment bytes block resume before optimizer')
            phase_path = root/'synthetic0.seed42.uncached.attempt.json'
            partial = Path(folder)/'partial'; partial.mkdir()
            (partial/phase_path.name).write_text('{}')
            rejected(starts[:1], partial)
            checks.append('orphan attempt marker blocks partial pair retry')
            later = Path(folder)/'later-partial'; later.mkdir()
            (later/'synthetic2.seed44.cached.json').write_text('{}')
            rejected(starts, later)
            checks.append('partial pair later in schedule blocks even earlier untouched pairs before any call')

            state['fail'] = True; failed = Path(folder)/'interrupted'
            try: module.execute_pairs(starts[:1], failed, 'protocol', sources)
            except RuntimeError: pass
            else: raise AssertionError('Interruption swallowed')
            assert (failed/'synthetic0.seed42.attempt.json').exists()
            state['fail'] = False
            rejected(starts[:1], failed)
            checks.append('interrupted pair preserves exclusive marker and blocks retry')

            state['sleep'] = True; badenv = Path(folder)/'sleep'
            before = len(state['calls'])
            try: module.execute_pairs(starts, badenv, 'protocol', sources)
            except AssertionError: pass
            else: raise AssertionError('Compromised environment accepted')
            assert len(state['calls']) == before+2
            assert json.loads(next(badenv.glob('environment.*.json')).read_text())['status'] == 'ENVIRONMENT_FAILURE'
            assert not list(badenv.glob('*.pair.json'))
            state['sleep'] = False
            rejected(starts, badenv)
            checks.append('sleep invalidates pair, preserves phases, stops before next pair and blocks retry')
    finally:
        module.ClockGuard, module.sleep_events_since, module.run_phase, module.audit_phase, module.load = original
    files = ['experiments/lkh/'+name for name in ('cached30_pairs.py', 'verify_cached30_pairs.py',
        'cached30_phase.py', 'cached30_phase_audit.py', 'cached30_starts.py',
        'raw_route_audit.py', 'awake300_train.py', 'paired70_train.py')]
    report = dict(status='PASS', scope='Synthetic pair orchestration only; scalar phase audit verified '
                  'separately; zero benchmark reads and native optimizer launches.', checks=checks,
                  source_sha256={p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in files})
    Path(__file__).with_suffix('.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(status='PASS', checks=len(checks), native_optimizer_calls=0)))


if __name__ == '__main__': main()
