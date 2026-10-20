"""Conditional pair orchestration, intentionally without a launch entrypoint.

A future collection owner MUST hold paired70_train.study_lock on the same
paired70_train.queue.lock inode and verify prerequisites before calling this.
No active collection imports this module. Independent batch auditing and its
source-bound fixtures remain prerequisites for a real launch.
"""
import hashlib
import json
import random
import sys
import time
from pathlib import Path

# The frozen awake guard imports its original worker by local module name.
# Support importing this adapter from the repository root as well as scripts.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from experiments.lkh.awake300_train import ClockGuard, sleep_events_since
from experiments.lkh.cached30_phase import run_phase
from experiments.lkh.cached30_phase_audit import audit_phase
from experiments.lkh.raw_route_audit import load


def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def execute_pairs(starts, out, protocol_hash, phase_sources, budget=30.):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    assert len({(s['instance'], s['source_seed']) for s in starts}) == len(starts)
    # Scan the entire schedule before any new optimizer call, including a
    # partial pair later in the shuffled order than untouched pairs.
    for suffix in ('*.attempt.json', '*.failed.json', '*.tmp'):
        assert not list(out.glob(suffix)), 'Unreviewed interrupted artifact; no optimizer launched'
    for start in starts:
        prefix = f"{start['instance']}.seed{start['source_seed']}"
        if not (out/(prefix+'.pair.json')).exists():
            assert not list(out.glob(prefix+'.*')) and not (out/('environment.'+prefix+'.json')).exists(), \
                'Unreviewed partial pair anywhere in schedule; no optimizer launched'
    for start in starts:
        name, seed = start['instance'], start['source_seed']
        prefix = f'{name}.seed{seed}'
        target = out/(prefix+'.pair.json')
        environment = out/('environment.'+prefix+'.json')
        marker = out/(prefix+'.attempt.json')
        order = ['uncached', 'cached']
        random.Random(prefix+'cached30').shuffle(order)
        data = load(name)
        assert data['input_hash'] == start['input_sha256']
        if target.exists():
            record = json.loads(target.read_text())
            assert not marker.exists(), 'Review pair committed before its attempt marker was removed'
            assert record['status'] == 'PAIR_COMPLETED' and record['start'] == start
            assert record['protocol_sha256'] == protocol_hash and record['arm_order'] == order
            assert record['budget_seconds'] == budget
            assert record['environment_file'] == environment.name
            assert digest(environment) == record['environment_sha256']
            assert json.loads(environment.read_text())['status'] == 'PASS'
            assert set(record['phases']) == {'uncached', 'cached'}
            for arm in order:
                path = out/(prefix+'.'+arm+'.json')
                descriptor = record['phases'][arm]
                assert descriptor['file'] == path.name and descriptor['sha256'] == digest(path)
                phase = json.loads(path.read_text())
                audit_phase(phase, data, start, arm, budget, phase_sources)
                assert descriptor['cost'] == phase['cost']
            continue
        artifacts = list(out.glob(prefix+'.*'))
        if environment.exists(): artifacts.append(environment)
        assert not artifacts, 'Preserve and review partial pair before any new optimizer call'
        with marker.open('x') as handle:
            json.dump(dict(start=start, protocol_sha256=protocol_hash), handle)
        phases = {}
        with ClockGuard() as guard:
            for arm in order:
                path = out/(prefix+'.'+arm+'.json')
                phase = run_phase(start, seed+49, arm, budget, path)
                audit_phase(phase, data, start, arm, budget, phase_sources)
                phases[arm] = dict(file=path.name, sha256=digest(path), cost=phase['cost'])
        environmental = dict(status='PASS', started_unix=guard.started, finished_unix=time.time(),
                             max_monitor_gap_seconds=guard.max_gap,
                             max_clock_difference_seconds=guard.clock_difference,
                             sleep_events=sleep_events_since(guard.started))
        if guard.max_gap > 5 or guard.clock_difference > 3 or environmental['sleep_events']:
            environmental['status'] = 'ENVIRONMENT_FAILURE'
        environment.write_text(json.dumps(environmental, indent=2)+'\n')
        assert environmental['status'] == 'PASS', 'Environmental failure preserved; no next pair'
        record = dict(status='PAIR_COMPLETED', start=start, protocol_sha256=protocol_hash,
                      arm_order=order, phases=phases, budget_seconds=budget,
                      environment_file=environment.name, environment_sha256=digest(environment))
        temp = target.with_suffix('.tmp')
        temp.write_text(json.dumps(record, indent=2)+'\n'); temp.replace(target)
        marker.unlink()
        print(json.dumps(dict(instance=name, source_seed=seed,
                              costs={arm: phases[arm]['cost'] for arm in order})), flush=True)
