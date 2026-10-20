"""Independent file/environment/phase reconstruction for the conditional study."""
import hashlib
import json
import math
import random
from fractions import Fraction
from pathlib import Path

from experiments.lkh.cached30_phase_audit import audit_phase
from experiments.lkh.raw_route_audit import load


def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def audit_batch(out, starts, protocol_hash, phase_sources, data_loader=load, require_complete=True):
    out = Path(out)
    keys = {(s['instance'], s['source_seed']) for s in starts}
    geometries = {s['geometry'] for s in starts}
    assert len(starts) == len(keys) == 27 and len(geometries) == 9
    for geometry in geometries:
        group = [s for s in starts if s['geometry'] == geometry]
        assert len(group) == 3 and {s['source_seed'] for s in group} == {42, 43, 44}
        assert len({s['instance'] for s in group}) == 1
    expected = {f"{s['instance']}.seed{s['source_seed']}.pair.json" for s in starts}
    actual = {p.name for p in out.glob('*.pair.json')}
    assert actual <= expected, 'Unexpected or duplicate pair identity'
    if require_complete: assert actual == expected, 'All 27 pairs required'
    rows, hashes, routes, totals = [], {}, 0, {g: {'uncached': [], 'cached': []} for g in geometries}
    for start in starts:
        prefix = f"{start['instance']}.seed{start['source_seed']}"
        path = out/(prefix+'.pair.json')
        if not path.exists(): continue
        assert not any(out.glob(prefix+'*.attempt.json'))
        assert not any(out.glob(prefix+'*.failed.json'))
        assert not any(out.glob(prefix+'*.tmp'))
        pair = json.loads(path.read_text())
        assert pair['status'] == 'PAIR_COMPLETED' and pair['start'] == start
        assert pair['protocol_sha256'] == protocol_hash and pair['budget_seconds'] == 30.
        order = ['uncached', 'cached']; random.Random(prefix+'cached30').shuffle(order)
        assert pair['arm_order'] == order and set(pair['phases']) == set(order)
        envpath = out/('environment.'+prefix+'.json')
        assert pair['environment_file'] == envpath.name and pair['environment_sha256'] == digest(envpath)
        env = json.loads(envpath.read_text())
        assert env['status'] == 'PASS' and env['sleep_events'] == []
        for field in ('started_unix', 'finished_unix', 'max_monitor_gap_seconds', 'max_clock_difference_seconds'):
            assert math.isfinite(env[field])
        assert env['finished_unix'] >= env['started_unix']
        assert 0 <= env['max_monitor_gap_seconds'] <= 5
        assert 0 <= env['max_clock_difference_seconds'] <= 3
        data = data_loader(start['instance'])
        assert data['input_hash'] == start['input_sha256']
        phases = {}; previous_end = None
        for arm in order:
            phasepath = out/(prefix+'.'+arm+'.json')
            descriptor = pair['phases'][arm]
            assert descriptor['file'] == phasepath.name and descriptor['sha256'] == digest(phasepath)
            phase = json.loads(phasepath.read_text())
            verified = audit_phase(phase, data, start, arm, 30., phase_sources)
            routes += verified['routes_checked']
            assert descriptor['cost'] == phase['cost']
            if previous_end is not None: assert phase['parent_started'] >= previous_end-1e-5
            previous_end = phase['parent_started']+phase['wall_seconds']
            phases[arm] = phase
            totals[start['geometry']][arm].append(Fraction(str(phase['cost'])))
            hashes[phasepath.name] = digest(phasepath)
        assert sum(p['wall_seconds'] for p in phases.values()) <= env['finished_unix']-env['started_unix']+1e-3
        hashes[path.name], hashes[envpath.name] = digest(path), digest(envpath)
        rows.append(dict(instance=start['instance'], geometry=start['geometry'], source_seed=start['source_seed'],
                         costs={a: p['cost'] for a, p in phases.items()},
                         counters={a: {k: p['raw_phase'][k] for k in ('attempts', 'feasible_repairs',
                            'lkh_calls', 'lkh_completed', 'lkh_timeouts', 'lkh_failed')} for a, p in phases.items()},
                         cache_retained_python_bytes=phases['cached']['cache_retained_python_bytes'],
                         cache_build_events=phases['cached']['cache_build_events']))
    complete = actual == expected
    if complete:
        assert not list(out.glob('*.attempt.json')) and not list(out.glob('*.failed.json'))
        assert not list(out.glob('*.tmp'))
        allowed = expected | {'environment.'+p.replace('.pair.json', '.json') for p in expected}
        allowed |= {p.replace('.pair.json', '.'+arm+'.json') for p in expected for arm in ('cached', 'uncached')}
        assert {p.name for p in out.glob('*.json')} <= allowed | {'protocol.json', 'completion.json'}
    effects = []
    for geometry, arms in sorted(totals.items()):
        if len(arms['cached']) != 3: continue
        u, c = sum(arms['uncached'])/3, sum(arms['cached'])/3
        assert u > 0
        gain = 100*(u-c)/u
        effects.append(dict(geometry=geometry, mean_uncached=float(u), mean_cached=float(c),
                            gain_percent=float(gain), exact_gain_fraction=str(gain)))
    exact = [Fraction(e['exact_gain_fraction']) for e in effects]
    mean = sum(exact)/len(exact) if exact else None
    wins, losses = sum(g > 0 for g in exact), sum(g < 0 for g in exact)
    gate = complete and mean >= Fraction(1, 5) and wins >= 6 and losses <= 1
    return dict(status='PASS' if complete else 'PARTIAL', batch_complete=complete,
                protocol_sha256=protocol_hash,
                scheduled_pairs=27, scheduled_phases=54, audited_pairs=len(rows), audited_phases=2*len(rows),
                routes_checked=routes, rows=rows, geometries=effects,
                mean_gain_percent=float(mean) if mean is not None else None,
                exact_mean_gain_fraction=str(mean) if mean is not None else None,
                wins=wins, losses=losses, ties=len(exact)-wins-losses, gate_passed=bool(gate),
                artifact_sha256=hashes, scope='Conditional reused TRAIN starts; generation excluded. '
                'Gate is exploratory; not a complete-search speedup or independent confirmation.')
