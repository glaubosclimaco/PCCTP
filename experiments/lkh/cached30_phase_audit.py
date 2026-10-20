"""Independent scalar reconstruction of conditional phase eligibility."""
import math
from experiments.lkh.raw_route_audit import check


def audit_phase(record, data, start, arm, budget, expected_sources):
    assert record['start'] == start and record['instance'] == start['instance']
    assert start['source_seed'] in (42, 43, 44)
    assert record['repair_seed'] == start['source_seed']+49
    assert record['arm'] == arm and arm in ('uncached', 'cached')
    assert record['budget_seconds'] == budget and math.isfinite(budget) and budget > 0
    assert record['source_sha256'] == expected_sources
    began, elapsed = record['parent_started'], record['wall_seconds']
    assert math.isfinite(began) and math.isfinite(elapsed) and elapsed >= 0
    check(data, start['route'], start['cost'])
    raw = record['raw_phase']
    assert raw['instance'] == start['instance'] and raw['arm'] == 'random'
    assert raw['source_sha256'] is None and raw['initial_cost'] == start['cost']
    native, duration = raw['native_started'], raw['native_elapsed']
    assert math.isfinite(native) and math.isfinite(duration) and duration >= 0
    assert began <= native <= began+elapsed+1e-5
    assert native+duration <= began+elapsed+1e-5
    for name in ('attempts', 'feasible_repairs', 'lkh_calls', 'lkh_completed', 'lkh_timeouts', 'lkh_failed'):
        assert type(raw[name]) is int and raw[name] >= 0
    assert raw['attempts'] >= raw['feasible_repairs'] >= raw['lkh_calls']
    assert raw['lkh_calls'] == raw['lkh_completed']+raw['lkh_timeouts']+raw['lkh_failed']
    cost, route, origin = start['cost'], start['route'], 'common_start'
    previous_cost, previous_time = start['cost'], native
    raw_cost, raw_route = start['cost'], start['route']
    checked = 1
    for event in raw['trace']:
        check(data, event['route'], event['cost']); checked += 1
        timestamp = event['native_timestamp']
        assert math.isfinite(timestamp) and previous_time <= timestamp <= native+duration+1e-5
        assert math.isfinite(event['native_seconds']) and event['native_seconds'] >= 0
        assert abs(timestamp-native-event['native_seconds']) < 1e-5
        assert event['cost'] < previous_cost
        previous_time, previous_cost = timestamp, event['cost']
        raw_cost, raw_route = event['cost'], event['route']
        if timestamp <= began+budget and event['cost'] < cost:
            cost, route, origin = event['cost'], event['route'], 'phase_event'
    assert (raw['cost'], raw['route']) == (raw_cost, raw_route)
    check(data, raw['route'], raw['cost']); checked += 1
    assert (record['cost'], record['route'], record['origin']) == (cost, route, origin)
    check(data, record['route'], record['cost']); checked += 1
    builds = record['cache_build_events']
    footprint = record['cache_retained_python_bytes']
    assert type(footprint) is int and footprint >= 0
    if arm == 'uncached':
        assert builds == [] and footprint == 0
    elif raw['lkh_calls']:
        assert len(builds) == 1 and footprint > 0
        event = builds[0]
        assert all(math.isfinite(event[k]) for k in ('call_started', 'call_finished', 'build_seconds', 'allowance_seconds'))
        assert native <= event['call_started'] <= event['call_finished'] <= began+elapsed+1e-5
        assert 0 <= event['build_seconds'] <= event['call_finished']-event['call_started']+1e-5
        assert 0 < event['allowance_seconds'] <= .2
    else:
        assert builds == [] and footprint == 0
    return dict(status='PASS', routes_checked=checked,
                eligible_phase_events=sum(e['native_timestamp'] <= began+budget for e in raw['trace']))
