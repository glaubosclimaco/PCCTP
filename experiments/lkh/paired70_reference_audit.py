"""Independent raw-export audit of a completed or censored70-iteration reference."""
import hashlib,json,math
from pathlib import Path
from raw_route_audit import check


def audit_reference(r,d,folder,seed):
    assert r['seed']==seed and r['input_sha256']==d['input_hash']
    native=r['native'];assert native['budget_seconds']==0 and native['elite_count_bounded']
    checked=0;payload=None
    if native['payload_source']:
        source=Path(folder)/native['payload_source']
        assert hashlib.sha256(source.read_bytes()).hexdigest()==r['payload_sha256']
        payload=json.loads(source.read_text())
        assert payload['seed']==seed and payload['budget_seconds']==0
        check(d,payload['route'],payload['cost']);checked+=1
        previous_iteration=0;previous_time=0.;previous_cost=float('inf')
        for event in payload['snapshots']:
            assert previous_iteration<event['iteration']<=70
            assert math.isfinite(event['seconds']) and math.isfinite(event['cost'])
            assert previous_time<=event['seconds'] and event['cost']<=previous_cost
            previous_iteration,previous_time,previous_cost=event['iteration'],event['seconds'],event['cost']
            check(d,event['route'],event['cost']);checked+=1
        assert r['observed_cost']==payload['cost'] and r['observed_iterations']==payload['iterations']
        assert r['mining_calls']==payload['mining_calls']
    else:assert r['payload_sha256'] is None
    complete=bool(payload and payload['completed'] and payload['iterations']==70
                  and not native['hard_timeout'] and native['exit_code']==0)
    assert (r['status']=='COMPLETED70')==complete
    if complete:
        assert payload['snapshots'][-1]['iteration']==70
        assert payload['snapshots'][-1]['cost']==payload['cost']
        assert (r['cost'],r['route'])==(payload['cost'],payload['route'])
        duration=native['process_wall_seconds']
        assert math.isfinite(duration) and duration>0 and r['wrapper_wall_seconds']>=duration
        assert r['reference_duration_seconds']==duration
        check(d,r['route'],r['cost']);checked+=1
    else:
        assert r['cost'] is None and r['route']==[] and r['reference_duration_seconds'] is None
    return dict(routes_checked=checked,completed70=complete,
                duration=r['reference_duration_seconds'],cost=r['cost'])
