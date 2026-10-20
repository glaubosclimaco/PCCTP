"""Synthetic positive/corruption checks; never launches an optimizer."""
import copy,hashlib,json,tempfile
from pathlib import Path
from paired70_reference_audit import audit_reference

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

if __name__=='__main__':
    data=dict(n=4,types=[1,0,0,0],prizes=[1]*4,quota=1,radius=0,
        matrix=[0 if i==j else 1 for i in range(4) for j in range(4)],input_hash='synthetic')
    payload=dict(seed=81,budget_seconds=0,completed=True,iterations=70,
        route=[0,1,2],cost=3,mining_calls=2,snapshots=[
            dict(iteration=1,seconds=.1,route=[0,1,2,3],cost=4),
            dict(iteration=70,seconds=1.8,route=[0,1,2],cost=3)])
    record=dict(seed=81,input_sha256='synthetic',status='COMPLETED70',
        native=dict(budget_seconds=0,elite_count_bounded=True,payload_source='native.raw',
            hard_timeout=False,exit_code=0,process_wall_seconds=2.),
        observed_cost=3,observed_iterations=70,mining_calls=2,
        cost=3,route=[0,1,2],reference_duration_seconds=2.,wrapper_wall_seconds=2.1)
    rejected=[]
    with tempfile.TemporaryDirectory(prefix='pcctp-ref-audit-',dir='/private/tmp') as tmp:
        folder=Path(tmp);source=folder/'native.raw'
        def verify(r,p):
            source.write_text(json.dumps(p));r['payload_sha256']=digest(source)
            return audit_reference(r,data,folder,81)
        assert verify(copy.deepcopy(record),copy.deepcopy(payload))['completed70']
        incomplete=copy.deepcopy(record);censored=copy.deepcopy(payload)
        censored.update(completed=False,iterations=35)
        censored['snapshots'][-1]['iteration']=35
        incomplete.update(status='REFERENCE_INCOMPLETE',cost=None,route=[],
            reference_duration_seconds=None,observed_iterations=35)
        assert not verify(incomplete,censored)['completed70']
        mutations={
            'wrong_seed':lambda r,p:p.update(seed=82),
            'positive_reference_budget':lambda r,p:r['native'].update(budget_seconds=300),
            'wrong_duration':lambda r,p:r.update(reference_duration_seconds=3.),
            'timeout_as_complete':lambda r,p:r['native'].update(hard_timeout=True),
            'nonzero_exit_as_complete':lambda r,p:r['native'].update(exit_code=1),
            'partial_iterations_as_complete':lambda r,p:p.update(iterations=35),
            'missing_final_snapshot':lambda r,p:p['snapshots'].pop(),
            'negative_event_time':lambda r,p:p['snapshots'][0].update(seconds=-1.),
            'nonfinite_event_time':lambda r,p:p['snapshots'][0].update(seconds=float('nan')),
            'wrong_route_cost':lambda r,p:p.update(cost=4),
            'missing_mandatory':lambda r,p:p.update(route=[1,2,3]),
            'short_wrapper_time':lambda r,p:r.update(wrapper_wall_seconds=1.),
            'false_mining_count':lambda r,p:r.update(mining_calls=0),
        }
        for label,mutate in mutations.items():
            r,p=copy.deepcopy(record),copy.deepcopy(payload);mutate(r,p)
            try:verify(r,p)
            except AssertionError:rejected.append(label)
            else:raise AssertionError('Corruption accepted: '+label)
        r=copy.deepcopy(record);verify(r,copy.deepcopy(payload));source.write_text('{}')
        try:audit_reference(r,data,folder,81)
        except AssertionError:rejected.append('changed_payload_hash')
        else:raise AssertionError('Changed raw payload accepted')
    report=dict(status='PASS',valid_cases=2,rejected_corruptions=rejected,
        source_sha256=digest(Path(__file__).with_name('paired70_reference_audit.py')),
        scope='synthetic audit fixtures only; no native solver or PCCTP benchmark execution')
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report))
