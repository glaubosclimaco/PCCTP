"""Mocked adapters exercise reference/candidate boundary without native solvers."""
import hashlib,json,tempfile
from pathlib import Path
import paired70_reference as ref
import paired70_hybrid as hybrid
from paired70_candidate_audit import audit_candidate
from paired70_reference_audit import audit_reference
from raw_route_audit import check

if __name__=='__main__':
    d=dict(n=4,types=[1,0,0,0],prizes=[1]*4,quota=1,radius=0,
           matrix=[0 if i==j else 1 for i in range(4) for j in range(4)],
           visits=list(range(4)),input_hash='synthetic')
    clock=dict(now=100.,complete=True);original_clock=ref.time.clock_gettime
    saved=(ref.run_one,ref.load,hybrid.run_one,hybrid.read_instance,
           hybrid.construct,hybrid.validate,hybrid.search)
    def native(name,seed,policy,budget,target,bounded,parent_started,hard_timeout_seconds):
        assert name.startswith('synthetic') and policy=='original' and bounded
        clock['now']=parent_started+2.1
        if budget==0:
            raw=Path(target).with_suffix('.raw');iteration=70 if clock['complete'] else 35
            payload=dict(seed=seed,budget_seconds=0,completed=clock['complete'],iterations=iteration,route=[0,1,2,3],
                         cost=4,mining_calls=1,snapshots=[dict(iteration=iteration,seconds=1.9,route=[0,1,2,3],cost=4)])
            raw.write_text(json.dumps(payload))
            return dict(budget_seconds=0,elite_count_bounded=True,payload_source=raw.name,
                        process_wall_seconds=2.,hard_timeout=False,exit_code=0)
        assert budget==2.
        return dict(budget_seconds=budget,feasible_before_deadline=True,cost=3,route=[0,1,2],
                    anytime=[dict(seconds=.5,cost=3,route=[0,1,2])])
    def fallback(data):
        clock['now']+=.1
        return dict(cost=4,route=[0,1,2,3])
    def forbidden_search(*args,**kwargs):raise AssertionError('No search allowed after spent horizon')
    try:
        ref.time.clock_gettime=lambda _:clock['now']
        ref.run_one=native;ref.load=lambda _:d
        hybrid.run_one=native;hybrid.read_instance=lambda _:d;hybrid.construct=fallback
        hybrid.validate=lambda _,route:check(d,route,len(route));hybrid.search=forbidden_search
        with tempfile.TemporaryDirectory(prefix='pcctp-adapters-',dir='/private/tmp') as tmp:
            folder=Path(tmp);reference=ref.run('synthetic_complete',81,folder)
            assert reference['status']=='COMPLETED70' and reference['cost']==4
            assert audit_reference(reference,d,folder,81)['completed70']
            budget=reference['reference_duration_seconds'];assert budget==2.
            candidate=hybrid.run('synthetic_complete',81,'hybrid10',budget,folder)
            assert candidate['initialization_budget_seconds']==2.
            assert candidate['cost']==3 and candidate['wall_seconds']>budget
            audit_candidate(candidate,d,budget)
            clock['complete']=False
            incomplete=ref.run('synthetic_incomplete',82,folder)
            assert incomplete['status']=='REFERENCE_INCOMPLETE'
            assert incomplete['cost'] is None and incomplete['reference_duration_seconds'] is None
            assert not audit_reference(incomplete,d,folder,82)['completed70']
    finally:
        ref.time.clock_gettime=original_clock
        (ref.run_one,ref.load,hybrid.run_one,hybrid.read_instance,
         hybrid.construct,hybrid.validate,hybrid.search)=saved
    report=dict(status='PASS',scope='mocked adapter integration only; no native solver, LKH or benchmark run',
        checks=['full70 duration transferred','initial phase clamped to short T',
                'timely incumbent retained after late process return','no search after horizon',
                'incomplete reference supplies no comparison budget'],
        source_hashes={n:hashlib.sha256(Path(__file__).with_name(n).read_bytes()).hexdigest()
                       for n in ('paired70_reference.py','paired70_hybrid.py','paired70_candidate_audit.py',
                                 'paired70_reference_audit.py')})
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report))
