"""Synthetic nonmetric route/late-event audit; no optimizer or benchmark run."""
import copy,hashlib,json
from pathlib import Path
from paired70_candidate_audit import audit_candidate

if __name__=='__main__':
    n=5;m=[[0 if i==j else 1 for j in range(n)] for i in range(n)]
    for i,j,c in ((0,1,5),(1,2,10),(2,3,10),(3,4,10),(4,0,10)):
        m[i][j]=m[j][i]=c
    d=dict(n=n,types=[1,0,0,0,0],prizes=[1]*n,quota=1,radius=0,
           matrix=[v for row in m for v in row])
    r=dict(mode='hybrid10',budget_seconds=100.,initialization_budget_seconds=10.,
        native_started=100.,wall_seconds=100.4,
        fallback=dict(seconds=.1,cost=45,route=[0,1,2,3,4]),
        native=dict(budget_seconds=10.,feasible_before_deadline=True,cost=16,route=[0,1,2],
                    anytime=[dict(seconds=5.,cost=16,route=[0,1,2])]),
        pre_cost=16,pre_route=[0,1,2],post_offset=10.1,
        post=dict(source_sha256=None,initial_cost=16,native_started=110.1,native_elapsed=90.2,
            lkh_calls=3,lkh_completed=3,lkh_timeouts=0,lkh_failed=0,attempts=3,feasible_repairs=3,
            trace=[dict(native_timestamp=150.,native_seconds=39.9,cost=12,route=[0,2,4]),
                   dict(native_timestamp=200.1,native_seconds=90.,cost=7,route=[0,1,3])],
            cost=7,route=[0,1,3]),cost=12,route=[0,2,4],origin='post_event',feasible=True)
    checked=audit_candidate(r,d,100.);assert checked['cost']==12
    bad=[]
    x=copy.deepcopy(r);x.update(cost=7,route=[0,1,3]);bad.append(('accepted_late_incumbent',x))
    x=copy.deepcopy(r);x['native']['anytime'][0]['seconds']=11.;bad.append(('native_phase_overrun',x))
    x=copy.deepcopy(r);x['budget_seconds']=99.;bad.append(('wrong_reference_duration',x))
    x=copy.deepcopy(r);x['post']['source_sha256']='offline';bad.append(('offline_pool_leak',x))
    x=copy.deepcopy(r);x['post']['trace'][0]['cost']=11.;bad.append(('false_route_cost',x))
    rejected=[]
    for label,record in bad:
        try:audit_candidate(record,d,100.)
        except AssertionError:rejected.append(label)
    assert len(rejected)==len(bad)
    source=Path(__file__).with_name('paired70_candidate_audit.py')
    report=dict(status='PASS',valid_late_return_case=True,corruptions_rejected=rejected,
        routes_checked=checked['routes_checked'],source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        scope='finite synthetic audit checks only; ten-seed main study not launched')
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report))
