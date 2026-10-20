"""Verify exact subset-DP repair by independent SCF and exhaustive tiny tours."""
import json
import hashlib
from pathlib import Path
import numpy as np
from dp_repair import DPRepair
from pilot import Repair

def main():
    rng=np.random.default_rng(20261018); rows=[]
    for trial in range(30):
        n=8; points=rng.normal(size=(n,2)); costs=np.ceil(10*np.linalg.norm(points[:,None]-points[None,:],axis=2))
        d=dict(visits=np.arange(n),mandatory={0,1},costs=costs,prizes=np.arange(n)%3,
               quota=4,neighborhoods=[{2,3},{4,5,6}])
        route=list(range(7)); released=set(map(int,rng.choice(np.arange(n),size=5,replace=False)))
        a=Repair(d).solve(route,released,5,'route_context'); b=DPRepair(d).solve(route,released,5)
        assert a['proven_optimal'] and b['proven_optimal'] and abs(a['cost']-b['cost'])<1e-6,(trial,a,b)
        rows.append(dict(trial=trial,cost=b['cost'],components=b['contracted_nodes'],states=b['dp_processed_states']))
    import verify_repair as exhaustive
    exhaustive.Repair=DPRepair
    for released in ({2},{2,3,4},{2,3,4,5,6}):
        for operator in ('incident','random_extra','route_context','geometric_extra','route_random','fragment_clique'):
            exhaustive.check(released,operator)
    # Vertex zero is selected here: a clock check tied to mask low bits would
    # never trigger. Native DP must still stop and retain a feasible incumbent.
    n=14; costs=np.ones((n,n)); np.fill_diagonal(costs,0)
    d=dict(visits=np.arange(n),mandatory=set(range(n)),costs=costs,prizes=np.zeros(n),quota=0,neighborhoods=[])
    timed=DPRepair(d).solve(list(range(n)),set(range(n)),1e-6)
    assert timed['status']==1 and timed['reason']=='time_limit',timed
    assert timed['cost']==n
    result=dict(status='PASS',independent_SCF_comparisons=30,exhaustive_neighborhoods=18,
                includes_fixed_mandatory_routing_anchors=True,timeout_retains_feasible_route=True,rows=rows,
                native_clock='clock_gettime(CLOCK_MONOTONIC)',
                source_hashes={n:hashlib.sha256(Path(__file__).with_name(n).read_bytes()).hexdigest()
                               for n in ['verify_dp.py','path_dp.cpp','dp_repair.py']})
    Path(__file__).with_name('dp_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='rows'}))

if __name__=='__main__':main()
