"""Reproduce native HiGHS outer-cut checks against the full SCF formulation."""
import json
from pathlib import Path
import numpy as np
from fragment_repair import FragmentRepair
from pilot import Repair

def main():
    rng=np.random.default_rng(671); cuts=0
    for trial in range(30):
        points=rng.normal(size=(8,2))
        costs=np.ceil(10*np.linalg.norm(points[:,None]-points[None,:],axis=2))
        d=dict(visits=np.arange(8),mandatory={0},costs=costs,prizes=np.ones(8),
               quota=4,neighborhoods=[{2,3}])
        route=[0,1,2,3,4,5]
        u=set(map(int,rng.choice(np.arange(1,8),size=5,replace=False)))
        reference=Repair(d).solve(route,u,10,'route_context')
        for warm in (False,True):
            row=FragmentRepair(d,'highspy',warm,'cuts').solve(route,u,10,'route_context')
            assert row['proven_optimal'] and reference['proven_optimal']
            assert abs(row['cost']-reference['cost'])<1e-6,(trial,row,reference)
            cuts+=row['cuts_added']
    result=dict(native_cut_backend_checks=60,cuts_added=cuts,status='PASS')
    Path(__file__).with_name('native_cut_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))

if __name__=='__main__': main()
