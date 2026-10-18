"""Check contraction against independently enumerated tours and original SCF."""
import json
import numpy as np
from fragment_repair import FragmentRepair
from pilot import Repair,validate
import verify_repair as exhaustive

def main():
    exhaustive.Repair=FragmentRepair
    count=0
    for released in ({2},{2,3,4},{2,3,4,5,6}):
        for operator in ('incident','random_extra','route_context','geometric_extra','route_random'):
            exhaustive.check(released,operator); count+=1
    rng=np.random.default_rng(193)
    for trial in range(20):
        points=rng.normal(size=(7,2))
        costs=np.ceil(10*np.linalg.norm(points[:,None]-points[None,:],axis=2))
        d=dict(visits=np.arange(7),mandatory={0},costs=costs,prizes=np.ones(7),quota=3,
               neighborhoods=[{0,1},{2,3,4}])
        route=[0,1,2,3,4]
        released=set(map(int,rng.choice(np.arange(1,7),size=3,replace=False)))
        operator=('incident','route_context','random_extra','geometric_extra','route_random')[trial%5]
        full=Repair(d).solve(route,released,5.,operator,seed=trial)
        small=FragmentRepair(d).solve(route,released,5.,operator,seed=trial)
        assert full['proven_optimal'] and small['proven_optimal']
        assert abs(full['cost']-small['cost'])<1e-7,(trial,full,small)
        count+=1
    # Freeing only absent vertices leaves an immutable incumbent cycle.
    for operator in ('incident','route_context'):
        a=Repair(d).solve(route,{5,6},5.,operator)
        b=FragmentRepair(d).solve(route,{5,6},5.,operator)
        assert abs(a['cost']-b['cost'])<1e-7
        validate(d,[int(v) for v in b['route']]); count+=1
    result=dict(status='PASS',contracted_cases=count,exhaustive_cases=15,
                random_full_SCF_comparisons=20,absent_only_cases=2)
    (__import__('pathlib').Path(__file__).with_name('fragment_verification.json')).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))

if __name__=='__main__': main()
