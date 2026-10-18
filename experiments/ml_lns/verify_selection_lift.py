"""Independent subset enumeration for the route_context selection-lifting lemma."""
from itertools import combinations,permutations
import json
from pathlib import Path
import numpy as np
from selection_lift import lift
from pilot import validate

def main():
    rng=np.random.default_rng(812); checked=0; outside_sizes=set()
    for trial in range(80):
        n=8; incumbent=[0,1,2,3]
        u=({1,2,3,4} if trial%3==0 else {2,3,4,5} if trial%3==1
           else set(map(int,rng.choice(np.arange(1,n),size=4,replace=False))))
        points=rng.normal(size=(n,2)); costs=np.ceil(10*np.linalg.norm(points[:,None]-points[None,:],axis=2))
        prizes=rng.integers(0,5,size=n)
        neighborhoods=[]
        for _ in range(3):
            ns=set(map(int,rng.choice(n,size=3,replace=False))); ns.add(int(rng.choice(incumbent)))
            neighborhoods.append(ns)
        d=dict(visits=np.arange(n),mandatory={0},costs=costs,prizes=prizes,
               quota=float(sum(prizes[incumbent])*.5),neighborhoods=neighborhoods)
        validate(d,incumbent); outside=set(incumbent)-u
        outside_sizes.add(len(outside))
        for k in range(len(u)+1):
            for subset in combinations(sorted(u),k):
                selection=outside|set(subset)
                if len(selection)<3 or sum(prizes[v] for v in selection)<d['quota']-1e-6: continue
                if not all(selection&ns for ns in neighborhoods): continue
                for budget in (len(u),2*len(u)):
                    route=lift(d,incumbent,u,set(subset),seed=trial,extra_edges=budget)
                    assert set(route)==selection
                    validate(d,route); checked+=1
    assert {1,2,3}<=outside_sizes
    tightness=[]
    for k in range(1,7):
        incumbent=[]
        for i in range(k): incumbent.extend([i,k+2+i])
        incumbent.extend([k,k+1])
        old={tuple(sorted((incumbent[i-1],v))) for i,v in enumerate(incumbent)}
        retained=set(range(k+2))
        fixed={e for e in old if set(e)<=retained}
        assert len(fixed)==2
        minimum=float('inf'); tours=0
        for order in permutations(range(1,k+2)):
            if order[0]>order[-1]: continue
            route=[0,*order]
            edges={tuple(sorted((route[i-1],v))) for i,v in enumerate(route)}
            if not fixed<=edges: continue
            minimum=min(minimum,len(edges-old)); tours+=1
        assert minimum==k,(k,minimum)
        tightness.append(dict(k=k,minimum_required_extra_edges=int(minimum),retained_tours_enumerated=tours))
    result=dict(status='PASS',random_contexts=80,feasible_assignments_lifted=checked,
                retained_outside_sizes=sorted(outside_sizes),extra_edge_budgets=['k','2k'],
                worst_case_lower_bound_checks=tightness,
                scope='Finite independent assignment tests; proof reviewed separately')
    Path(__file__).with_name('selection_lift_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))

if __name__=='__main__': main()
