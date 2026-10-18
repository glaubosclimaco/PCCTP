"""Independent set-based verification of the conditional assignment descriptors."""
from itertools import combinations
import json
import numpy as np
from exchange_signature import signature

def main():
    rng=np.random.default_rng(1901)
    checked=0
    for trial in range(60):
        n=8
        points=rng.normal(size=(n,2))
        costs=np.ceil(10*np.linalg.norm(points[:,None]-points[None,:],axis=2))
        route=list(range(6)); selected=set(route)
        u=set(map(int,rng.choice(np.arange(1,n),size=4,replace=False)))
        neighborhoods=[]
        for _ in range(4):
            ns=set(map(int,rng.choice(n,size=3,replace=False)))
            ns.add(int(rng.choice(route)))
            neighborhoods.append(ns)
        prizes=rng.integers(0,5,size=n)
        quota=float(sum(prizes[route])*.75)
        d=dict(visits=np.arange(n),mandatory={0},costs=costs,prizes=prizes,
               quota=quota,neighborhoods=neighborhoods)
        actual=signature(d,route,u)
        outside=selected-u; states=[]
        for size in range(len(u)+1):
            for subset in combinations(sorted(u),size):
                chosen=outside|set(subset)
                cov=all(chosen & ns for ns in neighborhoods)
                prize=sum(prizes[v] for v in chosen)>=quota-1e-8
                valid=cov and prize and len(chosen)>=3
                states.append((set(subset),cov,prize,valid,len(chosen^selected)))
        feasible=[s for s in states if s[3]]
        assert feasible
        expected=[len(feasible)/16,sum(s[1] for s in states)/16,
                  sum(s[2] for s in states)/16,min(len(s[0]) for s in feasible)/4,
                  sum(s[4]==1 for s in feasible)/4,
                  sum(s[4]==2 for s in feasible)/6]
        alternatives=[s[4] for s in feasible if s[4]]
        expected.append(min(alternatives)/4 if alternatives else 1.)
        entropy=0.
        for v in u:
            p=sum(v in s[0] for s in feasible)/len(feasible)
            if 0<p<1: entropy-=p*np.log2(p)+(1-p)*np.log2(1-p)
        expected.append(entropy/4)
        assert np.allclose(actual[:8],expected,atol=1e-12),(actual,expected)
        # Every feasible assignment's additive saving evaluated independently.
        values={}
        for v in u:
            if v in selected:
                i=route.index(v); a,b=route[i-1],route[(i+1)%len(route)]
                values[v]=costs[a,v]+costs[v,b]-costs[a,b]
            else:
                values[v]=min(costs[a,v]+costs[v,b]-costs[a,b]
                              for a,b in zip(route[-1:]+route[:-1],route))
        saving=max(sum(values[v] for v in (selected&u)-s[0])-
                   sum(values[v] for v in s[0]-selected) for s in feasible)
        route_value=sum(costs[a,b] for a,b in zip(route[-1:]+route[:-1],route))
        assert abs(actual[8]-saving/route_value)<1e-12
        assert abs(actual[9]-expected[0]/(expected[1]*expected[2]))<1e-12
        checked+=1
    print(json.dumps(dict(status='PASS',random_cases=checked,assignments=checked*16,
                         descriptors_per_case=10)))

if __name__=='__main__': main()
