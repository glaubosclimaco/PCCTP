"""Enumerate every small feasible tour independently of the solver model."""
import sys,json,random
from pathlib import Path
from itertools import combinations,permutations
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'experiments/ml_lns'))
import numpy as np
from highs_warm_start import solve

def edges(route):return {tuple(sorted((route[i-1],v))) for i,v in enumerate(route)}
if __name__=='__main__':
    rng=random.Random(20261019);n=6;c=np.zeros((n,n))
    for i in range(n):
        for j in range(i):c[i,j]=c[j,i]=rng.randint(1,100)
    d=dict(visits=np.arange(n),mandatory={0},costs=c,prizes=np.ones(n),quota=3,neighborhoods=[{1,2},{3,4,5}])
    initial=list(range(n));other=[0,2,4,1,5,3];pools=[edges(initial),edges(initial)|edges(other),set(combinations(range(n),2))];rows=[]
    for pool in pools:
        values=[]
        for k in range(3,n+1):
            for chosen in combinations(range(1,n),k-1):
                S={0,*chosen}
                if not all(S&cover for cover in d['neighborhoods']):continue
                for order in permutations(chosen):
                    route=[0,*order]
                    if edges(route)<=pool:values.append(sum(c[route[i-1],v] for i,v in enumerate(route)))
        optimum=min(values);r=solve(d,initial,True,seconds=5.,allowed_edges=pool)
        assert r['proven_optimal'] and r['proof_scope']=='restricted edge pool' and r['cost']==optimum and edges(r['route'])<=pool
        rows.append(dict(pool_edges=len(pool),enumerated_optimum=optimum,solver_cost=r['cost']))
    try:solve(d,initial,True,seconds=5.,allowed_edges=set())
    except ValueError:pass
    else:raise AssertionError('Invalid start pool accepted')
    (ROOT/'experiments/lkh/edge_pool_verification.json').write_text(json.dumps(dict(status='PASS',scope='three six-vertex edge pools exhaustively enumerated; invalid pool rejected',rows=rows),indent=2)+'\n');print('PASS three pools and start rejection')
