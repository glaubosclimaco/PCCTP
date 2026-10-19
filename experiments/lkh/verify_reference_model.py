"""Finite equivalence check: edge-subset model versus independently enumerated cycles."""
from itertools import combinations,permutations
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[2]
if __name__=='__main__':
    rows=[];graphs=0
    for n in range(3,7):
        edges=list(combinations(range(n),2));index={e:i for i,e in enumerate(edges)}
        candidates=[]
        for mask in range(1<<len(edges)):
            degree=[0]*n;chosen=[]
            for k,(a,b) in enumerate(edges):
                if mask>>k&1:degree[a]+=1;degree[b]+=1;chosen.append((a,b))
            graphs+=1
            if degree[0]!=2 or any(d not in (0,2) for d in degree):continue
            selected={i for i,d in enumerate(degree) if d==2};connected=True
            for submask in range(1,1<<(n-1)):
                S={i+1 for i in range(n-1) if submask>>i&1}
                if S&selected and sum((a in S)!=(b in S) for a,b in chosen)<2:connected=False;break
            candidates.append((mask,selected,connected))
        cycles={}
        for size in range(3,n+1):
            for rest in combinations(range(1,n),size-1):
                for order in permutations(rest):
                    route=(0,)+order
                    mask=sum(1<<index[tuple(sorted((route[i-1],v)))] for i,v in enumerate(route))
                    cycles[mask]=set(route)
        configs=[({0},0,[{1,2}],list(range(1,n+1))),({0,1},n,[{n-1}],list(range(1,n+1))),({0},0,[set()],[0]*n),(set(range(n)),0,[],[0]*n)]
        for mandatory,quota,cover,prizes in configs:
            feasible=lambda selected:mandatory<=selected and sum(prizes[i] for i in selected)>=quota and all(selected&S for S in cover)
            model={mask for mask,S,connected in candidates if connected and len(S)>=3 and feasible(S)}
            tour={mask for mask,S in cycles.items() if feasible(S)}
            assert model==tour
            rows.append(dict(n=n,mandatory=sorted(mandatory),quota=quota,coverage_sets=[sorted(S) for S in cover],feasible_cycles=len(model)))
    p=ROOT/'experiments/lkh/reference_model_verification.json';p.write_text(json.dumps(dict(status='PASS',scope='finite check supports but does not replace general integer-equivalence proof',graphs_enumerated=graphs,configurations=len(rows),rows=rows),indent=2)+'\n');print('PASS',graphs,'graphs,',len(rows),'configurations')
