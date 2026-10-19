"""Deterministic all-visit cycle followed by feasible optional deletions."""
import time
from experiments.ml_lns.pilot import validate

def construct(d,seconds=.25):
    began=time.perf_counter();end=began+seconds
    remaining=set(range(len(d['visits'])));route=[min(d['mandatory'])];remaining.remove(route[0])
    while remaining:
        nxt=min(remaining,key=lambda v:(d['costs'][route[-1],v],v));route.append(nxt);remaining.remove(nxt)
    # Always finish and validate the complete cycle, even if construction exceeds
    # its soft cap; callers must check their outer acceptance deadline.
    best=validate(d,route)
    while time.perf_counter()<end:
        options=[]
        for i,v in enumerate(route):
            if v in d['mandatory'] or len(route)<=3:continue
            saving=d['costs'][route[i-1],v]+d['costs'][v,route[(i+1)%len(route)]]-d['costs'][route[i-1],route[(i+1)%len(route)]]
            if saving>0:options.append((-saving,v))
        changed=False
        for _,v in sorted(options):
            if time.perf_counter()>=end:break
            candidate=[u for u in route if u!=v]
            try:value=validate(d,candidate)
            except ValueError:continue
            if value<best:route,best,changed=candidate,value,True;break
        if not changed:break
    return dict(route=route,cost=validate(d,route),seconds=time.perf_counter()-began)
