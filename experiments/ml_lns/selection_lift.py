"""Constructively lift a feasible conditional selection into route_context.

Requires a complete undirected graph, optional U, at least one mandatory vertex,
and the current route_context expansion (all incident edges plus 2|U| extras).
Does not assert optimal routing cost or a new result in the literature.
"""
from itertools import combinations
import numpy as np
from pilot import free_edge_set,validate

def lift(d,incumbent,released,assignment,seed=42,extra_edges=None):
    released=set(released); assignment=set(assignment)
    if not released or released & d['mandatory']: raise ValueError('Release optional vertices only')
    if not assignment<=released: raise ValueError('Assignment outside release set')
    route=[v for v in incumbent if v not in released]
    selected=set(route)|assignment
    if len(selected)<3 or not d['mandatory']<=selected: raise ValueError('Invalid selection size or mandatory set')
    if sum(d['prizes'][list(selected)])<d['quota']-1e-6: raise ValueError('Invalid prize selection')
    if not all(selected & ns for ns in d['neighborhoods']): raise ValueError('Invalid coverage selection')
    free=free_edge_set(d,incumbent,released,'route_context',seed=seed,extra_edges=extra_edges)
    remaining=set(assignment)
    if len(route)==1:
        a=route[0]
        b,c=min(combinations(sorted(remaining),2),key=lambda pair:
                (d['costs'][a,pair[0]]+d['costs'][pair[0],pair[1]]+d['costs'][pair[1],a],pair))
        route=[a,b,c]; remaining-={b,c}
    elif len(route)==2:
        a,b=route
        v=min(remaining,key=lambda v:(d['costs'][a,v]+d['costs'][v,b],v))
        route=[a,v,b]; remaining.remove(v)
    while remaining:
        choices=[]
        for i,b in enumerate(route):
            a=route[i-1]
            if tuple(sorted((a,b))) not in free: continue
            for v in remaining:
                choices.append((float(d['costs'][a,v]+d['costs'][v,b]-d['costs'][a,b]),v,i))
        if not choices: raise AssertionError('No free insertion gap for a feasible selection')
        _,v,i=min(choices); route.insert(i,v); remaining.remove(v)
    validate(d,route)
    original={tuple(sorted((incumbent[i-1],v))) for i,v in enumerate(incumbent)}
    edges={tuple(sorted((route[i-1],v))) for i,v in enumerate(route)}
    if set(route)!=selected: raise AssertionError('Lift changed selection')
    if any(e not in free and ((e in edges)!=(e in original)) for e in combinations(range(len(d['visits'])),2)):
        raise AssertionError('Lift violated an edge fixing')
    return route

def best_lift(d,incumbent,released,limit=3,seed=42,keep_incumbent=True):
    """Lift top additive-proxy assignments; retain incumbent as fallback.

    Exact feasibility enumeration, heuristic choice and routing. No claim of
    optimal repair, optimal selection under routing cost, or additive cost bound.
    """
    u=sorted(released); k=len(u)
    if not 1<=k<=12: raise ValueError('Enumeration limited to 1..12 released vertices')
    outside=set(incumbent)-set(u); selected=set(incumbent)
    masks=((np.arange(1<<k)[:,None]>>np.arange(k))&1).astype(np.int8)
    valid=np.ones(len(masks),dtype=bool)
    for ns in d['neighborhoods']:
        if not ns&outside: valid &= masks@np.array([v in ns for v in u],dtype=int)>=1
    valid &= masks@d['prizes'][u]+sum(d['prizes'][list(outside)])>=d['quota']-1e-6
    valid &= masks.sum(axis=1)+len(outside)>=3
    contribution=[]
    for v in u:
        if v in selected:
            i=incumbent.index(v); a,b=incumbent[i-1],incumbent[(i+1)%len(incumbent)]
            contribution.append(d['costs'][a,v]+d['costs'][v,b]-d['costs'][a,b])
        else:
            contribution.append(min(d['costs'][incumbent[i-1],v]+d['costs'][v,b]-d['costs'][incumbent[i-1],b]
                                    for i,b in enumerate(incumbent)))
    current=np.array([v in selected for v in u],dtype=int)
    saving=(current-masks)@np.array(contribution)
    indices=sorted(np.flatnonzero(valid),key=lambda i:(-saving[i],int(i)))[:limit]
    best=incumbent[:]; cost=validate(d,best) if keep_incumbent else float('inf')
    for i in indices:
        assignment={v for j,v in enumerate(u) if masks[i,j]}
        candidate=lift(d,incumbent,set(u),assignment,seed)
        value=validate(d,candidate)
        if value<cost: best=candidate; cost=value
    if not np.isfinite(cost): raise AssertionError('No feasible lifted assignment')
    return dict(route=best,cost=cost,feasible_assignments=int(valid.sum()),assignments_lifted=len(indices))
