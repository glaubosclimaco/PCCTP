"""Exact conditional selection signatures, not routing-feasibility guarantees.

For k<=12 optional vertices U, enumerate assignments A subset U, holding every
outside visit at its incumbent value. Coverage and quota are jointly checked.
Routing reconnection may rule out additional assignments in the actual repair.
"""
from itertools import combinations
import numpy as np
from pilot import route_cost

NAMES=['joint_feasible_fraction','coverage_feasible_fraction','prize_feasible_fraction',
       'minimum_selected_fraction','single_flip_feasible_fraction',
       'pair_flip_feasible_fraction','minimum_exchange_fraction',
       'conditional_visit_entropy','joint_independent_saving','coverage_prize_dependence']


def signature(d,route,released):
    u=sorted(released); k=len(u)
    if not 1<=k<=12: raise ValueError('Exact signature supports 1<=k<=12')
    selected=set(route); outside=selected-released
    assignments=((np.arange(1<<k)[:,None]>>np.arange(k))&1).astype(np.int8)
    current=np.array([v in selected for v in u],dtype=np.int8)
    residual=[ns for ns in d['neighborhoods'] if not ns & outside]
    if residual:
        incidence=np.array([[v in ns for v in u] for ns in residual],dtype=np.int16)
        coverage=(assignments@incidence.T>=1).all(axis=1)
    else: coverage=np.ones(len(assignments),dtype=bool)
    remaining=max(0.,d['quota']-float(sum(d['prizes'][list(outside)])))
    prize=assignments@d['prizes'][u]>=remaining-1e-8
    size=len(outside)+assignments.sum(axis=1)>=3
    valid=coverage & prize & size
    if not valid.any() or not np.any(np.all(assignments[valid]==current,axis=1)):
        raise AssertionError('Incumbent must belong to its conditional selection space')
    flips=np.abs(assignments-current).sum(axis=1)
    alternatives=flips[valid & (flips>0)]
    marginals=assignments[valid].mean(axis=0)
    entropy=0.
    for p in marginals:
        if 0<p<1: entropy+=-p*np.log2(p)-(1-p)*np.log2(1-p)
    # Label-free additive insertion/deletion surrogate, evaluated only on
    # jointly feasible assignments. It deliberately ignores reconnection overlap.
    contribution=[]
    for v in u:
        if v in selected:
            i=route.index(v); a,b=route[i-1],route[(i+1)%len(route)]
            contribution.append(d['costs'][a,v]+d['costs'][v,b]-d['costs'][a,b])
        else:
            contribution.append(min(d['costs'][route[i-1],v]+d['costs'][v,b]-d['costs'][route[i-1],b]
                                    for i,b in enumerate(route)))
    saving=(current-assignments)@np.array(contribution)
    cov_fraction=coverage.mean(); prize_fraction=prize.mean(); fraction=valid.mean()
    return [float(fraction),float(cov_fraction),float(prize_fraction),
            float(assignments[valid].sum(axis=1).min()/k),
            float(np.count_nonzero(valid & (flips==1))/k),
            float(np.count_nonzero(valid & (flips==2))/max(1,k*(k-1)/2)),
            float(alternatives.min()/k) if len(alternatives) else 1.,
            float(entropy/k),float(saving[valid].max()/max(1.,route_cost(d,route))),
            float(fraction/max(cov_fraction*prize_fraction,1e-12))]
