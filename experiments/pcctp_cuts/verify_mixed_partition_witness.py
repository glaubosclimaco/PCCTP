"""Exact rational verification; Python standard library only.

Two-vertex blocks prevent confusing partition rank with vertex cardinality.
Checks even the stronger correct degree-substituted blossom inequalities.
"""
from fractions import Fraction as Q
from itertools import combinations, permutations, product
import json
from pathlib import Path

V=set(range(9)); root=0
blocks=[{0}]+[{2*i+1,2*i+2} for i in range(4)]
group={v:k for k,B in enumerate(blocks) for v in B}
E=list(combinations(sorted(V),2))
x={e:Q(5,8) if group[e[0]]==group[e[1]] else
   (Q(1,4) if root in e else Q(1,16)) for e in E}
y={v:Q(1) if v==root else Q(5,8) for v in V}
demands=[blocks[i]|blocks[j] for i,j in combinations(range(1,5),2) if (i,j)!=(3,4)]
prize={v:Q() if v==root else (Q(1,2) if group[v] in (1,2) else Q(1)) for v in V}; quota=Q(3)
cut=lambda S:sum((x[e] for e in E if (e[0] in S)!=(e[1] in S)),Q())
counts=dict(sec=0,coverage_forced=0,prize_forced=0,odd_F_blossoms=0,
            knapsack_cover=0,extended_knapsack_cover=0,residual_knapsack_cover=0,
            cover_intersection=0,feasible_cycles=0,feasible_selections=0)
cost={e:Q(1,10) if group[e[0]]==group[e[1]] else Q(1) for e in E}
metric=dict(cost)
for w,D in enumerate(demands,9):
    for v in V:metric[min(v,w),max(v,w)]=Q(1 if v in D else 2)
for a,b in combinations(range(9,14),2):metric[a,b]=Q(2)
dist=lambda a,b:Q() if a==b else metric[min(a,b),max(a,b)]
triangle_checks=0
for a,b,c in product(range(14),repeat=3):
    triangle_checks+=1; assert dist(a,c)<=dist(a,b)+dist(b,c)
for w,D in enumerate(demands,9):assert {v for v in V if dist(v,w)<=1}==D
hull={v:Q() for v in V}; hull_weight=Q();hull_terms=0
for active in combinations(range(1,5),3):
    active_vertices=set().union(*(blocks[b] for b in active))
    for omitted in active_vertices:
        weight=Q(1,24);selected={root}|(active_vertices-{omitted})
        assert all(selected&D for D in demands)
        assert sum(prize[v] for v in selected)>=quota
        hull_weight+=weight;hull_terms+=1
        for v in selected:hull[v]+=weight
assert hull_weight==1 and hull==y
minimum_cycle_cost=None
ranks={'coverage_only':None,'prize_only':None,'joint':None}
for v in V:
    assert sum(x[e] for e in E if v in e)==2*y[v]
for i,j in E: assert x[i,j]<=min(y[i],y[j])
for D in demands: assert sum(y[v] for v in D)>=1
assert sum(prize[v]*y[v] for v in V)>=quota
for mask in range(1,1<<9):
    H={v for v in V if mask>>v&1}
    if H==V:continue
    U=H if root not in H else V-H
    for v in U:
        counts['sec']+=1; assert cut(H)>=2*y[v]
    if any(D<=U for D in demands):
        counts['coverage_forced']+=1; assert cut(H)>=2
    if any(D<=H for D in demands) or root in H:
        for D in demands+[{root}]:
            A=H&D
            val=sum(x[e] for e in E if ((e[0] in H)!=(e[1] in H)) or ((e[0] in A)!=(e[1] in A)))
            counts['cover_intersection']+=1
            assert val>=2,('CI',H,D,val)
    if sum(prize[v] for v in V-U)<quota:
        counts['prize_forced']+=1; assert cut(H)>=2
    cross=[e for e in E if (e[0] in H)!=(e[1] in H)]
    # Maximize sum(2*x_e-1) over ALL odd subsets F by parity DP.
    # This proves all such inequalities, including nonmatching F.
    dp=[Q(),None]
    for e in cross:
        weight=2*x[e]-1; old=dp[:]
        dp=[max(z for z in (old[0],None if old[1] is None else old[1]+weight) if z is not None),
            max(z for z in (old[1],old[0]+weight) if z is not None)]
    counts['odd_F_blossoms']+=1<<(len(cross)-1)
    assert cut(H)>=1+dp[1],(H,cut(H),dp)
for mask in range(1<<8):
    S={i+1 for i in range(8) if mask>>i&1}
    if sum(prize[v] for v in (V-{root})-S)<quota:
        counts['knapsack_cover']+=1
        assert sum(y[v] for v in S)>=1
        threshold=max(prize[v] for v in S)
        extension={v for v in V-{root}-S if prize[v]>=threshold}
        counts['extended_knapsack_cover']+=1
        assert sum(y[v] for v in S|extension)>=1+len(extension)
    if sum(prize[v] for v in S)<quota:
        residual=quota-sum(prize[v] for v in S)
        counts['residual_knapsack_cover']+=1
        assert sum(min(prize[v],residual)*y[v] for v in V-{root}-S)>=residual
    visited=S|{root}
    nblocks=len({group[v] for v in S})
    cover_ok=all(visited&D for D in demands)
    prize_ok=sum(prize[v] for v in visited)>=quota
    for key,ok in [('coverage_only',cover_ok),('prize_only',prize_ok),('joint',cover_ok and prize_ok)]:
        if ok:ranks[key]=nblocks if ranks[key] is None else min(ranks[key],nblocks)
    if len(visited)<3 or not all(visited&D for D in demands):continue
    if sum(prize[v] for v in visited)<quota:continue
    counts['feasible_selections']+=1
    for order in permutations(sorted(S)):
        if order[0]>order[-1]:continue # reversal symmetry
        route=(root,)+order
        crossings=sum(group[route[i]]!=group[route[(i+1)%len(route)]] for i in range(len(route)))
        counts['feasible_cycles']+=1
        assert crossings>=4
        cycle_cost=sum(dist(route[i],route[(i+1)%len(route)]) for i in range(len(route)))
        minimum_cycle_cost=cycle_cost if minimum_cycle_cost is None else min(minimum_cycle_cost,cycle_cost)
lhs=sum(x[e] for e in E if group[e[0]]!=group[e[1]])
assert lhs==Q(7,2)<4
fractional_cost=sum(x[e]*cost[e] for e in E)
assert fractional_cost==Q(15,4)<minimum_cycle_cost==Q(41,10)
assert ranks=={'coverage_only':2,'prize_only':2,'joint':3}
result={'status':'all assertions passed','vertices':9,'blocks':[sorted(B) for B in blocks],
        'quota':str(quota),'partition_rank':3,'partition_lhs':str(lhs),'partition_rhs':4,
        'separate_and_joint_ranks':ranks,
        'optional_vertex_mass':str(sum(y[v] for v in V-{root})),
        'full_selection_hull':{'terms':hull_terms,'weight':str(hull_weight),'exact_marginals':True},
        'metric':{'triangle_checks':triangle_checks,'all_passed':True,'fractional_cost':str(fractional_cost),
                  'enumerated_optimal_cycle_cost':str(minimum_cycle_cost)},
        'blossom_verification':'all odd F certified by exact parity DP, not enumerated individually',
        'checks':counts,'x':{f'{i}-{j}':str(x[i,j]) for i,j in E},
        'y':{str(v):str(y[v]) for v in V},'coverage_neighborhoods':[sorted(D) for D in demands]}
out=Path(__file__).with_name('mixed_rank_witness_result.json')
out.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ('x','y')},indent=2))
