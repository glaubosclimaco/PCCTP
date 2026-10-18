"""Exact-rational audit of a mixed coverage/prize rank witness; no solver needed."""
from fractions import Fraction as F
from itertools import combinations, permutations
import json
from pathlib import Path


def subsets(items):
    items = sorted(items)
    for n in range(len(items) + 1):
        yield from (set(s) for s in combinations(items, n))


T = {0, 1, 2, 3}
R = {4, 5, 6}
V = T | R
root = 0
prize = {4: 1, 5: 1, 6: 2}
quota = 2
demand_neighbors = {4, 5}
y = {v: F(1) if v in T else F(1, 2) for v in V}
x = {frozenset((i, j)): F(5, 12) if i in T and j in T
     else F(1, 4) if (i in T) != (j in T) else F(0)
     for i, j in combinations(sorted(V), 2)}
checks = {}


def check(name, predicate):
    checks[name] = bool(predicate)
    assert predicate, name


check('degree', all(sum(value for edge, value in x.items() if i in edge) == 2*y[i] for i in V))
check('edge_count', sum(x.values()) == len(T) + sum(y[i] for i in R))
check('edge_linking', all(0 <= value <= min(y[i] for i in edge) for edge, value in x.items()))
check('prize', sum(prize[i]*y[i] for i in R) >= quota)
check('coverage', sum(y[i] for i in demand_neighbors) >= 1)

cut_count = 0
for S in subsets(V - {root}):
    if not S:
        continue
    cut = sum(value for edge, value in x.items() if len(edge & S) == 1)
    rhs = F(2) if S & T or demand_neighbors <= S else max(2*y[i] for i in S)
    check(f'cutset_{sorted(S)}', cut >= rhs)
    cut_count += 1

kc_count = extended_count = 0
for A in subsets(R):
    residual = quota - sum(prize[i] for i in A)
    if residual > 0:
        check(f'KC_{sorted(A)}', sum(min(prize[i], residual)*y[i] for i in R-A) >= residual)
        kc_count += 1
for C in subsets(R):
    if C and sum(prize[i] for i in R-C) < quota:
        extension = {i for i in R-C if prize[i] >= max(prize[j] for j in C)}
        check(f'cover_{sorted(C)}', sum(y[i] for i in C) >= 1)
        check(f'extended_cover_{sorted(C)}', sum(y[i] for i in C | extension) >= 1+len(extension))
        extended_count += 1

# Strong corrected blossom bound, stronger than both manuscript and code forms.
blossom_count = 0
for H in subsets(V):
    crossing = [edge for edge in x if len(edge & H) == 1]
    inside = sum(value for edge, value in x.items() if edge <= H)
    for n in range(1, len(crossing)+1, 2):
        # For fixed cardinality the largest x values maximize the left side.
        maximum = inside + sum(sorted((x[edge] for edge in crossing), reverse=True)[:n])
        check(f'blossom_H{sorted(H)}_teeth{n}', maximum <= sum(y[i] for i in H) + n//2)
        blossom_count += 1

feasible_selections = []
cycle_count = 0
for S in subsets(R):
    if sum(prize[i] for i in S) < quota or not (S & demand_neighbors):
        continue
    feasible_selections.append(sorted(S))
    for order in permutations(sorted((T | S)-{root})):
        if order[0] > order[-1]:
            continue
        cycle = (root,) + order + (root,)
        assert len(S) >= 2  # Proposed valid mixed rank.
        edges = {frozenset(e) for e in zip(cycle, cycle[1:])}
        assert len(edges) == len(T | S)
        cycle_count += 1

result = {
    'status': 'PASS',
    'witness': {'mandatory': sorted(T), 'optional': sorted(R), 'prizes': prize,
                'quota': quota, 'demand_neighbors': sorted(demand_neighbors),
                'optional_y': '1/2', 'terminal_optional_x': '1/4',
                'terminal_terminal_x': '5/12', 'optional_optional_x': '0'},
    'existing_families_passed': {'root_excluding_cutsets': cut_count,
                                'all_knapsack_covers': kc_count,
                                'all_extended_prize_covers': extended_count,
                                'strong_blossom_maxima_by_handle_and_odd_size': blossom_count},
    'feasible_integer_selections': feasible_selections,
    'feasible_undirected_cycles_enumerated': cycle_count,
    'mixed_rank': {'lhs': str(sum(y[i] for i in R)), 'rhs': 2, 'violation': '1/2'},
    'checks': len(checks),
    'interpretation': 'Mixed rank is a rank-1 CG consequence of three pair obligations; generic novelty is not claimed.'
}
output = Path(__file__).with_name('joint_witness_results.json')
output.write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
