"""Compare SCF repair to exhaustive feasible-cycle enumeration on a tiny case."""
from itertools import combinations, permutations
import numpy as np
from pilot import Repair, validate, free_edge_set

points = np.array([[0, 0], [9, 0], [2, 5], [7, 5], [1, 2], [8, 2], [5, 8]])
costs = np.ceil(np.linalg.norm(points[:, None]-points[None, :], axis=-1))
d = dict(visits=np.arange(7), mandatory={0, 1}, prizes=np.array([0, 0, 1, 2, 3, 1, 2]),
         quota=5, neighborhoods=[{2, 3}, {4, 5}, {1, 6}], costs=costs)
incumbent = [0, 1, 2, 3, 4, 5, 6]
inc_edges = {tuple(sorted((incumbent[i-1], v))) for i, v in enumerate(incumbent)}
def check(released, operator):
    free = free_edge_set(d, incumbent, released, operator)
    optimum = float('inf'); checked = 0
    for k in range(2, 6):
        for opt in combinations(range(2, 7), k):
            selected = {0, 1, *opt}
            if not (set(incumbent)-released) <= selected:
                continue
            for order in permutations(sorted(selected-{0})):
                if order[0] > order[-1]:
                    continue
                route = [0, *order]
                edges = {tuple(sorted((route[i-1], v))) for i, v in enumerate(route)}
                if any(e not in free and ((e in edges) != (e in inc_edges))
                       for e in combinations(range(7), 2)):
                    continue
                try:
                    value = validate(d, route)
                except ValueError:
                    continue
                checked += 1; optimum = min(optimum, value)
    result = Repair(d).solve(incumbent, released, 5., operator)
    assert result['proven_optimal'], result
    assert abs(result['cost']-optimum) < 1e-7, (result, optimum)
    print(dict(operator=operator, released=sorted(released), feasible_cycles=checked,
               enumerated_optimum=optimum, repair_optimum=result['cost'], status='PASS'))
    return result['cost']


single = {}
for released in ({2}, {2, 3, 4}, {2, 3, 4, 5, 6}):
    for operator in ('incident', 'random_extra', 'route_context', 'geometric_extra', 'route_random'):
        value = check(released, operator)
        if released == {2}:
            single[operator] = value
assert single['route_context'] < single['incident'], single
print({'shortcut_mechanism_witness': single, 'status': 'PASS'})
