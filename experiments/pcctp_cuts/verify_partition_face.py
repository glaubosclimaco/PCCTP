"""Exact modular affine-rank certificates for two finite PCCTP polytopes.

A nonzero minor modulo a prime is a nonzero integer minor, hence a rational
rank lower bound. Equal upper bounds from explicit equations make ranks exact.
No claim about arbitrary instances or a general facet family is made.
"""
from itertools import combinations, permutations
from pathlib import Path
import json
import time

PRIME = 1000003
E = list(combinations(range(9), 2))
edge_index = {frozenset(edge): k for k, edge in enumerate(E)}
blocks = [{0}] + [{2*i+1, 2*i+2} for i in range(4)]
group = {v: k for k, block in enumerate(blocks) for v in block}
dimension = len(E)+8


class ModularBasis:
    def __init__(self):
        self.rows = {}
    def add(self, point):
        row = [v % PRIME for v in point]
        for column in range(len(row)):
            value = row[column]
            if not value:
                continue
            if column in self.rows:
                basis = self.rows[column]
                for j in range(column, len(row)):
                    row[j] = (row[j]-value*basis[j]) % PRIME
            else:
                inverse = pow(value, -1, PRIME)
                self.rows[column] = [(v*inverse) % PRIME for v in row]
                return True
        return False
    @property
    def rank(self):
        return len(self.rows)


def affine_add(basis, point, first):
    if first is None:
        return point
    basis.add([a-b for a, b in zip(point, first)])
    return first


degree_rows = []
for v in range(9):
    row = [int(v in edge) for edge in E] + [-2*int(i == v) for i in range(1, 9)]
    degree_rows.append(row)
partition_row = [int(group[i] != group[j]) for i,j in E]+[0]*8
equation_basis = ModularBasis()
for row in degree_rows:
    equation_basis.add(row)
degree_rank = equation_basis.rank
equation_basis.add(partition_row)
face_equation_rank = equation_basis.rank


def experiment(mixed):
    pairs = list(combinations(range(1, 5), 2))
    if mixed:
        pairs.remove((3, 4))
    demands = [blocks[i] | blocks[j] for i,j in pairs]
    prize = {v: (1 if group[v] in (1, 2) else 2) for v in range(1,9)} if mixed else {v:1 for v in range(1,9)}
    quota = 6 if mixed else 1
    # Check validity on every feasible selection before early stopping rank work.
    valid_selection_count = 0
    for mask in range(256):
        S = {i+1 for i in range(8) if mask >> i & 1}
        if len(S)>=2 and sum(prize[i] for i in S)>=quota and all(S & D for D in demands):
            assert len({group[i] for i in S}) >= 3
            valid_selection_count += 1
    all_basis, tight_basis = ModularBasis(), ModularBasis()
    first_all = first_tight = None
    count = tight_count = selections = 0
    started = time.monotonic()
    all_generators, tight_generators = [], []
    def result():
        full_upper = dimension-degree_rank
        face_upper = dimension-face_equation_rank
        exact = all_basis.rank == full_upper and tight_basis.rank == face_upper
        return {
            'case': 'mixed_prize_coverage' if mixed else 'all_pair_coverage',
            'coordinate_count': dimension, 'prime': PRIME,
            'degree_equations_rank': degree_rank,
            'degree_plus_partition_equations_rank': face_equation_rank,
            'polytope_affine_rank_lower': all_basis.rank,
            'polytope_affine_dimension_upper': full_upper,
            'tight_face_affine_rank_lower': tight_basis.rank,
            'tight_face_affine_dimension_upper': face_upper,
            'facet_certified_for_this_instance': exact,
            'cycles_processed': count, 'tight_cycles_processed': tight_count,
            'feasible_selections_processed': selections,
            'all_feasible_selections_checked_for_partition_validity': valid_selection_count,
            'validity_argument': 'Every feasible selection meets at least three optional blocks; its cycle also meets the mandatory block, hence crosses the partition at least four times.',
            'early_stop_at_matching_bounds': exact,
            'runtime_seconds': round(time.monotonic()-started, 3),
            'all_affine_generating_routes': all_generators,
            'face_affine_generating_routes': tight_generators,
            'prizes': prize, 'quota': quota,
            'coverage_neighborhoods': [sorted(D) for D in demands]
        }
    for mask in range(256):
        selected = [i+1 for i in range(8) if mask >> i & 1]
        S = set(selected)
        if len(selected)<2 or sum(prize[i] for i in selected)<quota or not all(S & D for D in demands):
            continue
        selections += 1
        for order in permutations(selected):
            if order[0]>order[-1]:
                continue
            route = (0,)+order
            point = [0]*dimension
            crossing = 0
            for i in range(len(route)):
                a,b = route[i],route[(i+1) % len(route)]
                point[edge_index[frozenset((a,b))]] = 1
                crossing += group[a] != group[b]
            for i in selected:
                point[len(E)+i-1] = 1
            assert crossing >= 4
            count += 1
            previous = all_basis.rank
            first_all = affine_add(all_basis, point, first_all)
            if previous != all_basis.rank or len(all_generators)==0:
                all_generators.append(list(route))
            if crossing == 4:
                tight_count += 1
                previous = tight_basis.rank
                first_tight = affine_add(tight_basis, point, first_tight)
                if previous != tight_basis.rank or len(tight_generators)==0:
                    tight_generators.append(list(route))
            if all_basis.rank == dimension-degree_rank and tight_basis.rank == dimension-face_equation_rank:
                return result()
    return result()

if __name__ == '__main__':
    results = [experiment(False), experiment(True)]
    output = Path(__file__).with_name('partition_face_results.json')
    output.write_text(json.dumps(results, indent=2)+'\n')
    print(json.dumps([{k:v for k,v in r.items() if not k.endswith('routes')} for r in results], indent=2))
