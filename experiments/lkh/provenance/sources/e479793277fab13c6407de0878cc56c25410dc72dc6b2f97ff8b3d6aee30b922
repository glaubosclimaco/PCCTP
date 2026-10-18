"""Fixed-state neighborhood oracle for PCCTP; not a trained LNS algorithm.

Same k optional vertices and same number of free x/y for all candidates.
SCF MILP is integer-exact; its relaxation differs from the paper's branch-and-cut.
Budget labels mean best observed gain, NOT proven neighborhood optimality unless
HiGHS reports optimal. Solver bounds are local and never reported as global.
"""
from __future__ import annotations

import argparse
from itertools import combinations
import json
from pathlib import Path
import random
import sys
import time

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from experiments.ml_branching.oracle_probe import read_instance, sparse_rows


def route_cost(d, route):
    return float(sum(d['costs'][route[i-1], v] for i, v in enumerate(route)))


def validate(d, route):
    selected = set(route)
    if len(route) < 3 or len(selected) != len(route):
        raise ValueError('Not a simple cycle')
    if not selected <= set(range(len(d['visits']))) or not d['mandatory'] <= selected:
        raise ValueError('Missing mandatory or invalid vertex')
    if sum(d['prizes'][list(selected)]) < d['quota']-1e-6:
        raise ValueError('Insufficient prize')
    if not all(selected & ns for ns in d['neighborhoods']):
        raise ValueError('Uncovered demand')
    return route_cost(d, route)


def initial_route(d):
    remaining = set(range(len(d['visits'])))
    route = [min(d['mandatory'])]; remaining.remove(route[0])
    while remaining:
        nxt = min(remaining, key=lambda v: (d['costs'][route[-1], v], v))
        route.append(nxt); remaining.remove(nxt)
    # Deterministic feasible deletion and 2-opt; not the historical GRASP.
    changed = True
    while changed:
        changed = False
        for v in list(route):
            if v in d['mandatory'] or len(route) <= 3:
                continue
            candidate = [u for u in route if u != v]
            try:
                value = validate(d, candidate)
            except ValueError:
                continue
            if value < route_cost(d, route)-1e-8:
                route = candidate; changed = True
        for i in range(len(route)-1):
            for j in range(i+2, len(route)):
                candidate = route[:i+1] + list(reversed(route[i+1:j+1])) + route[j+1:]
                if route_cost(d, candidate) < route_cost(d, route)-1e-8:
                    route = candidate; changed = True
    validate(d, route)
    return route


def free_edge_set(d, route, released, operator='incident', extra_edges=None, seed=42):
    """Matched extra-edge budget for context and randomized expansion controls.

    Context first permits shortcuts across released runs and then frees current
    edges near unvisited candidates to permit insertion. Remaining slots are
    filled by proximity. Every edge choice uses only the current state, no labels.
    """
    if operator == 'fragment_clique':
        if not released: raise ValueError('Fragment clique needs a nonempty release set')
        if extra_edges is not None: raise ValueError('Fragment clique defines its own complete boundary')
        current={tuple(sorted((route[i-1],v))) for i,v in enumerate(route)}
        base={edge for edge in current if released.intersection(edge)}
        if not set(route)&released:
            # Break the retained cycle once, so absent candidates can be inserted.
            gap=min(current,key=lambda e:(min(d['costs'][e[0],v]+d['costs'][v,e[1]]-d['costs'][e[0],e[1]] for v in released),e))
            base.add(gap)
        fixed=current-base
        degree={v:0 for v in route}
        for a,b in fixed:degree[a]+=1;degree[b]+=1
        boundary=set(released)|{v for v in route if degree[v]<2}
        # A retained two-vertex path has both endpoints in the boundary. Its
        # internal edge must stay fixed, rather than being inadvertently freed.
        return set(combinations(sorted(boundary),2))-fixed
    edges = list(combinations(range(len(d['visits'])), 2))
    base = {edge for edge in edges if released.intersection(edge)}
    if operator == 'incident':
        return base
    if operator not in ('random_extra', 'route_context', 'geometric_extra', 'route_random'):
        raise ValueError(f'Unknown operator: {operator}')
    eligible = [edge for edge in edges if edge not in base]
    budget = min(len(eligible), 2*len(released) if extra_edges is None else extra_edges)
    if budget < 0:
        raise ValueError('Negative extra-edge budget')
    if operator == 'random_extra':
        rng = random.Random(seed + sum((i+1)*(v+1)*1009 for i, v in enumerate(sorted(released))))
        return base | set(rng.sample(eligible, budget))
    priorities = {}
    # A released run may disappear only if the edge between its retained
    # endpoints is free. Previously this shortcut was fixed to zero.
    for i, a in enumerate(route):
        if a in released or route[(i+1) % len(route)] not in released:
            continue
        offset = 1
        while offset < len(route) and route[(i+offset) % len(route)] in released:
            offset += 1
        b = route[(i+offset) % len(route)]
        if a != b:
            edge = tuple(sorted((a, b)))
            if edge not in base:
                priorities[edge] = (0, 0., edge)
    absent = released-set(route)
    if absent:
        for i, a in enumerate(route):
            b = route[i-1]; edge = tuple(sorted((a, b)))
            if edge in base or edge in priorities:
                continue
            insertion = min(d['costs'][a, v]+d['costs'][v, b]-d['costs'][a, b] for v in absent)
            priorities[edge] = (1, float(insertion), edge)
    for a, b in eligible:
        edge = (a, b)
        if edge not in priorities:
            proximity = min(d['costs'][a, v] for v in released) + min(d['costs'][b, v] for v in released)
            priorities[edge] = (2, float(proximity), edge)
    if operator == 'geometric_extra':
        proximity = lambda e: (float(min(d['costs'][e[0], v] for v in released) +
                                    min(d['costs'][e[1], v] for v in released)), e)
        return base | set(sorted(eligible, key=proximity)[:budget])
    if operator == 'route_random':
        rng = random.Random(seed + sum((i+1)*(v+1)*1009 for i, v in enumerate(sorted(released))))
        current = {tuple(sorted((route[i-1], v))) for i, v in enumerate(route)}
        relevant = sorted(edge for edge in eligible if edge in current or priorities[edge][0] == 0)
        selected = set(rng.sample(relevant, min(budget, len(relevant))))
        rest = [edge for edge in eligible if edge not in selected]
        selected.update(rng.sample(rest, budget-len(selected)))
        return base | selected
    chosen = sorted(eligible, key=priorities.__getitem__)[:budget]
    return base | set(chosen)


class Repair:
    def __init__(self, d):
        self.d = d; self.n = len(d['visits'])
        self.edges = list(combinations(range(self.n), 2)); self.m = len(self.edges)
        self.f0 = self.m+self.n; self.width = self.f0+2*self.m
        self.obj = np.r_[[d['costs'][i, j] for i, j in self.edges], np.zeros(self.n+2*self.m)]
        self.integrality = np.r_[np.ones(self.f0), np.zeros(2*self.m)]
        self.lb = np.zeros(self.width)
        self.ub = np.r_[np.ones(self.f0), np.full(2*self.m, self.n-1.)]
        for v in d['mandatory']:
            self.lb[self.m+v] = 1.
        rows, lower, upper = [], [], []

        def constraint(row, lo=-np.inf, hi=np.inf):
            rows.append(row); lower.append(lo); upper.append(hi)

        root = min(d['mandatory'])
        for v in range(self.n):
            row = {e: 1. for e, edge in enumerate(self.edges) if v in edge}
            row[self.m+v] = -2.
            constraint(row, 0., 0.)
            # Root supplies one unit to every other selected vertex.
            flow = {}
            for e, (i, j) in enumerate(self.edges):
                a, b = self.f0+2*e, self.f0+2*e+1
                if v == i:
                    flow[a] = -1.; flow[b] = 1.
                if v == j:
                    flow[a] = 1.; flow[b] = -1.
            if v != root:
                flow[self.m+v] = -1.
            else:
                for u in range(self.n):
                    if u != root:
                        flow[self.m+u] = 1.
            constraint(flow, 0., 0.)
        for e, (i, j) in enumerate(self.edges):
            for arc in (self.f0+2*e, self.f0+2*e+1):
                constraint({arc: 1., e: -(self.n-1.)}, hi=0.)
            for v in (i, j):
                constraint({e: 1., self.m+v: -1.}, hi=0.)
        constraint({self.m+i: float(p) for i, p in enumerate(d['prizes'])}, lo=d['quota'])
        constraint({self.m+i: 1. for i in range(self.n)}, lo=3.)
        for ns in d['neighborhoods']:
            constraint({self.m+i: 1. for i in ns}, lo=1.)
        self.constraints = LinearConstraint(sparse_rows(rows, self.width), lower, upper)

    def solve(self, route, released, seconds, operator='incident', extra_edges=None, seed=42):
        start = time.perf_counter()
        selected = set(route)
        incumbent_edges = {tuple(sorted((route[i-1], v))) for i, v in enumerate(route)}
        free_set = free_edge_set(self.d, route, released, operator, extra_edges, seed)
        lb, ub = self.lb.copy(), self.ub.copy()
        free_edges = 0
        for e, edge in enumerate(self.edges):
            if edge in free_set:
                free_edges += 1
            else:
                lb[e] = ub[e] = float(edge in incumbent_edges)
        for v in range(self.n):
            if v not in released:
                lb[self.m+v] = ub[self.m+v] = float(v in selected)
        began_solve = time.perf_counter()
        result = milp(self.obj, integrality=self.integrality,
                      bounds=Bounds(lb, ub), constraints=self.constraints,
                      options={'time_limit': seconds, 'mip_rel_gap': 0., 'disp': False})
        solve_seconds = time.perf_counter()-began_solve
        # The incumbent is mathematically feasible but scipy.milp has no MIP start.
        if result.status not in (0, 1):
            raise RuntimeError(f'Feasible neighborhood reported failure: {result.message}')
        best = route
        returned = result.x is not None
        if returned:
            vector = result.x
            if np.max(np.abs(vector[:self.f0]-np.rint(vector[:self.f0]))) > 1e-5:
                raise ValueError('Nonintegral MILP solution')
            adjacency = {v: [] for v in range(self.n) if vector[self.m+v] > .5}
            for e, (i, j) in enumerate(self.edges):
                if vector[e] > .5:
                    if i not in adjacency or j not in adjacency:
                        raise ValueError('Selected edge without selected endpoints')
                    adjacency[i].append(j); adjacency[j].append(i)
            if any(len(ns) != 2 for ns in adjacency.values()):
                raise ValueError('Invalid degree')
            decoded = [min(self.d['mandatory'])]
            previous = None
            while True:
                nxt = next(v for v in adjacency[decoded[-1]] if v != previous)
                if nxt == decoded[0]:
                    break
                if nxt in decoded:
                    raise ValueError('Repeated vertex before cycle closure')
                previous = decoded[-1]; decoded.append(nxt)
            if set(decoded) != set(adjacency):
                raise ValueError('Disconnected cycles')
            value = validate(self.d, decoded)
            if abs(value-float(result.fun)) > 1e-5:
                raise ValueError('Solver objective mismatch')
            # Check neighborhood fixing independently on decoded edges and vertices.
            decoded_edges = {tuple(sorted((decoded[i-1], v))) for i, v in enumerate(decoded)}
            for edge in self.edges:
                if edge not in free_set and ((edge in decoded_edges) != (edge in incumbent_edges)):
                    raise ValueError('Fixed edge changed')
            for v in range(self.n):
                if v not in released and ((v in decoded) != (v in route)):
                    raise ValueError('Fixed visit changed')
            if value < route_cost(self.d, route):
                best = decoded
        return dict(operator=operator, status=int(result.status), returned_solution=returned,
                    proven_optimal=result.status == 0, solve_seconds=solve_seconds,
                    total_seconds=time.perf_counter()-start, free_x=free_edges,
                    free_y=len(released), cost=validate(self.d, best),
                    removed_optional=[int(self.d['visits'][v]) for v in sorted(selected-set(best))],
                    added_optional=[int(self.d['visits'][v]) for v in sorted(set(best)-selected)],
                    route=[int(self.d['visits'][v]) for v in best])


FEATURES = ['visited_fraction', 'released_prize_fraction', 'residual_prize_fraction',
            'endangered_fraction', 'repairable_fraction', 'coverage_overlap',
            'deletion_savings_fraction', 'normalized_spread']


def features(d, route, released):
    selected = set(route); outside = selected-released
    affected = [ns for ns in d['neighborhoods'] if not ns & outside]
    repairable = sum(bool(ns & (released-selected)) for ns in affected)
    cover_sets = [{w for w, ns in enumerate(d['neighborhoods']) if v in ns} for v in released]
    overlap = [len(a & b)/max(1, len(a | b)) for a, b in combinations(cover_sets, 2)]
    savings = 0.
    for i, v in enumerate(route):
        if v in released:
            prev, nxt = route[i-1], route[(i+1) % len(route)]
            savings += max(0., d['costs'][prev, v]+d['costs'][v, nxt]-d['costs'][prev, nxt])
    spread = [d['costs'][u, v] for u, v in combinations(sorted(released), 2)]
    quota = max(d['quota'], 1.)
    return [len(released & selected)/len(released), float(sum(d['prizes'][list(released)]))/quota,
            max(0., d['quota']-float(sum(d['prizes'][list(outside)])))/quota,
            len(affected)/max(1, len(d['neighborhoods'])), repairable/max(1, len(affected)),
            float(np.mean(overlap)) if overlap else 0., savings/max(route_cost(d, route), 1.),
            float(np.mean(spread))/max(float(np.mean(d['costs'])), 1.) if spread else 0.]


def candidates(d, route, k, count, seed):
    rng = random.Random(seed)
    optional = sorted(set(range(len(d['visits'])))-d['mandatory'])
    k = min(k, len(optional))
    if not k:
        return []
    seen, result = set(), []
    sets = {v: {w for w, ns in enumerate(d['neighborhoods']) if v in ns} for v in optional}
    for attempt in range(count*30):
        kind = ['random', 'coverage', 'prize'][attempt % 3]
        seed_vertex = rng.choice(optional)
        if kind == 'random':
            group = set(rng.sample(optional, k))
        else:
            def score(v):
                overlap = len(sets[v] & sets[seed_vertex])/max(1, len(sets[v] | sets[seed_vertex]))
                if kind == 'coverage':
                    return (overlap, int((v in route) != (seed_vertex in route)), -d['costs'][seed_vertex, v])
                return (float(d['prizes'][v]), overlap, -d['costs'][seed_vertex, v])
            others = sorted((v for v in optional if v != seed_vertex), key=score, reverse=True)
            group = {seed_vertex, *others[:k-1]}
        key = tuple(sorted(group))
        if key in seen:
            continue
        seen.add(key); result.append((kind, group))
        if len(result) == count:
            break
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('instances', nargs='+')
    parser.add_argument('--seconds', type=float, default=.5)
    parser.add_argument('--k', type=int, default=5)
    parser.add_argument('--candidates', type=int, default=12)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if min(args.seconds, args.k, args.candidates) <= 0:
        parser.error('budgets and sizes must be positive')
    output = dict(scope='Fixed-state neighborhood oracle; no trained policy or publication-level evidence',
                  formulation='integer-exact single-commodity flow; no historical GRASP or branch-and-cut',
                  seed=args.seed, seconds_per_candidate=args.seconds, k=args.k,
                  feature_names=FEATURES, instances=[])
    for path in args.instances:
        start = time.perf_counter()
        d = read_instance(path); route = initial_route(d); initial = validate(d, route)
        model = Repair(d); records = []
        for kind, released in candidates(d, route, args.k, args.candidates, args.seed):
            description = features(d, route, released)
            result = model.solve(route, released, args.seconds)
            records.append(dict(kind=kind, released=[int(d['visits'][v]) for v in sorted(released)],
                                features=description, gain=initial-result['cost'],
                                normalized_gain=(initial-result['cost'])/max(initial, 1.), **result))
        if len({(r['free_x'], r['free_y']) for r in records}) > 1:
            raise AssertionError('Neighborhood sizes differ')
        row = dict(instance=d['name'], geometry=d['geometry'], initial_cost=initial,
                   initial_route=[int(d['visits'][v]) for v in route],
                   elapsed_seconds=time.perf_counter()-start, records=records)
        output['instances'].append(row)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(output, indent=2)+'\n')
        print(json.dumps(dict(instance=d['name'], initial=initial,
                              best=min((r['cost'] for r in records), default=initial),
                              queries=len(records), optimal=sum(r['proven_optimal'] for r in records),
                              elapsed_seconds=row['elapsed_seconds'])), flush=True)


if __name__ == '__main__':
    main()
