"""Exploratory branching oracle, not a solver-performance benchmark.

Python dependencies: numpy, scipy, networkx. Read the original instances without
reductions. LP: degrees, linking, quota, coverage, root GSEC and coverage-forced
cut separation. Prize-forced cuts, blossoms and general rank cuts are NOT included.
Each probed child inherits root cuts and is separated again. All LP calls and
oracle time are counted. No training or test-set selection occurs here.
"""
from __future__ import annotations

import argparse
from itertools import combinations
import json
from pathlib import Path
import time

import networkx as nx
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix


def read_instance(path):
    tokens = Path(path).read_text().split()
    n = int(tokens[0])
    if len(tokens) != 4 + 2*n + n*n:
        raise ValueError(f"Unexpected token count: {path}")
    types = np.array(list(map(int, tokens[4:4+n])))
    prizes = np.array(list(map(float, tokens[4+n:4+2*n])))
    costs = np.array(list(map(float, tokens[4+2*n:]))).reshape(n, n)
    if not np.allclose(costs, costs.T):
        raise ValueError("Undirected costs must be symmetric")
    visits = np.flatnonzero(types != 2)
    mandatory = set(np.flatnonzero(types[visits] == 1).tolist())
    if not mandatory:
        raise ValueError("A mandatory root is required")
    neighborhoods = []
    for w in np.flatnonzero(types == 2):
        neighbors = set(np.flatnonzero(costs[visits, w] <= float(tokens[2]) + 1e-9).tolist())
        if not neighbors:
            raise ValueError(f"Uncoverable demand {w}")
        neighborhoods.append(neighbors)
    return dict(name=Path(path).name, geometry=Path(path).name.split('_VT')[0],
                n=n, visits=visits, mandatory=mandatory,
                prizes=prizes[visits], quota=float(tokens[1]),
                costs=costs[np.ix_(visits, visits)], neighborhoods=neighborhoods)


def sparse_rows(rows, width):
    ii, jj, vv = [], [], []
    for r, row in enumerate(rows):
        for col, value in row.items():
            ii.append(r); jj.append(col); vv.append(value)
    return coo_matrix((vv, (ii, jj)), shape=(len(rows), width)).tocsr()


class Relaxation:
    def __init__(self, instance):
        self.d = instance
        self.n = len(instance['visits'])
        self.edges = list(combinations(range(self.n), 2))
        self.m = len(self.edges)
        self.width = self.m+self.n
        self.root = min(instance['mandatory'])
        self.obj = np.r_[[instance['costs'][i, j] for i, j in self.edges], np.zeros(self.n)]
        self.bounds = [(0., 1.)]*self.m + [
            (1., 1.) if i in instance['mandatory'] else (0., 1.) for i in range(self.n)]
        eq = []
        for v in range(self.n):
            row = {e: 1. for e, edge in enumerate(self.edges) if v in edge}
            row[self.m+v] = -2.
            eq.append(row)
        self.eq = sparse_rows(eq, self.width)
        self.rows = [{self.m+i: -float(p) for i, p in enumerate(instance['prizes']) if p}]
        self.rhs = [-instance['quota']]
        for neighbors in instance['neighborhoods']:
            if not neighbors & instance['mandatory']:
                self.rows.append({self.m+i: -1. for i in neighbors}); self.rhs.append(-1.)
        for e, edge in enumerate(self.edges):
            for v in edge:
                if v not in instance['mandatory']:
                    self.rows.append({e: 1., self.m+v: -1.}); self.rhs.append(0.)
        # Explicitly retain the minimum size of a simple cycle.
        self.rows.append({self.m+i: -1. for i in range(self.n)}); self.rhs.append(-3.)
        self.calls = 0

    def solve(self, initial_cuts=None, fix=None, max_rounds=60):
        cuts = dict(initial_cuts or {})
        bounds = self.bounds.copy()
        if fix is not None:
            bounds[fix[0]] = (float(fix[1]), float(fix[1]))
        began = time.perf_counter()
        for round_id in range(max_rounds):
            rows = self.rows + [row for row, rhs in cuts.values()]
            rhs = self.rhs + [rhs for row, rhs in cuts.values()]
            result = linprog(self.obj, A_ub=sparse_rows(rows, self.width), b_ub=rhs,
                             A_eq=self.eq, b_eq=np.zeros(self.n), bounds=bounds,
                             method='highs', options={'time_limit': 20.})
            self.calls += 1
            if result.status != 0:
                if result.status != 2:
                    raise RuntimeError(f"LP did not finish: {result.message}")
                return result, cuts, time.perf_counter()-began, round_id+1
            graph = nx.Graph()
            graph.add_nodes_from(range(self.n))
            for e, (i, j) in enumerate(self.edges):
                graph.add_edge(i, j, capacity=max(0., float(result.x[e])))
            new = 0

            def add_cut(side, vertex):
                nonlocal new
                side = frozenset(side)
                # Mandatory endpoints all produce the same constant-RHS cut.
                vertex = None if vertex in self.d['mandatory'] else vertex
                key = (tuple(sorted(side)), vertex)
                if key in cuts:
                    return
                row = {e: -1. for e, (i, j) in enumerate(self.edges) if (i in side) != (j in side)}
                rhs_value = -2.
                if vertex is not None:
                    row[self.m+vertex] = 2.
                    rhs_value = 0.
                cuts[key] = (row, rhs_value)
                new += 1

            for v in range(self.n):
                if v == self.root or result.x[self.m+v] < 1e-8:
                    continue
                value, partition = nx.minimum_cut(graph, self.root, v)
                if value < 2*result.x[self.m+v] - 1e-7:
                    add_cut(partition[1], v)
            for neighbors in self.d['neighborhoods']:
                if neighbors & self.d['mandatory']:
                    continue
                augmented = graph.copy()
                sink = self.n
                for v in neighbors:
                    augmented.add_edge(v, sink, capacity=2.)
                value, partition = nx.minimum_cut(augmented, self.root, sink)
                if value < 2. - 1e-7:
                    side = set(partition[1]) - {sink}
                    if not neighbors <= side:
                        raise AssertionError("Invalid coverage contraction")
                    add_cut(side, self.root)
            if not new:
                return result, cuts, time.perf_counter()-began, round_id+1
        raise RuntimeError("Separation round limit reached; no converged result reported")


def probe(path, max_candidates=12):
    began = time.perf_counter()
    d = read_instance(path)
    model = Relaxation(d)
    root, cuts, root_time, rounds = model.solve()
    if root.status != 0:
        raise ValueError("Original instance LP is infeasible")
    frac = [i for i, value in enumerate(root.x) if 1e-6 < value < 1-1e-6]
    optional = [i for i in frac if i >= model.m]
    edge_vars = [i for i in frac if i < model.m]
    # A predeclared balanced shortlist, not exhaustive strong branching.
    order = lambda i: (abs(root.x[i]-.5), i)
    half = max_candidates//2
    chosen = sorted(optional, key=order)[:half] + sorted(edge_vars, key=order)[:half]
    remaining = sorted(set(frac)-set(chosen), key=order)
    chosen += remaining[:max_candidates-len(chosen)]
    # Feasible upper bound: nearest-neighbor cycle visiting ALL visitable vertices.
    remaining_vertices = set(range(model.n)) - {model.root}
    route = [model.root]
    while remaining_vertices:
        nxt = min(remaining_vertices, key=lambda v: (d['costs'][route[-1], v], v))
        route.append(nxt); remaining_vertices.remove(nxt)
    upper = float(sum(d['costs'][route[i-1], route[i]] for i in range(len(route))))
    if sum(d['prizes']) < d['quota']-1e-8 or root.fun > upper+1e-6:
        raise AssertionError("Inconsistent primal/dual bound")
    gain_cap = max(upper-float(root.fun), 1e-6)
    records = []
    for index in chosen:
        gains, child_status, child_time = [], [], 0.
        for value in (0, 1):
            child, _, elapsed, _ = model.solve(cuts, (index, value))
            child_time += elapsed
            child_status.append(int(child.status))
            if child.status == 2:
                gains.append(gain_cap)
            else:
                if child.fun < root.fun-1e-6:
                    raise AssertionError("Child bound weaker despite inherited root cuts")
                gains.append(min(gain_cap, max(0., float(child.fun-root.fun))))
        record = dict(variable=('y' if index >= model.m else 'x'), index=index,
                      value=float(root.x[index]), child_gains=gains,
                      child_status=child_status, probe_seconds=child_time,
                      score=min(gains)+.1*max(gains))
        if index >= model.m:
            v = index-model.m
            relevant = [ns for ns in d['neighborhoods'] if v in ns and not ns & d['mandatory']]
            record.update(vertex=int(d['visits'][v]), prize=float(d['prizes'][v]),
                          residual_demands_covered=len(relevant),
                          coverage_scarcity=sum(1/len(ns) for ns in relevant))
        else:
            record['edge'] = [int(d['visits'][v]) for v in model.edges[index]]
        records.append(record)
    records.sort(key=lambda row: (-row['score'], row['index']))
    return dict(instance=d['name'], geometry=d['geometry'], visitable=model.n,
                mandatory=len(d['mandatory']), optional=model.n-len(d['mandatory']),
                demands=len(d['neighborhoods']), root_lp=float(root.fun),
                all_visit_upper=upper, root_seconds=root_time, root_rounds=rounds,
                root_cuts=len(cuts), fractional_x=len(edge_vars), fractional_y=len(optional),
                probed_candidates=len(records), lp_calls=model.calls,
                total_seconds=time.perf_counter()-began, oracle=records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('instances', nargs='+')
    parser.add_argument('--max-candidates', type=int, default=12)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.max_candidates < 2:
        parser.error('--max-candidates must be >=2')
    result = {'scope': 'root-only diagnostic; no trained policy, no solver benchmark',
              'relaxation': 'degrees/linking/quota/coverage/GSEC/coverage-forced cuts; no blossoms/prize-forced/rank',
              'candidate_selection': 'balanced most-fractional shortlist, then fill; not exhaustive oracle',
              'branch_score': 'min(clipped_child_gain)+0.1*max(clipped_child_gain); infeasible clipped at all-visit UB',
              'instances': []}
    for path in args.instances:
        row = probe(path, args.max_candidates)
        result['instances'].append(row)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2)+'\n')
        print(json.dumps({k: v for k, v in row.items() if k != 'oracle'}), flush=True)


if __name__ == '__main__':
    main()
