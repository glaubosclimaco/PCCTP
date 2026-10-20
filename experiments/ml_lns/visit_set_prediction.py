"""Predict optimal optional visits on development labels; score completed training runs.

The decision rule is frozen in research_notes/visit_set_prediction.md.
Validation and test geometries are not read.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from experiments.ml_branching.oracle_probe import read_instance
from experiments.ml_lns.pilot import Repair, initial_route

HERE = Path(__file__).resolve().parent
DATA = ROOT / 'data'
PROFILES = ROOT / 'experiments/original_mdm/move_profiles_monotonic_elite_bounded'
ELITES = ROOT / 'experiments/original_mdm/elite_banks'
OPTIMA = {
    'g_att48_VT16_T16_W16_75': 2650,
    'g_berlin52_VT16_T18_W18_50': 5186,
    'g_brazil58_VT11_T12_W35_75': 20466,
    'g_eil51_VT17_T17_W17_25': 273,
    'g_eil76_VT14_T16_W46_75': 332,
    'g_gr48_VT16_T16_W16_25': 3571,
    'g_hk48_VT16_T16_W16_50': 7593,
    'g_st70_VT14_T14_W42_25': 419,
}
EVALUATION = [
    'g_rd100_VT20_T20_W60_75',
    'g_pr144_VT28_T29_W87_50',
    'g_u159_VT31_T32_W96_50',
    'g_kroB150_VT30_T30_W90_25',
    'g_gr137_VT26_T28_W83_25',
    'g_kroC100_VT20_T20_W60_75',
    'g_pr136_VT26_T28_W82_25',
    'g_bier127_VT24_T26_W77_50',
    'g_lin105_VT21_T21_W63_25',
    'g_gr96_VT18_T20_W58_75',
    'g_kroA150_VT30_T30_W90_50',
]
ALPHA = 1.0
THRESHOLD = 0.5
SOLVE_SECONDS = 15.0


def median_scale(d):
    n = len(d['visits'])
    values = [d['costs'][i, j] for i in range(n) for j in range(i)]
    positive = [v for v in values if v > 0]
    return float(np.median(positive if positive else values or [1.]))


def unique_closest(path, d):
    tokens = Path(path).read_text().split()
    n = int(tokens[0])
    types = np.array(list(map(int, tokens[4:4+n])))
    radius = float(tokens[2])
    costs = np.array(list(map(float, tokens[4+2*n:]))).reshape(n, n)
    closest = np.zeros(len(d['visits']))
    for demand in np.flatnonzero(types == 2):
        covered = [(costs[int(d['visits'][i]), demand], i) for i in range(len(d['visits']))
                   if costs[int(d['visits'][i]), demand] <= radius+1e-9]
        best = min(distance for distance, _ in covered)
        winners = [i for distance, i in covered if abs(distance-best) <= 1e-9]
        if len(winners) == 1:
            closest[winners[0]] += 1
    return closest


def vertex_features(d, path):
    scale = median_scale(d)
    prize_max = float(np.max(d['prizes'])) or 1.
    demands = max(1, len(d['neighborhoods']))
    closest = unique_closest(path, d)
    mandatory = sorted(d['mandatory'])
    rows = []
    optional = []
    for v in range(len(d['visits'])):
        if v in d['mandatory']:
            continue
        cover = sum(v in ns for ns in d['neighborhoods']) / demands
        nearest = min(d['costs'][v, u] for u in mandatory) / scale
        mean = float(np.mean([d['costs'][v, u] for u in mandatory])) / scale
        rows.append([
            float(d['prizes'][v]) / prize_max,
            cover,
            closest[v] / demands,
            nearest,
            mean,
        ])
        optional.append(v)
    return optional, np.array(rows, dtype=float), closest


def feasible(d, chosen):
    if not d['mandatory'] <= set(chosen):
        return False
    if float(sum(d['prizes'][v] for v in chosen)) < d['quota'] - 1e-8:
        return False
    return all(set(chosen) & ns for ns in d['neighborhoods'])


def learned_set(d, optional, probability):
    order = [optional[i] for i in np.argsort(-probability)]
    chosen = set(d['mandatory'])
    prefix = 0
    for count, vertex in enumerate(order, start=1):
        if feasible(d, chosen | set(order[:count])):
            chosen |= set(order[:count])
            prefix = count
            break
    else:
        return None, None
    for vertex, p in zip(optional, probability):
        if p >= THRESHOLD:
            chosen.add(vertex)
    return chosen, prefix


def heuristic_set(d, optional, features, size):
    score = features[:, 1] + features[:, 0] - features[:, 3]
    order = [optional[i] for i in np.argsort(-score)]
    chosen = set(d['mandatory']) | set(order[:size])
    if len(chosen - d['mandatory']) != size or not feasible(d, chosen):
        return None
    return chosen


def fit(x, y):
    design = np.c_[np.ones(len(x)), x]
    weights = np.zeros(design.shape[1])
    positive = float(y.sum())
    sample = np.where(y > 0.5, len(y) / (2 * positive), len(y) / (2 * (len(y) - positive)))
    penalty = ALPHA * np.diag(np.r_[0., np.ones(design.shape[1] - 1)])
    for _ in range(40):
        linear = np.clip(design @ weights, -30, 30)
        prob = 1 / (1 + np.exp(-linear))
        residual = sample * (prob - y)
        curvature = sample * prob * (1 - prob)
        hessian = design.T @ (design * curvature[:, None]) + penalty
        gradient = design.T @ residual + penalty @ weights
        step = np.linalg.solve(hessian, gradient)
        weights -= step
        if np.max(np.abs(step)) < 1e-8:
            break
    return weights


def predict(weights, x):
    linear = np.clip(np.c_[np.ones(len(x)), x] @ weights, -30, 30)
    return 1 / (1 + np.exp(-linear))


def standardize(train, test):
    mean = train.mean(axis=0)
    scale = train.std(axis=0)
    scale[scale < 1e-10] = 1.
    return (train - mean) / scale, (test - mean) / scale, mean, scale


def subset(d, chosen):
    chosen = sorted(chosen)
    index = {v: i for i, v in enumerate(chosen)}
    neighborhoods = []
    for ns in d['neighborhoods']:
        mapped = {index[v] for v in ns if v in index}
        if not mapped:
            return None
        neighborhoods.append(mapped)
    return dict(
        name=d['name'], geometry=d['geometry'], n=d['n'],
        visits=np.array([d['visits'][v] for v in chosen]),
        mandatory={index[v] for v in d['mandatory']},
        prizes=d['prizes'][chosen], quota=d['quota'],
        costs=d['costs'][np.ix_(chosen, chosen)], neighborhoods=neighborhoods,
    )


def solve_subset(d, chosen):
    reduced = subset(d, chosen)
    if reduced is None:
        return None
    route = initial_route(reduced)
    repair = Repair(reduced)
    result = repair.solve(route, set(range(repair.n)), SOLVE_SECONDS, operator='incident')
    return result


def labels_for(name):
    d = read_instance(DATA / name)
    elite = json.loads((ELITES / f'{name}.seed42.json').read_text())
    if int(elite['cost']) != OPTIMA[name] or not elite['completed']:
        raise ValueError(f'Elite route is not the proven optimum: {name}')
    mapping = {int(v): i for i, v in enumerate(d['visits'])}
    route = {mapping[int(v)] for v in elite['route']}
    optional, features, _ = vertex_features(d, DATA / name)
    y = np.array([1. if v in route else 0. for v in optional])
    return d, optional, features, y


def main():
    began = time.perf_counter()
    train_x, train_y = [], []
    label_rows = []
    for name in OPTIMA:
        d, optional, features, y = labels_for(name)
        train_x.append(features)
        train_y.append(y)
        label_rows.append(dict(instance=name, optional=len(optional), positives=int(y.sum())))
    train_x = np.vstack(train_x)
    train_y = np.concatenate(train_y)
    if train_y.sum() == 0 or train_y.sum() == len(train_y):
        raise ValueError('Labels have a single class')
    scaled, _, mean, scale = standardize(train_x, train_x)
    weights = fit(scaled, train_y)
    rows = []
    for name in EVALUATION:
        profile = json.loads((PROFILES / f'{name}.seed42.profiled_baseline.json').read_text())
        if not profile['completed'] or profile['iterations'] != 70:
            rows.append(dict(instance=name, status='excluded_incomplete_baseline'))
            continue
        d = read_instance(DATA / name)
        optional, features, _ = vertex_features(d, DATA / name)
        scaled_test = (features - mean) / scale
        probability = predict(weights, scaled_test)
        chosen, prefix = learned_set(d, optional, probability)
        mdm = float(profile['cost'])
        row = dict(instance=name, geometry=d['geometry'], mdm_cost=mdm,
                   optional=len(optional), prefix=prefix)
        if chosen is None:
            row.update(status='learned_infeasible', learned_better_than_mdm=False)
            rows.append(row)
            continue
        heuristic = heuristic_set(d, optional, features, len(chosen - d['mandatory']))
        try:
            learned = solve_subset(d, chosen)
        except Exception as error:
            learned = None
            row['learned_error'] = str(error)
        row.update(learned_optional=len(chosen - d['mandatory']),
                   learned_probability_min=float(min(probability[optional.index(v)] for v in chosen - d['mandatory'])) if chosen - d['mandatory'] else None)
        if learned is None:
            row.update(status='learned_unsolved', learned_better_than_mdm=False)
        else:
            row.update(status='solved', learned_cost=learned['cost'], learned_proven=learned['proven_optimal'],
                       learned_seconds=learned['solve_seconds'],
                       learned_better_than_mdm=learned['cost'] < mdm - 1e-6)
        if heuristic is None:
            row.update(heuristic_status='infeasible')
        else:
            try:
                heuristic_result = solve_subset(d, heuristic)
            except Exception as error:
                heuristic_result = None
                row['heuristic_error'] = str(error)
            if heuristic_result is None:
                row.update(heuristic_status='unsolved')
            else:
                row.update(heuristic_status='solved', heuristic_cost=heuristic_result['cost'],
                           heuristic_proven=heuristic_result['proven_optimal'])
                if learned is not None:
                    row['learned_better_than_heuristic'] = learned['cost'] < heuristic_result['cost'] - 1e-6
                    row['learned_worse_than_heuristic'] = learned['cost'] > heuristic_result['cost'] + 1e-6
        rows.append(row)
        print(json.dumps({k: row[k] for k in row if k != 'route'}), flush=True)
    comparable = [r for r in rows if r.get('status') == 'solved']
    learned_wins = sum(bool(r.get('learned_better_than_mdm')) for r in rows if 'learned_better_than_mdm' in r)
    heuristic_beats = sum(bool(r.get('learned_better_than_heuristic')) for r in rows)
    heuristic_losses = sum(bool(r.get('learned_worse_than_heuristic')) for r in rows)
    improves = learned_wins >= 3 and heuristic_beats > heuristic_losses
    source = Path(__file__).read_bytes()
    result = dict(
        scope='Application of L2 logistic visit prediction; rule frozen before training costs were compared',
        alpha=ALPHA, threshold=THRESHOLD, solve_seconds=SOLVE_SECONDS,
        code_sha256=hashlib.sha256(source).hexdigest(),
        label_rows=label_rows, label_positives=int(train_y.sum()), label_rows_total=int(len(train_y)),
        weights=weights.tolist(),
        learned_strictly_better_than_mdm=learned_wins,
        learned_better_than_heuristic=heuristic_beats,
        learned_worse_than_heuristic=heuristic_losses,
        solved=len(comparable),
        improves=improves,
        seconds=time.perf_counter() - began,
        rows=rows,
    )
    target = HERE / 'visit_set_prediction.json'
    target.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('rows', 'weights', 'label_rows')}, indent=2))


if __name__ == '__main__':
    main()
