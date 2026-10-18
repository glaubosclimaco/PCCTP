"""Paired operator study with matched extra x-variable budget, no learned selector.

Frozen development geometries from geometry_pilot.json. Default seeds 42,43,44;
all operators evaluate exactly the same candidates from the same incumbent.
Operator execution order is shuffled independently of objective outcomes.
"""
import argparse
from collections import defaultdict
import json
from pathlib import Path
import platform
import random
import time

import numpy as np
import scipy
from pilot import (FEATURES, Repair, candidates, features, free_edge_set,
                   initial_route, read_instance, route_cost, validate)


ROUTE_FEATURES = ['available_shortcut_savings', 'released_runs_fraction',
                  'minimum_open_insertion_cost', 'open_incumbent_edges_fraction']


def route_features(d, route, released, free):
    cost = max(1., route_cost(d, route))
    savings = 0.; runs = 0
    for i, v in enumerate(route):
        prev, nxt = route[i-1], route[(i+1) % len(route)]
        if v in released:
            if prev not in released:
                runs += 1
            if tuple(sorted((prev, nxt))) in free:
                savings += max(0., d['costs'][prev, v]+d['costs'][v, nxt]-d['costs'][prev, nxt])
    opened = [(route[i-1], v) for i, v in enumerate(route)
              if tuple(sorted((route[i-1], v))) in free]
    absent = released-set(route)
    insertion = [d['costs'][a, v]+d['costs'][v, b]-d['costs'][a, b]
                 for a, b in opened for v in absent if v not in (a, b)]
    return [savings/cost, runs/len(released),
            float(min(insertion))/cost if insertion else 0., len(opened)/len(route)]


def summarize(states):
    # Average across seeds and candidates WITHIN geometry before across geometries.
    geometry = defaultdict(lambda: defaultdict(list))
    paired = dict(context_better=0, random_better=0, ties=0, pairs=0, all_proven_pairs=0)
    optimal = defaultdict(int); queries = defaultdict(int); substitutions = defaultdict(int)
    for state in states:
        for record in state['records']:
            for op, result in record['results'].items():
                geometry[state['geometry']][op].append(result['normalized_gain'])
                queries[op] += 1; optimal[op] += result['proven_optimal']
                substitutions[op] += bool(result['removed_optional'] or result['added_optional'])
            ctx, rnd = record['results']['route_context'], record['results']['random_extra']
            difference = rnd['cost']-ctx['cost']
            key = 'context_better' if difference > 1e-6 else 'random_better' if difference < -1e-6 else 'ties'
            paired[key] += 1; paired['pairs'] += 1
            paired['all_proven_pairs'] += ctx['proven_optimal'] and rnd['proven_optimal']
            for op in ('random_extra', 'route_context'):
                base = record['results']['incident']; expanded = record['results'][op]
                if base['proven_optimal'] and expanded['proven_optimal'] and expanded['cost'] > base['cost']+1e-6:
                    raise AssertionError('Expansion cannot worsen an optimally solved nested neighborhood')
    controls = {}
    for op in queries:
        if op in ('incident', 'route_context'):
            continue
        stats = dict(context_better=0, control_better=0, ties=0, pairs=0, all_proven_pairs=0)
        for state in states:
            for record in state['records']:
                ctx, other = record['results']['route_context'], record['results'][op]
                diff = other['cost']-ctx['cost']
                stats['context_better' if diff > 1e-6 else 'control_better' if diff < -1e-6 else 'ties'] += 1
                stats['pairs'] += 1
                stats['all_proven_pairs'] += ctx['proven_optimal'] and other['proven_optimal']
        controls[op] = stats
    return dict(geometries=len(geometry), states=len(states), queries=dict(queries),
                proven_optimal=dict(optimal), changed_optional_selection=dict(substitutions),
                matched_pairs=paired, matched_controls=controls,
                mean_normalized_gain={op: float(np.mean([np.mean(v[op]) for v in geometry.values()]))
                                      for op in queries})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path(__file__).with_name('geometry_pilot.json'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--ranker-input', type=Path, required=True)
    parser.add_argument('--seeds', nargs='+', type=int, default=[42, 43, 44])
    parser.add_argument('--seconds', type=float, default=.5)
    parser.add_argument('--k', type=int, default=5)
    parser.add_argument('--candidates', type=int, default=12)
    parser.add_argument('--strong-controls', action='store_true')
    args = parser.parse_args()
    if min(args.seconds, args.k, args.candidates) <= 0:
        parser.error('budgets and sizes must be positive')
    source = json.loads(args.source.read_text())
    data_dir = Path(__file__).resolve().parents[2]/'data'
    output = dict(scope='Exploratory paired operator study on development geometries; not confirmatory',
                  selection='Same frozen geometry set as source; all results retained',
                  source=str(args.source), seeds=args.seeds, k=args.k,
                  strong_controls=args.strong_controls,
                  seconds_per_query=args.seconds, extra_edge_budget='min(2*k, available_nonincident_edges)',
                  versions=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__),
                  states=[])
    ranker = dict(scope='Fixed-state cross-geometry exploratory ranking; same candidate pool',
                  feature_names=FEATURES+ROUTE_FEATURES,
                  generic_feature_indices=[0, 6, 7, 8, 9, 10, 11], instances=[])
    began = time.perf_counter()
    for old in source['instances']:
        d = read_instance(data_dir/old['instance']); route = initial_route(d)
        initial = validate(d, route); model = Repair(d)
        if abs(initial-old['initial_cost']) > 1e-6:
            raise AssertionError('Initial state changed from frozen pilot')
        for seed in args.seeds:
            state = dict(instance=d['name'], geometry=d['geometry'], seed=seed,
                         initial_cost=initial, records=[])
            ranking_state = dict(instance=d['name'], geometry=d['geometry'], seed=seed, records=[])
            for candidate_index, (kind, released) in enumerate(candidates(d, route, args.k, args.candidates, seed)):
                ops = ['incident', 'random_extra', 'route_context']
                if args.strong_controls:
                    ops += ['geometric_extra', 'route_random']
                random.Random(seed+candidate_index*7919).shuffle(ops)
                row = dict(kind=kind, released=[int(d['visits'][v]) for v in sorted(released)], results={})
                for op in ops:
                    result = model.solve(route, released, args.seconds, operator=op, seed=seed)
                    result['normalized_gain'] = (initial-result['cost'])/max(initial, 1.)
                    row['results'][op] = result
                ctx, rnd = row['results']['route_context'], row['results']['random_extra']
                if (ctx['free_x'], ctx['free_y']) != (rnd['free_x'], rnd['free_y']):
                    raise AssertionError('Matched-budget comparison violated')
                for op in ops:
                    r = row['results'][op]
                    if op != 'incident' and (r['free_x'], r['free_y']) != (ctx['free_x'], ctx['free_y']):
                        raise AssertionError('Strong control budget differs')
                free = free_edge_set(d, route, released, 'route_context', seed=seed)
                ranking_state['records'].append(dict(kind=kind,
                    features=features(d, route, released)+route_features(d, route, released, free),
                    normalized_gain=ctx['normalized_gain']))
                state['records'].append(row)
            output['states'].append(state); ranker['instances'].append(ranking_state)
            output['summary'] = summarize(output['states'])
            output['elapsed_seconds'] = time.perf_counter()-began
            args.output.write_text(json.dumps(output, indent=2)+'\n')
            args.ranker_input.write_text(json.dumps(ranker, indent=2)+'\n')
            print(json.dumps(dict(instance=d['name'], seed=seed,
                                  elapsed=output['elapsed_seconds'], summary=output['summary'])), flush=True)


if __name__ == '__main__':
    main()
