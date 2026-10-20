"""Finite, independent repair checks; no benchmark data or frozen-search changes."""
import hashlib
import json
import random
from itertools import product
from pathlib import Path

import numpy as np
from insertion_mip_repair import repair


def enumerate_assignments(d, base, slots):
    absent = sorted(set(range(len(d['visits']))) - set(base))
    best = None
    checked = 0
    for choices in product(range(-1, len(base)), repeat=len(absent)):
        checked += 1
        used = [i for i in choices if i >= 0]
        if len(used) != len(set(used)):
            continue
        if slots is not None and any(i >= 0 and i not in slots.get(v, [])
                                     for v, i in zip(absent, choices)):
            continue
        at = {i: v for v, i in zip(absent, choices) if i >= 0}
        route = []
        for i, v in enumerate(base):
            if i in at:
                route.append(at[i])
            route.append(v)
        selected = set(route)
        # Deliberately independent of pilot.validate and the MIP row builder.
        if len(selected) != len(route) or not d['mandatory'] <= selected:
            continue
        if sum(float(d['prizes'][v]) for v in route) < d['quota']:
            continue
        if any(not selected.intersection(c) for c in d['neighborhoods']):
            continue
        value = sum(float(d['costs'][route[i-1], v]) for i, v in enumerate(route))
        best = value if best is None else min(best, value)
    return best, checked


if __name__ == '__main__':
    rng = random.Random(20261020)
    rows = []
    total = 0
    for trial in range(60):
        b = 3 + trial % 3
        k = trial % 4
        n = b + k
        base = list(range(b))
        rng.shuffle(base)
        costs = np.zeros((n, n))
        for i in range(n):
            for j in range(i):
                costs[i, j] = costs[j, i] = rng.randrange(0, 100)
        prizes = np.array([rng.randrange(0, 20) for _ in range(n)])
        quota = rng.randrange(0, int(sum(prizes)) + 2)
        covers = [{v for v in range(n) if rng.random() < .4}
                  for _ in range(trial % 4)]
        d = dict(visits=np.arange(n), mandatory={base[0]}, costs=costs,
                 prizes=prizes, quota=quota, neighborhoods=covers)
        slots = None if trial % 3 == 0 else {
            v: [i for i in range(b) if rng.random() < .5] for v in range(b, n)}
        if trial % 10 == 1:
            slots = {}  # Missing candidates means no insertion options.
        optimum, count = enumerate_assignments(d, base, slots)
        total += count
        result = repair(d, base, seconds=2., slots=slots)
        if optimum is None:
            assert result['route'] is None and not result['accepted'], result
            assert 'Infeasible' in result['status'], result
        else:
            assert result['accepted'] and result['proven_local_optimal'], result
            assert abs(result['cost'] - optimum) < 1e-8, (result, optimum)
        rows.append(dict(trial=trial, retained=b, absent=k,
                         restricted=slots is not None, optimum=optimum,
                         solver_cost=result['cost']))

    # Feasible identity repair after the deadline must not be eligible either.
    d = dict(visits=np.arange(3), mandatory={0}, costs=np.ones((3, 3)),
             prizes=np.ones(3), quota=1, neighborhoods=[{0}])
    expired = repair(d, [0, 1, 2], seconds=0.)
    assert expired['route'] == [0, 1, 2] and not expired['accepted']
    source = Path(__file__).with_name('insertion_mip_repair.py')
    report = dict(status='PASS', models=len(rows), assignments_checked=total,
                  expired_identity_rejected=True,
                  source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                  verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  scope='finite independent enumeration only; no benchmark performance claim',
                  rows=rows)
    Path(__file__).with_suffix('.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k != 'rows'}))
