"""Exploratory leave-one-geometry-out pairwise ridge ranking; CPU/numpy only.

All candidates for a held-out geometry are excluded from training and scaling.
Evaluates one fixed state per instance, not an end-to-end LNS trajectory. The
oracle is only the best of the sampled, budget-limited neighborhood queries.
"""
import argparse
import json
from pathlib import Path
import time
import numpy as np


def fit_predict(train, test, indices, alpha=1.):
    x = np.array([r['features'] for s in train for r in s['records']])[:, indices]
    mean, scale = x.mean(axis=0), x.std(axis=0)
    scale[scale < 1e-10] = 1.
    xx, yy = [], []
    for state in train:
        records = state['records']
        for i in range(len(records)):
            for j in range(i):
                # No between-instance comparison or leakage of objective scale.
                delta = records[i]['normalized_gain']-records[j]['normalized_gain']
                if abs(delta) < 1e-12:
                    continue
                xx.append((np.array(records[i]['features'])[indices]-
                           np.array(records[j]['features'])[indices])/scale)
                yy.append(delta)
    if not xx:
        weights = np.zeros(len(indices))
    else:
        design = np.array(xx); labels = np.array(yy)
        weights = np.linalg.solve(design.T@design+alpha*np.eye(len(indices)), design.T@labels)
    candidates = (np.array([r['features'] for r in test['records']])[:, indices]-mean)/scale
    return candidates@weights


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    began = time.perf_counter()
    data = json.loads(args.input.read_text())
    states = data['instances']
    if len({s['geometry'] for s in states}) < 3:
        raise ValueError('At least three independent geometries required')
    rows = []
    for test in states:
        train = [s for s in states if s['geometry'] != test['geometry']]
        records = test['records']; gains = np.array([r['normalized_gain'] for r in records])
        if not records:
            continue
        full = fit_predict(train, test, list(range(len(data['feature_names']))))
        # Same linear ranking machinery with only generic route/geometric features.
        generic = fit_predict(train, test, data.get('generic_feature_indices', [0, 6, 7]))
        oracle = float(max(gains))
        selected = {'uniform_candidate_expectation': float(np.mean(gains)),
                    'deletion_savings': float(gains[np.argmax([r['features'][6] for r in records])]),
                    'generic_pairwise_ridge': float(gains[np.argmax(generic)]),
                    'coverage_prize_pairwise_ridge': float(gains[np.argmax(full)])}
        rows.append(dict(instance=test['instance'], geometry=test['geometry'],
                         sampled_oracle_normalized_gain=oracle,
                         normalized_gain=selected,
                         sampled_oracle_regret={k: oracle-v for k, v in selected.items()},
                         train_geometries=sorted({s['geometry'] for s in train})))
    keys = list(rows[0]['normalized_gain'])
    summary = {k: dict(mean_normalized_gain=float(np.mean([r['normalized_gain'][k] for r in rows])),
                       mean_sampled_oracle_regret=float(np.mean([r['sampled_oracle_regret'][k] for r in rows])))
               for k in keys}
    result = dict(scope='Exploratory fixed-state cross-geometry ranking; no significance or runtime superiority claim',
                  ridge_alpha=1., alpha_selected_without_test_tuning=True,
                  geometries=len({r['geometry'] for r in rows}), states=len(rows),
                  informative_states=sum(r['sampled_oracle_normalized_gain'] > 1e-10 for r in rows),
                  fit_and_evaluation_seconds=time.perf_counter()-began, summary=summary, rows=rows)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'rows'}, indent=2))


if __name__ == '__main__':
    main()
