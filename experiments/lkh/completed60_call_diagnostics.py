"""Passive counters from completed audited hybrids; no optimizer invocation.

The historical counters classify calls to polish(), not necessarily launches
of an LKH subprocess. In particular, the 'failed' bucket also admits a call
whose matrix setup exhausted the allowance before launching the subprocess.
"""
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FIELDS = ('attempts', 'feasible_repairs', 'lkh_calls', 'lkh_completed',
          'lkh_timeouts', 'lkh_failed')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def inspect(prefix, expected_geometries, seeds):
    folder = HERE/prefix
    protocol_path = folder/'protocol.json'
    protocol = json.loads(protocol_path.read_text())
    # Verify that the code semantics being interpreted are those collected.
    recorded = {}
    for section in ('source_hashes', 'dependency_hashes'):
        for name, sha in protocol.get(section, {}).items():
            path = ROOT/name if '/' in name else HERE/name
            recorded[path.resolve()] = sha
    for name in ('ruin_repair_absolute.py', 'screen.py'):
        path = HERE/name
        assert digest(path) == recorded[path.resolve()]
    audit_path = HERE/(prefix+'_audit.json')
    audit = json.loads(audit_path.read_text())
    assert audit['status'] == 'PASS' and audit['batch_complete']
    selected = [row for row in audit['rows'] if row['mode'] == 'hybrid']
    assert len(selected) == expected_geometries*len(seeds)
    rows, hashes, grouped, seen = [], {}, {}, set()
    for row in selected:
        name, seed = row['instance'], row['seed']
        assert seed in seeds and (name, seed) not in seen
        seen.add((name, seed))
        path = folder/f'{name}.seed{seed}.hybrid.json'
        sha = digest(path)
        assert sha == row['result_sha256']
        hashes[str(path.relative_to(ROOT))] = sha
        result = json.loads(path.read_text())
        assert result['instance'] == name and result['seed'] == seed
        post = result['post']
        counters = {key: post[key] if post else 0 for key in FIELDS}
        assert all(type(v) is int and v >= 0 for v in counters.values())
        assert counters['lkh_calls'] == sum(counters[k] for k in
               ('lkh_completed', 'lkh_timeouts', 'lkh_failed'))
        assert counters['attempts'] >= counters['feasible_repairs'] >= counters['lkh_calls']
        if post:
            assert post['instance'] == name and post['arm'] == 'random'
        item = dict(instance=name, seed=seed, post_phase_present=post is not None,
                    **counters)
        rows.append(item)
        grouped.setdefault(name, []).append(item)
    assert len(grouped) == expected_geometries
    geometry = []
    for name, items in sorted(grouped.items()):
        assert {x['seed'] for x in items} == set(seeds)
        total = {key: sum(x[key] for x in items) for key in FIELDS}
        calls = total['lkh_calls']
        geometry.append(dict(instance=name, runs=len(items), **total,
                             mean_polish_invocations=calls/len(items),
                             timeout_percent=100*total['lkh_timeouts']/calls if calls else None,
                             normal_exit_percent=100*total['lkh_completed']/calls if calls else None))
    totals = {key: sum(row[key] for row in rows) for key in FIELDS}
    return dict(runs=len(rows), geometries=expected_geometries, totals=totals,
                geometry=geometry, rows=rows,
                result_sha256=hashes,
                audit_sha256=digest(audit_path), protocol_sha256=digest(protocol_path))


if __name__ == '__main__':
    studies = {label: inspect(prefix, count, seeds) for label, prefix, count, seeds in
               (('TRAIN', 'ruin_total60_absolute', 9, (42, 43, 44)),
                ('Reused validation', 'ruin_validation60', 7, tuple(range(45, 50))))}
    report = dict(status='PASS', scope='Completed 60-second studies only; passive cumulative counters, '
                  'not per-call timing, causal attribution, or successful-tour counts.',
                  category_semantics={
                      'lkh_calls': 'Invocations of polish(), including possible setup-only calls.',
                      'lkh_completed': 'LKH exit code zero; does not prove a timely accepted or improved tour.',
                      'lkh_timeouts': 'Subprocess timeout reported; not an infeasible returned incumbent.',
                      'lkh_failed': 'Other or absent exit code, including possible setup exhaustion without launch.'},
                  limitations=['No per-call logs retained in the aggregate search record.',
                               'Counters can include the last call returning after the outer deadline.',
                               'Do not identify parent process CPU time as total solver CPU time.'],
                  source_sha256={name: digest(HERE/name) for name in
                                 ('completed60_call_diagnostics.py', 'ruin_repair_absolute.py', 'screen.py')},
                  studies=studies)
    (HERE/'completed60_call_diagnostics.json').write_text(json.dumps(report, indent=2)+'\n')
    for label, study in studies.items():
        print(label, json.dumps(study['totals']))
    lines = ['# Passive call diagnostics of completed 60-second studies', '', report['scope'], '',
             'Counters classify polish invocations; normal exit is not a successful or improving tour.', '',
             '| Study / geometry | Invocations | Mean per run | Normal exits | Timeouts | Other/absent exit |',
             '|---|---:|---:|---:|---:|---:|']
    for label, study in studies.items():
        for row in study['geometry']:
            name = row['instance'].split('_VT')[0].removeprefix('g_')
            lines.append(f"| {label}: {name} | {row['lkh_calls']} | {row['mean_polish_invocations']:.1f} | "
                         f"{row['lkh_completed']} | {row['lkh_timeouts']} | {row['lkh_failed']} |")
    lines += ['', 'No optimizer, active 300-second outcome or reserved test input was read.', '',
              'The other/absent-exit bucket cannot distinguish crashes from setup-only calls. '
              'Per-call output validity, improvements and CPU totals cannot be reconstructed '
              'from these counters. No operator or budget is changed in response to this report.']
    (ROOT/'research_notes/completed60_call_diagnostics.md').write_text('\n'.join(lines)+'\n')
