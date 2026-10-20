"""Independent batch audit fixtures with 27 synthetic pairs; no native solver."""
import copy
import hashlib
import json
import random
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from experiments.lkh.cached30_batch_audit import audit_batch, digest


def main():
    data = dict(n=5, types=[1, 0, 0, 0, 0], prizes=[1]*5, quota=1,
                radius=0, matrix=[int(i != j) for i in range(5) for j in range(5)], input_hash='toy')
    starts = [dict(instance='synthetic'+str(g), geometry='geometry'+str(g), source_seed=s,
                   cost=5, route=[0, 1, 2, 3, 4], input_sha256='toy')
              for g in range(9) for s in (42, 43, 44)]
    sources = dict(synthetic='source hash')
    checks = []
    with tempfile.TemporaryDirectory(prefix='pcctp-cached30-audit-') as folder:
        out = Path(folder)
        def write(path, record): path.write_text(json.dumps(record, indent=2)+'\n')
        def fixture(costs=None):
            for p in out.iterdir(): p.unlink()
            for i, start in enumerate(starts):
                prefix = f"{start['instance']}.seed{start['source_seed']}"
                order = ['uncached', 'cached']; random.Random(prefix+'cached30').shuffle(order)
                descriptors = {}
                for j, arm in enumerate(order):
                    parent = 1000.+i*10+j*2
                    cost = 4 if arm == 'uncached' else (costs or {}).get(start['geometry'], 3)
                    route = list(range(cost))
                    trace = [] if cost == 5 else [dict(cost=cost, route=route, native_timestamp=parent+.5, native_seconds=.4)]
                    raw = dict(instance=start['instance'], arm='random', initial_cost=5,
                               cost=cost, route=route, native_started=parent+.1, native_elapsed=.8,
                               trace=trace, source_sha256=None, attempts=1, feasible_repairs=1,
                               lkh_calls=1, lkh_completed=1, lkh_timeouts=0, lkh_failed=0)
                    builds = [dict(call_started=parent+.2, call_finished=parent+.3,
                                   build_seconds=.01, allowance_seconds=.2)] if arm == 'cached' else []
                    phase = dict(instance=start['instance'], start=start, repair_seed=start['source_seed']+49,
                                 arm=arm, budget_seconds=30., parent_started=parent, wall_seconds=1.,
                                 raw_phase=raw, cache_build_events=builds,
                                 cache_retained_python_bytes=128 if arm == 'cached' else 0,
                                 cost=cost, route=route, origin='common_start' if cost == 5 else 'phase_event',
                                 source_sha256=sources)
                    path = out/(prefix+'.'+arm+'.json'); write(path, phase)
                    descriptors[arm] = dict(file=path.name, sha256=digest(path), cost=cost)
                environment = out/('environment.'+prefix+'.json')
                write(environment, dict(status='PASS', started_unix=100., finished_unix=103.,
                     max_monitor_gap_seconds=1., max_clock_difference_seconds=.01, sleep_events=[]))
                write(out/(prefix+'.pair.json'), dict(status='PAIR_COMPLETED', start=start,
                     protocol_sha256='protocol', arm_order=order, phases=descriptors,
                     budget_seconds=30., environment_file=environment.name, environment_sha256=digest(environment)))

        def audited(require=True, study=None):
            return audit_batch(out, starts if study is None else study, 'protocol', sources,
                               data_loader=lambda _: data, require_complete=require)

        fixture()
        result = audited()
        assert result['audited_pairs'] == 27 and result['audited_phases'] == 54
        assert result['wins'] == 9 and result['mean_gain_percent'] == 25. and result['gate_passed']
        checks.append('27 pairs/54 phases independently reconstructed, exact mean25%, nine wins and gate pass')
        saved = {p.name: p.read_text() for p in out.iterdir()}
        def restore():
            for p in out.iterdir(): p.unlink()
            for name, text in saved.items(): (out/name).write_text(text)
        def rejected(label, mutation):
            restore(); mutation()
            try: audited()
            except (AssertionError, FileNotFoundError): pass
            else: raise AssertionError('Accepted corruption: '+label)
            checks.append('rejects '+label)
        pairpath = out/'synthetic0.seed42.pair.json'
        def pairchange(mutator):
            pair = json.loads(pairpath.read_text()); mutator(pair); write(pairpath, pair)
        def envchange(mutator):
            pair = json.loads(pairpath.read_text()); path = out/pair['environment_file']
            env = json.loads(path.read_text()); mutator(env); write(path, env)
            pair['environment_sha256'] = digest(path); write(pairpath, pair)
        def phasechange(arm, mutator):
            pair = json.loads(pairpath.read_text()); path = out/pair['phases'][arm]['file']
            phase = json.loads(path.read_text()); mutator(phase); write(path, phase)
            pair['phases'][arm]['sha256'] = digest(path); write(pairpath, pair)
        rejected('wrong protocol', lambda: pairchange(lambda p: p.update(protocol_sha256='changed')))
        rejected('wrong arm order', lambda: pairchange(lambda p: p['arm_order'].reverse()))
        rejected('sleep despite updated file hash', lambda: envchange(lambda e: e.update(sleep_events=['sleep'])))
        rejected('monitor gap beyond five seconds', lambda: envchange(lambda e: e.update(max_monitor_gap_seconds=5.01)))
        rejected('nonfinite environment clock', lambda: envchange(lambda e: e.update(finished_unix=float('nan'))))
        rejected('environment duration shorter than phase time', lambda: envchange(lambda e: e.update(finished_unix=101.)))
        rejected('phase input/source mismatch', lambda: phasechange('cached', lambda p: p.update(source_sha256={})))
        rejected('changed phase cost despite updated hash', lambda: phasechange('cached', lambda p: p.update(cost=2)))
        rejected('unreviewed attempt marker', lambda: (out/'synthetic0.seed42.cached.attempt.json').write_text('{}'))
        rejected('extra result identity', lambda: (out/'extra.pair.json').write_text('{}'))
        rejected('missing scheduled pair', lambda: pairpath.unlink())
        restore(); pairpath.unlink()
        partial = audited(False)
        assert partial['status'] == 'PARTIAL' and partial['audited_pairs'] == 26 and not partial['gate_passed']
        checks.append('partial audit retains denominator54 and cannot pass gate')
        restore(); changed = copy.deepcopy(starts); changed[0]['source_seed'] = 43
        try: audited(study=changed)
        except AssertionError: pass
        else: raise AssertionError('Duplicate start identity accepted')
        checks.append('duplicate start slot rejected')
        fixture({g: 4 for g in {s['geometry'] for s in starts}})
        ties = audited()
        assert ties['ties'] == 9 and ties['mean_gain_percent'] == 0 and not ties['gate_passed']
        checks.append('all ties retained and gate fails')
        fixture({'geometry'+str(g): 4 if g >= 5 else 3 for g in range(9)})
        five = audited()
        assert five['wins'] == 5 and five['mean_gain_percent'] > .2 and not five['gate_passed']
        checks.append('positive mean with only five winning geometries fails prespecified gate')
    files = ['experiments/lkh/'+n for n in ('cached30_batch_audit.py', 'cached30_phase_audit.py',
             'raw_route_audit.py', 'verify_cached30_batch.py')]
    report = dict(status='PASS', scope='Synthetic nine-geometry/27-pair audit and gate fixtures only; '
                  'zero benchmark inputs and native optimizer calls.', checks=checks,
                  source_sha256={p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in files})
    Path(__file__).with_suffix('.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(status='PASS', checks=len(checks), native_optimizer_calls=0)))


if __name__ == '__main__': main()
