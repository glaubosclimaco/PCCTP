"""Matched, conditional TRAIN-only repair screen; frozen protocol in research_notes."""
import hashlib
import json
import random
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from experiments.ml_lns.pilot import read_instance, validate
from experiments.lkh.ruin_repair_absolute import repair as greedy
from experiments.lkh.insertion_mip_repair import repair as joint
from experiments.lkh.screen import polish, BIN
import highspy

OUT = HERE / 'joint_insertion_train'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rng_for(label):
    return random.Random(int(hashlib.sha256(label.encode()).hexdigest(), 16))


def execute(d, original, retained, arm, label):
    began = time.clock_gettime(time.CLOCK_MONOTONIC)
    deadline = began + 1.
    fallback_cost = validate(d, original)
    best, best_cost = original[:], fallback_cost
    events = []
    candidate = None
    details = {}
    rng = rng_for(label + ':' + arm)
    if arm == 'greedy':
        # Existing greedy implementation uses perf_counter for its deadline.
        candidate = greedy(d, retained, rng, time.perf_counter() + max(
            0., deadline - time.clock_gettime(time.CLOCK_MONOTONIC)))
    else:
        slots = {}
        absent = sorted(set(range(len(d['visits']))) - set(retained))
        for v in absent:
            gaps = list(range(len(retained)))
            if arm == 'mip_nearest':
                gaps.sort(key=lambda i: (float(d['costs'][retained[i-1], v]
                    + d['costs'][v, retained[i]]
                    - d['costs'][retained[i-1], retained[i]]), i))
            else:
                rng.shuffle(gaps)
            slots[v] = gaps[:3]
        # Slot selection is charged to the 0.5-second repair sub-budget too.
        remaining = min(began + .5, deadline) - time.clock_gettime(time.CLOCK_MONOTONIC)
        details = dict(status='Slot selection exhausted repair budget', accepted=False)
        if remaining > 0:
            details = joint(d, retained, seconds=remaining, slots=slots)
            if details['accepted']:
                candidate = details['route']
        details = dict(details, slot_variables=sum(map(len, slots.values())),
                       slots=[dict(vertex=int(d['visits'][v]), gaps=gaps)
                              for v, gaps in slots.items()])
    repair_time = time.clock_gettime(time.CLOCK_MONOTONIC)
    feasible_repair = candidate is not None

    def capture(route, phase):
        nonlocal best, best_cost
        cost = validate(d, route)
        timestamp = time.clock_gettime(time.CLOCK_MONOTONIC)
        events.append(dict(phase=phase, native_timestamp=timestamp, cost=cost,
                           route=[int(d['visits'][v]) for v in route]))
        if timestamp <= deadline and cost < best_cost:
            best, best_cost = route[:], cost

    polishing = None
    if candidate is not None:
        capture(candidate, 'repair')
        remaining = deadline - time.clock_gettime(time.CLOCK_MONOTONIC)
        if remaining > 0:
            polishing = polish(d, candidate, seconds=min(.2, remaining))
            mapping = {int(v): i for i, v in enumerate(d['visits'])}
            capture([mapping[v] for v in polishing['route']], 'polish')
    return dict(arm=arm, cost=best_cost,
                route=[int(d['visits'][v]) for v in best],
                fallback_cost=fallback_cost, native_started=began,
                deadline=deadline, repair_returned=repair_time,
                native_elapsed=time.clock_gettime(time.CLOCK_MONOTONIC)-began,
                feasible_repair_returned=feasible_repair,
                events=events, repair_details=details, polishing=polishing)


if __name__ == '__main__':
    # No optimization or reserved instance load before the current batch ends.
    validation = json.loads((HERE/'ruin_validation60_summary.json').read_text())
    assert validation['status'] == 'COMPLETE' and not validation['gate_passed']
    audit = json.loads((HERE/'ruin_validation60_audit.json').read_text())
    assert audit['status'] == 'PASS' and audit['batch_complete']
    training = json.loads((HERE/'ruin_total60_absolute_audit.json').read_text())
    assert training['status'] == 'PASS' and training['batch_complete']
    split_path = ROOT/'experiments/ml_lns/decoder_initialized_states/protocol.json'
    split = json.loads(split_path.read_text())
    sources = sorted((HERE/'ruin_total60_absolute').glob('g_*.hybrid.json'))
    assert len(sources) == 27
    names = sorted({json.loads(p.read_text())['instance'] for p in sources})
    assert len(names) == 9
    for name in names:
        geometry = [g for g in split['roles'] if name.startswith(g+'_VT')]
        assert len(geometry) == 1 and split['roles'][geometry[0]] == 'train'
    paths = [Path(__file__), HERE/'insertion_mip_repair.py',
             HERE/'ruin_repair_absolute.py', HERE/'screen.py',
             ROOT/'experiments/ml_lns/highs_backend.py',
             ROOT/'experiments/ml_lns/pilot.py',
             ROOT/'research_notes/joint_insertion_train_protocol.md']
    protocol = dict(scope='conditional TRAIN only; MDM generation excluded',
        runs=324, source_hashes={str(p.relative_to(ROOT)): digest(p) for p in paths},
        input_hashes={p.name: digest(p) for p in sources},
        split_sha256=digest(split_path), lkh_sha256=digest(BIN),
        highs_version=highspy.Highs().version(),
        validation_summary_sha256=digest(HERE/'ruin_validation60_summary.json'))
    OUT.mkdir(exist_ok=True)
    target = OUT/'protocol.json'
    if target.exists():
        assert json.loads(target.read_text()) == protocol
    else:
        target.write_text(json.dumps(protocol, indent=2)+'\n')
    for source in sources:
        initial = json.loads(source.read_text())
        name, seed = initial['instance'], initial['seed']
        assert seed in (42, 43, 44)
        d = read_instance(ROOT/'data'/name)
        mapping = {int(v): i for i, v in enumerate(d['visits'])}
        original = [mapping[v] for v in initial['pre_route']]
        assert validate(d, original) == initial['pre_cost']
        for sample in range(4):
            label = f'{name}:{seed}:{sample}:20261020'
            rng = rng_for(label)
            optional = sorted(set(original)-d['mandatory'])
            k = min(6, len(optional), len(original)-3)
            removed = set(rng.sample(optional, k))
            retained = [v for v in original if v not in removed]
            arms = ['greedy', 'mip_nearest', 'mip_random']
            rng.shuffle(arms)
            for arm in arms:
                target = OUT/f'{name}.seed{seed}.sample{sample}.{arm}.json'
                if target.exists():
                    continue
                result = execute(d, original, retained, arm, label)
                result.update(instance=name, seed=seed, sample=sample,
                    input_file=source.name, input_sha256=digest(source),
                    removed=[int(d['visits'][v]) for v in sorted(removed)],
                    retained=[int(d['visits'][v]) for v in retained])
                temp = target.with_suffix('.tmp')
                temp.write_text(json.dumps(result, indent=2)+'\n')
                temp.replace(target)
                print(json.dumps({k: result[k] for k in
                                  ('instance', 'seed', 'sample', 'arm', 'cost')}), flush=True)
    (OUT/'completion.json').write_text(json.dumps(dict(runs=324))+'\n')
