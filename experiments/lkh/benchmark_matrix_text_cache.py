"""Serialization microbenchmark on nine previously used TRAIN starts only.

The real polish writer is exercised with subprocess.run replaced by a capture
stub: no LKH, MIP or other optimizer is launched. This is not an end-to-end
algorithm experiment and cannot establish a solution-quality advantage.
"""
import hashlib
import json
import math
import random
import statistics
import sys
import time
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from experiments.lkh import screen
from experiments.lkh.matrix_text_cache import MatrixTextCache
from experiments.ml_lns.pilot import read_instance


def original_assembly(costs, support):
    chosen = sorted(support)
    matrix = ['NAME: PCCTP_SUPPORT', 'TYPE: TSP', f'DIMENSION: {len(chosen)}',
              'EDGE_WEIGHT_TYPE: EXPLICIT', 'EDGE_WEIGHT_FORMAT: FULL_MATRIX',
              'EDGE_WEIGHT_SECTION']
    for a in chosen:
        values = [costs[a, b] for b in chosen]
        assert all(v == int(v) for v in values)
        matrix.append(' '.join(str(int(v)) for v in values))
    return '\n'.join(matrix+['EOF'])+'\n'


if __name__ == '__main__':
    profile_path = HERE/'hard_train_instance_profile.json'
    profile = json.loads(profile_path.read_text())
    roles_path = ROOT/'experiments/ml_lns/decoder_initialized_states/protocol.json'
    roles = json.loads(roles_path.read_text())['roles']
    assert len(profile['rows']) == 9
    records, hashes = [], {}
    for item in profile['rows']:
        name = item['instance']
        assert roles[item['geometry']] == 'train'
        source = HERE/'ruin_total60_absolute'/f'{name}.seed42.hybrid.json'
        result = json.loads(source.read_text())
        hashes[str(source.relative_to(ROOT))] = hashlib.sha256(source.read_bytes()).hexdigest()
        data_path = ROOT/'data'/name
        assert hashlib.sha256(data_path.read_bytes()).hexdigest() == item['input_sha256']
        d = read_instance(data_path)
        mapping = {int(v): i for i, v in enumerate(d['visits'])}
        route = [mapping[v] for v in result['pre_route']]
        captured = []
        original_run = screen.subprocess.run
        def capture(*args, **kwargs):
            captured.append((Path(kwargs['cwd'])/'instance.tsp').read_bytes())
            return SimpleNamespace(returncode=1, stdout='capture stub; no solver launched')
        try:
            screen.subprocess.run = capture
            probe = screen.polish(d, route, seconds=10.)
        finally:
            screen.subprocess.run = original_run
        assert len(captured) == 1 and probe['route'] == result['pre_route']
        began = time.perf_counter()
        cache = MatrixTextCache(d['costs'])
        setup = time.perf_counter()-began
        old = original_assembly(d['costs'], route)
        assert old.encode() == captured[0] == cache.serialize(route).encode()
        samples = {'original': [], 'cached': []}
        order_rng = random.Random(name+'serialize-only')
        for _ in range(5):
            order = ['original', 'cached']; order_rng.shuffle(order)
            for method in order:
                began = time.perf_counter()
                value = original_assembly(d['costs'], route) if method == 'original' else cache.serialize(route)
                elapsed = time.perf_counter()-began
                assert value == old
                samples[method].append(elapsed)
        baseline = statistics.median(samples['original'])
        cached = statistics.median(samples['cached'])
        saving = baseline-cached
        records.append(dict(instance=name, visitable_locations=len(d['visits']),
                            support_size=len(route), output_bytes=len(captured[0]),
                            output_sha256=hashlib.sha256(captured[0]).hexdigest(),
                            exact_actual_writer_match=True, cache_setup_seconds=setup,
                            original_median_seconds=baseline, cached_median_seconds=cached,
                            serialization_speed_ratio=baseline/cached,
                            estimated_reuse_break_even_calls=math.ceil(setup/saving) if saving > 0 else None,
                            samples_seconds=samples))
    source_files = [Path(__file__), HERE/'matrix_text_cache.py', HERE/'screen.py', profile_path, roles_path]
    report = dict(status='PASS', scope='Serialization-only microbenchmark; no optimizer, active 300-second '
                  'outcome or reserved test input. Five interleaved timing repetitions per method on '
                  'nine fixed seed42 TRAIN supports. Not a controlled end-to-end speedup.',
                  limits=['Cache setup measured separately and must be charged in future search.',
                          'Cache retains O(number of visitable locations squared) string tokens.',
                          'Full matrix must contain finite integral values.',
                          'Disk I/O, route validation, process launch and LKH search excluded from timings.',
                          'Host also running the ongoing study; no cross-machine or search-speed claim.'],
                  source_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                                 for p in source_files}, start_record_sha256=hashes, rows=records)
    (HERE/'matrix_text_cache_benchmark.json').write_text(json.dumps(report, indent=2)+'\n')
    for row in records:
        print(json.dumps({key: value for key, value in row.items() if key != 'samples_seconds'}))
