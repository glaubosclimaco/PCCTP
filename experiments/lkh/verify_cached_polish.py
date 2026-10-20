"""Synthetic deadline and output-contract checks; no native optimizer launched."""
import hashlib
import json
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from experiments.lkh import cached_polish as module
from experiments.lkh import screen

if __name__ == '__main__':
    costs = np.array([[0,10,1,1], [10,0,1,1], [1,1,0,10], [1,1,10,0]], dtype=float)
    data = dict(visits=np.array([100,101,102,103]), costs=costs,
                mandatory={0}, prizes=np.ones(4), quota=1, neighborhoods=[])
    clock = {'now': 100., 'build': .03, 'action': 'valid', 'launches': 0}
    original = module.time, module.MatrixTextCache, module.subprocess.run
    written_file_sets = []
    def build(matrix):
        clock['now'] += clock['build']
        return original[1](matrix)
    def native(*args, **kwargs):
        clock['launches'] += 1
        assert kwargs['timeout'] > 0
        folder = Path(kwargs['cwd'])
        written_file_sets.append({name: (folder/name).read_text()
                                  for name in ('instance.tsp', 'initial.tour', 'run.par')})
        text = (folder/'instance.tsp').read_text()
        dimension = int(next(line.split(':')[1] for line in text.splitlines() if line.startswith('DIMENSION:')))
        action = clock['action']
        order = [1,3,2,4] if dimension == 4 else list(range(1,dimension+1))
        if action == 'invalid': order[1] = order[0]
        if action == 'unchanged': order = list(range(1,dimension+1))
        if action != 'missing':
            (folder/'result.tour').write_text('TOUR_SECTION\n'+'\n'.join(map(str,order))+'\n-1\nEOF\n')
        clock['now'] += .21 if action == 'late' else .02
        if action == 'timeout': raise subprocess.TimeoutExpired(args[0], kwargs['timeout'])
        return SimpleNamespace(returncode=1 if action == 'nonzero' else 0, stdout='synthetic')
    checks = []
    try:
        module.time = SimpleNamespace(perf_counter=lambda: clock['now'])
        module.MatrixTextCache = build
        module.subprocess.run = native
        polisher = module.CachedPolisher()
        clock['build'] = .3
        spent = polisher(data, [0,1,2,3], seconds=.2)
        assert clock['launches'] == 0 and spent['cost'] == 22 and abs(spent['seconds']-.3) < 1e-10
        assert spent['matrix_cache_created'] and abs(spent['matrix_cache_build_seconds']-.3) < 1e-10
        checks.append('cold cache construction charged; exhausted allowance blocks launch')
        good = polisher(data, [0,1,2,3], seconds=.2)
        assert good['cost'] == 4 and good['improved'] and not good['matrix_cache_created']
        assert good['route'] == [100,102,101,103] and polisher.builds == 1
        changed = polisher(data, [0,1,2], seconds=.2)
        assert changed['route'] == [100,101,102] and polisher.builds == 1
        checks += ['cache reused with same matrix and changed supports',
                   'timely improving tour accepted with original location IDs']
        try:
            polisher(dict(data, costs=costs.copy()), [0,1,2,3], seconds=.2)
        except AssertionError: pass
        else: raise AssertionError('Different instance matrix silently reused')
        checks.append('different matrix object rejected')
        clock['build'] = .03
        for action in ('late', 'timeout', 'nonzero', 'missing', 'unchanged'):
            clock['action'] = action
            result = module.CachedPolisher()(data, [0,1,2,3], seconds=.2)
            assert result['cost'] == 22 and result['route'] == [100,101,102,103]
            assert not result['improved']
            if action == 'timeout': assert result['timeout']
            checks.append(action+' output preserves original feasible route')
        clock['action'] = 'invalid'
        try:
            module.CachedPolisher()(data, [0,1,2,3], seconds=.2)
        except AssertionError: pass
        else: raise AssertionError('Duplicate tour entries accepted')
        checks.append('malformed tour rejected')
        before = clock['launches']
        clock['action'] = 'valid'
        result = module.CachedPolisher()(data, [0,1,2,3], began=clock['now']-.25, seconds=.2)
        assert result['cost'] == 22 and clock['launches'] == before
        checks.append('external start-time debit retained')
        before = len(written_file_sets)
        screen.polish(data, [0,1,2,3], seconds=.2)
        module.CachedPolisher()(data, [0,1,2,3], seconds=.2)
        assert len(written_file_sets) == before+2
        assert written_file_sets[-2] == written_file_sets[-1]
        checks.append('TSP, initial-tour and parameter files match the actual original writer')
    finally:
        module.time, module.MatrixTextCache, module.subprocess.run = original
    files = [Path(__file__), Path(module.__file__), Path(module.__file__).with_name('matrix_text_cache.py'),
             Path(screen.__file__), ROOT/'experiments/ml_lns/pilot.py']
    report = dict(status='PASS', scope='Synthetic clocks, tours and subprocess stubs; zero optimizer calls or benchmark inputs.',
                  checks=checks, source_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
    Path(__file__).with_suffix('.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report))
