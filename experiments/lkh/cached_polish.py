"""Isolated, unused ordering adapter with charged lazy matrix-text caching.

Make one CachedPolisher per search phase. Problem costs must remain unchanged;
object identity rejects reuse across instances but does not detect in-place edits.
No active experiment or prospective paired70 candidate imports this module.
"""
import subprocess
import tempfile
import time
from pathlib import Path

from experiments.lkh.matrix_text_cache import MatrixTextCache
from experiments.lkh.screen import BIN
from experiments.ml_lns.pilot import validate


class CachedPolisher:
    def __init__(self):
        self.cache = None
        self.costs = None
        self.builds = 0

    def __call__(self, d, route, began=None, seconds=1.):
        # The first call's clock starts BEFORE cache construction.
        began = time.perf_counter() if began is None else began
        old = validate(d, route)
        created = self.cache is None
        build_seconds = 0.
        if created:
            setup_began = time.perf_counter()
            self.cache = MatrixTextCache(d['costs'])
            build_seconds = time.perf_counter()-setup_began
            self.costs = d['costs']
            self.builds += 1
        assert d['costs'] is self.costs, 'Use a new cache for each unchanged instance matrix'
        chosen = sorted(route)
        index = {v: i+1 for i, v in enumerate(chosen)}
        with tempfile.TemporaryDirectory(dir='/private/tmp', prefix='pcctp-cached-lkh-') as tmp:
            folder = Path(tmp)
            (folder/'instance.tsp').write_text(self.cache.serialize(chosen))
            (folder/'initial.tour').write_text('\n'.join(
                ['NAME: INITIAL', 'TYPE: TOUR', f'DIMENSION: {len(chosen)}', 'TOUR_SECTION']
                + [str(index[v]) for v in route] + ['-1', 'EOF'])+'\n')
            (folder/'run.par').write_text(
                'PROBLEM_FILE = instance.tsp\nINITIAL_TOUR_FILE = initial.tour\n'
                'TOUR_FILE = result.tour\nRUNS = 1\nMAX_TRIALS = 100\nSEED = 42\n'
                'TIME_LIMIT = '+str(.8*seconds)+'\nTRACE_LEVEL = 0\n')
            timeout, code, log, best = False, None, '', route
            remaining = seconds-(time.perf_counter()-began)
            if remaining > 0:
                try:
                    result = subprocess.run([str(BIN), 'run.par'], cwd=folder,
                                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                            text=True, timeout=remaining)
                    code, log = result.returncode, result.stdout
                except subprocess.TimeoutExpired:
                    timeout = True
            result_path = folder/'result.tour'
            if code == 0 and result_path.exists() and time.perf_counter()-began <= seconds:
                tokens = result_path.read_text().split()
                start = tokens.index('TOUR_SECTION')+1
                order = []
                for token in tokens[start:]:
                    if token == '-1':
                        break
                    order.append(int(token)-1)
                assert len(order) == len(chosen) and set(order) == set(range(len(chosen)))
                candidate = [chosen[i] for i in order]
                cost = validate(d, candidate)
                assert set(candidate) == set(route)
                if time.perf_counter()-began <= seconds and cost < old:
                    best = candidate
            elapsed = time.perf_counter()-began
        cost = validate(d, best)
        return dict(initial_cost=old, cost=cost, improved=cost < old, seconds=elapsed,
                    timeout=timeout, exit_code=code, log=log,
                    route=[int(d['visits'][v]) for v in best],
                    matrix_cache_created=created, matrix_cache_build_seconds=build_seconds)
