"""Isolated conditional phase adapter; no collection entry point or queue.

The future batch owner must hold the shared study lock. The temporary polish
hook is process-local and phases must be sequential within that process.
"""
import hashlib
import json
import math
import sys
import time
from pathlib import Path

from experiments.lkh import ruin_repair_absolute as core
from experiments.lkh.cached_polish import CachedPolisher

ROOT = Path(__file__).resolve().parents[2]
ORIGINAL_POLISH = core.polish
SOURCES = ('experiments/lkh/cached30_phase.py', 'experiments/lkh/cached_polish.py',
           'experiments/lkh/matrix_text_cache.py', 'experiments/lkh/ruin_repair_absolute.py',
           'experiments/lkh/screen.py', 'experiments/ml_lns/pilot.py',
           'experiments/ml_branching/oracle_probe.py')


def source_hashes():
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in SOURCES}


def cache_bytes(cache):
    """Approximate retained Python object bytes, not process RSS or peak memory."""
    if cache is None:
        return 0
    owned = [cache, cache.__dict__, cache.rows]
    owned.extend(cache.rows)
    owned.extend(token for row in cache.rows for token in row)
    unique = {id(obj): obj for obj in owned}
    return sum(sys.getsizeof(obj) for obj in unique.values())


def run_phase(start, repair_seed, arm, budget, target):
    began = time.clock_gettime(time.CLOCK_MONOTONIC)
    assert arm in ('uncached', 'cached') and math.isfinite(budget) and budget > 0
    assert start['source_seed'] in (42, 43, 44)
    assert repair_seed == start['source_seed']+49
    target = Path(target)
    failure = target.with_suffix('.failed.json')
    attempt = target.with_suffix('.attempt.json')
    assert not any(p.exists() for p in (target, failure, attempt)), 'Review existing phase artifacts before reuse'
    assert core.polish is ORIGINAL_POLISH, 'Another phase hook is already installed'
    target.parent.mkdir(parents=True, exist_ok=True)
    polisher = CachedPolisher() if arm == 'cached' else None
    builds, calls = [], 0
    # Leave this marker on interruption, post-search assertion or write failure too.
    # A completed result is committed before the attempt marker is removed.
    with attempt.open('x') as handle:
        json.dump(dict(instance=start['instance'], arm=arm, source_seed=start['source_seed'],
                       repair_seed=repair_seed, parent_started=began), handle)

    def invoke(d, route, began=None, seconds=1.):
        nonlocal calls
        call_began = time.clock_gettime(time.CLOCK_MONOTONIC)
        result = (polisher if polisher else ORIGINAL_POLISH)(d, route, began=began, seconds=seconds)
        call_end = time.clock_gettime(time.CLOCK_MONOTONIC)
        calls += 1
        if result.get('matrix_cache_created'):
            builds.append(dict(call_started=call_began, call_finished=call_end,
                               build_seconds=result['matrix_cache_build_seconds'],
                               allowance_seconds=seconds))
        return result

    try:
        core.polish = invoke
        raw = core.run(start['instance'], 'random', budget,
                       initial_override=start['route'][:], seed=repair_seed)
    except Exception as error:
        # Preserve an explicit failed attempt; a future batch must not auto-reuse it.
        failure.write_text(json.dumps(dict(instance=start['instance'], arm=arm,
            source_seed=start['source_seed'], repair_seed=repair_seed,
            error_type=type(error).__name__, error=str(error)), indent=2)+'\n')
        raise
    finally:
        core.polish = ORIGINAL_POLISH
    elapsed = time.clock_gettime(time.CLOCK_MONOTONIC)-began
    assert raw['initial_cost'] == start['cost'] and raw['lkh_calls'] == calls
    cost, route, origin = start['cost'], start['route'][:], 'common_start'
    for event in raw['trace']:
        if event['native_timestamp'] <= began+budget and event['cost'] < cost:
            cost, route, origin = event['cost'], event['route'][:], 'phase_event'
    record = dict(instance=start['instance'], start=start, repair_seed=repair_seed,
                  arm=arm, budget_seconds=budget, parent_started=began,
                  wall_seconds=elapsed, raw_phase=raw, cache_build_events=builds,
                  cache_retained_python_bytes=cache_bytes(polisher.cache) if polisher else 0,
                  memory_scope='Unique owned cache tuple/string/object getsizeof sum; '
                               'excludes input matrix, allocator overhead and process RSS.',
                  cost=cost, route=route, origin=origin, source_sha256=source_hashes(),
                  scope='Conditional common start; generation excluded. Phase setup, input '
                        'loading, cache initialization and search charged to outer event horizon.')
    temp = target.with_suffix('.tmp')
    temp.write_text(json.dumps(record, indent=2)+'\n'); temp.replace(target)
    attempt.unlink()
    return record
