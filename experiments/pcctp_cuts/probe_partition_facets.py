"""Finite exact rank certificates for the complete pair-coverage construction.

Sampling only finds integer generators; all arithmetic proving independence is
modular. Failure to reach the upper bound is explicitly inconclusive.
No benchmark, held-out geometry, or performance comparison is involved.
"""
from itertools import combinations
from pathlib import Path
import json
import random
from verify_partition_face import ModularBasis, PRIME


def probe(sizes, seed=20261018, limit=5000):
    rng = random.Random(seed)
    blocks = [[0]]
    for size in sizes:
        start = sum(map(len, blocks))
        blocks.append(list(range(start, start+size)))
    n = sum(map(len, blocks))
    edges = list(combinations(range(n), 2))
    index = {e: i for i, e in enumerate(edges)}
    group = {v: i for i, b in enumerate(blocks) for v in b}
    d = len(edges)+n-1
    rows = [[int(v in e) for e in edges]+[-2*int(i == v) for i in range(1,n)] for v in range(n)]
    eq = ModularBasis()
    for row in rows:
        eq.add(row)
    full_upper = d-eq.rank
    eq.add([int(group[a] != group[b]) for a,b in edges]+[0]*(n-1))
    # For q=3, the sum of the three coverage-forced pair cuts equals
    # twice the partition cut. Equality in the latter forces all three.
    if len(sizes) == 3:
        for i,j in combinations(range(1,4),2):
            union = set(blocks[i]+blocks[j])
            eq.add([int((a in union)!=(b in union)) for a,b in edges]+[0]*(n-1))
    face_upper = d-eq.rank
    bases = [ModularBasis(), ModularBasis()]
    first = [None, None]
    generators = [[], []]
    q = len(sizes)
    for step in range(limit):
        for tight in (False, True):
            omitted = rng.randrange(q+1) if not tight else rng.randrange(1,q+1)
            active = [i for i in range(1,q+1) if i != omitted]
            pieces = []
            for i in active:
                selected = rng.sample(blocks[i], rng.randrange(1,len(blocks[i])+1))
                pieces.append(selected)
            rng.shuffle(pieces)
            order = [v for piece in pieces for v in piece]
            if not tight:
                rng.shuffle(order)
            route = [0]+order
            selected_set = set(route)
            assert all(selected_set.intersection(blocks[i]+blocks[j]) for i,j in combinations(range(1,q+1),2))
            point = [0]*d
            cross = 0
            for a,b in zip(route, route[1:]+route[:1]):
                point[index[tuple(sorted((a,b)))]] = 1
                cross += group[a] != group[b]
            for v in order:
                point[len(edges)+v-1] = 1
            assert cross >= q
            if tight:
                assert cross == q
            targets = (0,1) if cross == q else (0,)
            for t in targets:
                if first[t] is None:
                    first[t] = point
                    generators[t].append(route)
                elif bases[t].add([a-b for a,b in zip(point,first[t])]):
                    generators[t].append(route)
        if bases[0].rank == full_upper and bases[1].rank == face_upper:
            break
    certified = bases[0].rank == full_upper and bases[1].rank == face_upper
    return dict(block_sizes=sizes, vertices=n, prime=PRIME, seed=seed,
                samples_per_mode=step+1, full_rank_lower=bases[0].rank,
                full_dimension_upper=full_upper, face_rank_lower=bases[1].rank,
                face_dimension_upper=face_upper,
                status=('finite_facet_certified' if full_upper-face_upper == 1 else 'finite_nonfacet_certified') if certified else 'inconclusive',
                full_generating_routes=generators[0], face_generating_routes=generators[1])


if __name__ == '__main__':
    cases = [[2]*q for q in range(3,7)]+[[3,2,2,2],[1,2,2,2]]
    results = [probe(sizes) for sizes in cases]
    Path(__file__).with_name('partition_facet_probe.json').write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps([{k:v for k,v in r.items() if not k.endswith('routes')} for r in results],indent=2))
