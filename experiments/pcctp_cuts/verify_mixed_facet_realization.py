"""Exhaustive selection-level audit of a mixed realization of the facet family.

No solver or held-out instances. The general equivalence is proved separately;
these checks can falsify indexing and small-case rank claims.
"""
from itertools import combinations
from pathlib import Path
import json


def check(q):
    sizes = [1]*(q-2)+[2,2]
    blocks = []
    for size in sizes:
        start = sum(map(len,blocks))
        blocks.append(set(range(start,start+size)))
    pairs = list(combinations(range(q),2))
    mixed_pairs = [p for p in pairs if p != (q-2,q-1)]
    n = sum(sizes)
    quota = q-1
    ranks = {'coverage_only':q+1,'prize_only':q+1,'joint':q+1}
    counts = {'coverage_only':0,'prize_only':0,'joint':0,'original':0}
    for mask in range(1<<n):
        selected = {v for v in range(n) if mask>>v&1}
        active = sum(bool(selected & b) for b in blocks)
        coverage = all(selected & (blocks[i]|blocks[j]) for i,j in mixed_pairs)
        prize = len(selected)>=quota
        original = all(selected & (blocks[i]|blocks[j]) for i,j in pairs)
        joint = coverage and prize
        assert original == joint, (q,mask)
        for key, feasible in [('coverage_only',coverage),('prize_only',prize),('joint',joint),('original',original)]:
            if feasible:
                counts[key]+=1
                if key in ranks:
                    ranks[key] = min(ranks[key],active)
    assert ranks['coverage_only']==q-2
    assert ranks['joint']==q-1
    assert ranks['prize_only']<q-1
    return dict(q=q, block_sizes=sizes, unit_prizes=True, quota=quota,
                omitted_coverage_pair=[q-1,q], selection_checks=1<<n,
                selection_families_identical=True, block_ranks=ranks,
                feasible_selection_counts=counts)


if __name__=='__main__':
    result=[check(q) for q in range(4,9)]
    Path(__file__).with_name('mixed_facet_realization.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
