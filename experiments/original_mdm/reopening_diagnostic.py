"""Release two edges of the MDM's first mined pattern.

The score, the cardinality and the continuation rule are fixed in RULE
before any cost from this comparison is inspected. Development geometries
are the first sample. Training geometries are used only when fewer than
three development geometries have a proven gap.
"""
import hashlib
import json
from pathlib import Path
import random
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from experiments.ml_lns.pilot import Repair, read_instance, validate
from experiments.original_mdm.pattern_selector_diagnostic import compare, parse_patterns, solve_forced

K=2
RULE='Last mining call, seed 42. Release exactly 2 edges of the last mined pattern. Priority is edge cost over the median positive visit cost, plus the average redundancy of its endpoints. An optional vertex is redundant when removing it from the pattern union does not reduce covered demands. Ties prefer higher cost, then the smaller index pair. Random uses SHA-256 of instance, mining index, 2 and 20261018. A geometry counts only when the unrestricted model and both reopenings are proven, the released sets differ, and the fully fixed pattern is strictly worse than the unrestricted optimum. Informed release continues only if it is strictly better than random on more such geometries than the reverse. Fewer than three development geometries with a gap makes this sample inconclusive and sends the unchanged rule to training geometries.'


def compact_pattern(pattern, mapping):
    edges=[]
    for origin, destination in pattern['edges']:
        if origin not in mapping or destination not in mapping:
            raise ValueError('Pattern endpoint is not a visit vertex')
        edges.append(tuple(sorted((mapping[origin], mapping[destination]))))
    unique=list(dict.fromkeys(edges))
    if len(unique)!=len(edges):
        raise ValueError('Pattern repeats an undirected edge')
    return unique


def redundancy(d, vertices):
    covered=sum(bool(vertices & neighborhood) for neighborhood in d['neighborhoods'])
    flags={}
    for vertex in vertices:
        if vertex in d['mandatory']:
            flags[vertex]=0
            continue
        reduced=vertices-{vertex}
        flags[vertex]=int(sum(bool(reduced & neighborhood) for neighborhood in d['neighborhoods'])==covered)
    return flags


def choose_informed(d, edges, median):
    flags=redundancy(d, {vertex for edge in edges for vertex in edge})
    def priority(edge):
        a, b = edge
        return d['costs'][a, b]/median+0.5*(flags[a]+flags[b])
    return sorted(edges, key=lambda edge: (-priority(edge), -d['costs'][edge], edge))[:K]


def choose_random(edges, instance, mining_index):
    seed=int(hashlib.sha256(f'{instance}:{mining_index}:{K}:20261018'.encode()).hexdigest()[:8], 16)
    return sorted(random.Random(seed).sample(list(edges), K))


def evaluate(summary, seconds):
    sidecar=Path(summary['sidecar'])
    bank=json.loads(sidecar.read_text().splitlines()[-1])
    d=read_instance(ROOT/'data'/summary['instance'])
    mapping={int(v): i for i, v in enumerate(d['visits'])}
    for route, cost in zip(bank['routes'], bank['costs']):
        measured=validate(d, [mapping[v] for v in route])
        if abs(measured-cost)>1e-6:
            raise ValueError(f'Elite cost mismatch in {summary["instance"]}')
    patterns=parse_patterns(bank['patterns_text'], bank['n'])
    if not patterns:
        return dict(instance=summary['instance'], geometry=summary['geometry'], status='no_pattern')
    pattern=patterns[-1]
    edges=compact_pattern(pattern, mapping)
    if len(edges)<K+1:
        return dict(instance=summary['instance'], geometry=summary['geometry'], status='pattern_too_small', size=len(edges))
    positive=d['costs'][d['costs']>0]
    median=float(np.median(positive))
    informed=choose_informed(d, edges, median)
    random_edges=choose_random(edges, summary['instance'], bank['mining_index'])
    repair=Repair(d)
    fixed=solve_forced(repair, edges, seconds)
    free=solve_forced(repair, [], seconds)
    informed_forced=[edge for edge in edges if edge not in set(informed)]
    random_forced=[edge for edge in edges if edge not in set(map(tuple, random_edges))]
    informed_result=solve_forced(repair, informed_forced, seconds)
    random_result=solve_forced(repair, random_forced, seconds)
    room=bool(fixed['proven'] and free['proven'] and fixed['cost'] is not None and free['cost'] is not None
              and fixed['cost']>free['cost']+1e-6)
    same=set(map(tuple, informed))==set(map(tuple, random_edges))
    verdict=compare(informed_result, random_result)
    status='compared' if room and not same and verdict!='unavailable' else (
        'saturated' if fixed['proven'] and free['proven'] and not room else
        'same_release' if same else 'unavailable')
    return dict(instance=summary['instance'], geometry=summary['geometry'], status=status,
                mining_index=bank['mining_index'], pattern_size=len(edges), support=pattern['support'],
                room=room, same_release=same, informed_vs_random=verdict,
                costs=dict(fixed=fixed['cost'], unrestricted=free['cost'], informed=informed_result['cost'],
                           random=random_result['cost'], mdm=summary.get('final_cost')),
                proven=dict(fixed=fixed['proven'], unrestricted=free['proven'],
                            informed=informed_result['proven'], random=random_result['proven']))


def main():
    role=sys.argv[1] if len(sys.argv)>1 else 'development'
    directory=HERE/'elite_banks' if role=='development' else HERE/'elite_banks_train'
    rows=[evaluate(summary, 20.) for summary in json.loads((directory/'collection.json').read_text())]
    countable=[row for row in rows if row['status']=='compared']
    better=sum(row['informed_vs_random']=='left' for row in countable)
    worse=sum(row['informed_vs_random']=='right' for row in countable)
    tied=sum(row['informed_vs_random']=='tie' for row in countable)
    result=dict(rule=RULE, role=role, release_count=K,
                code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                geometries=len(rows), countable=len(countable), informed_better=better,
                random_better=worse, tied=tied,
                continues=len(countable)>=3 and better>worse, inconclusive=len(countable)<3,
                rows=rows)
    target=HERE/('reopening_diagnostic.json' if role=='development' else 'reopening_diagnostic_train.json')
    target.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k!='rows'}, indent=2))
    for row in rows:
        print(json.dumps({k: v for k, v in row.items() if k!='proven'}))


if __name__=='__main__':
    main()
