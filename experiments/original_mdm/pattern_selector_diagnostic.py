"""Matched-cardinality pattern selectors on a fixed elite bank.

The formula and the comparison rule are fixed in this file before the costs
are read. Coverage and prize are properties of the vertex union. A larger
pattern is not allowed to win merely by forcing more edges: every solved
pattern in a bank has the same number of edges, the size of the largest
pattern returned by the miner.

No selector weight was fit on these banks. Development geometries only.
"""
import hashlib
import json
import re
from pathlib import Path
import sys
import time

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from experiments.ml_lns.pilot import Repair, read_instance, validate

PATTERN=re.compile(r'(\d+);(\d+);([0-9 ]*)')
RULE='Last mining call of each development geometry, seed 42. Quality comparisons use only patterns whose edge count equals the largest pattern in that call. Utility is the coverage fraction of the vertex union, plus the prize fraction capped at 1, minus the fraction of optional vertices that can be removed without losing coverage, minus the pattern cost divided by its edge count times the median positive visit cost. Each score breaks ties by the earlier file index and by nothing else. The MDM control is the last pattern in the miner file, and its cost is compared only when that pattern has the same edge count. A proven infeasible model is worse than a proven finite cost. A pair counts only when both outcomes are proven. Random uses SHA-256 of instance, mining index and 20261018.'


def parse_patterns(text, size):
    patterns=[]
    for match in PATTERN.finditer(text):
        length, support = int(match.group(1)), int(match.group(2))
        items=[int(token) for token in match.group(3).split()]
        if length==0: continue
        if len(items)!=length: raise ValueError('Pattern length does not match its items')
        edges=[divmod(item, size) for item in items]
        patterns.append(dict(size=length, support=support, items=items, edges=edges))
    return patterns


def utility(d, compact_edges, median_cost):
    vertices=set(v for edge in compact_edges for v in edge)
    covered=sum(bool(vertices & neighborhood) for neighborhood in d['neighborhoods'])
    coverage=covered/max(1, len(d['neighborhoods']))
    prize=min(1., float(sum(d['prizes'][v] for v in vertices))/max(d['quota'], 1.))
    optional=[v for v in vertices if v not in d['mandatory']]
    redundant=0
    for vertex in optional:
        reduced=vertices-{vertex}
        if sum(bool(reduced & neighborhood) for neighborhood in d['neighborhoods'])==covered:
            redundant+=1
    redundancy=redundant/max(1, len(optional)) if optional else 0.
    cost=sum(d['costs'][a, b] for a, b in compact_edges)
    return coverage+prize-redundancy-cost/(len(compact_edges)*median_cost)


def choose(patterns, score, instance, mining_index, name):
    if name=='random':
        seed=int(hashlib.sha256(f'{instance}:{mining_index}:20261018'.encode()).hexdigest()[:8], 16)
        return __import__('random').Random(seed).randrange(len(patterns))
    return max(range(len(patterns)), key=lambda i: (score(patterns[i]), -i))


def solve_forced(repair, forced, seconds):
    degree={}
    for a, b in forced:
        degree[a]=degree.get(a, 0)+1
        degree[b]=degree.get(b, 0)+1
    if any(value>2 for value in degree.values()) or any(a==b for a, b in forced):
        return dict(outcome='structural_infeasible', cost=None, proven=True, seconds=0.)
    lower, upper = repair.lb.copy(), repair.ub.copy()
    forced_set=set(map(tuple, forced))
    for index, edge in enumerate(repair.edges):
        if edge in forced_set: lower[index]=upper[index]=1.
    for vertex in degree: lower[repair.m+vertex]=1.
    began=time.perf_counter()
    result=milp(repair.obj, integrality=repair.integrality, bounds=Bounds(lower, upper),
                constraints=repair.constraints, options={'time_limit': seconds, 'mip_rel_gap': 0., 'disp': False})
    elapsed=time.perf_counter()-began
    if result.status==2:
        return dict(outcome='infeasible', cost=None, proven=True, seconds=elapsed)
    if result.status!=0 or result.x is None:
        return dict(outcome='unproven', cost=None, proven=False, seconds=elapsed, status=int(result.status))
    vector=result.x
    if np.max(np.abs(vector[:repair.f0]-np.rint(vector[:repair.f0])))>1e-5:
        raise ValueError('Nonintegral forced solution')
    adjacency={v: [] for v in range(repair.n) if vector[repair.m+v]>.5}
    for index, (a, b) in enumerate(repair.edges):
        if vector[index]>.5:
            adjacency[a].append(b); adjacency[b].append(a)
    if any(len(neighbors)!=2 for neighbors in adjacency.values()):
        raise ValueError('Forced solution has an invalid degree')
    route=[min(repair.d['mandatory'])]; previous=None
    while True:
        nxt=next(v for v in adjacency[route[-1]] if v!=previous)
        if nxt==route[0]: break
        previous=route[-1]; route.append(nxt)
    decoded=set(route)
    if decoded!=set(adjacency): raise ValueError('Forced solution is disconnected')
    cost=validate(repair.d, route)
    used={tuple(sorted((route[i-1], v))) for i, v in enumerate(route)}
    if not forced_set<=used: raise ValueError('Forced edge missing from the solution')
    if abs(cost-float(result.fun))>1e-5: raise ValueError('Forced objective mismatch')
    return dict(outcome='optimal', cost=cost, proven=True, seconds=elapsed)


def compare(left, right):
    if not left['proven'] or not right['proven']: return 'unavailable'
    if left['cost'] is None and right['cost'] is None: return 'tie'
    if left['cost'] is None: return 'right'
    if right['cost'] is None: return 'left'
    if abs(left['cost']-right['cost'])<=1e-6: return 'tie'
    return 'left' if left['cost']<right['cost'] else 'right'


def main():
    directory=HERE/'elite_banks'
    rows=[]
    for summary in json.loads((directory/'collection.json').read_text()):
        sidecar=Path(summary['sidecar'])
        if not sidecar.exists() or summary['elite_lines']==0:
            rows.append(dict(instance=summary['instance'], geometry=summary['geometry'], status='no_mining'))
            continue
        bank=json.loads(sidecar.read_text().splitlines()[-1])
        d=read_instance(ROOT/'data'/summary['instance'])
        mapping={int(v): i for i, v in enumerate(d['visits'])}
        for route, cost in zip(bank['routes'], bank['costs']):
            measured=validate(d, [mapping[v] for v in route])
            if abs(measured-cost)>1e-6: raise ValueError(f'Elite cost mismatch in {summary["instance"]}')
        patterns=parse_patterns(bank['patterns_text'], bank['n'])
        if not patterns:
            rows.append(dict(instance=summary['instance'], geometry=summary['geometry'], status='no_pattern',
                             elite_routes=len(bank['routes']), mining_index=bank['mining_index']))
            continue
        positive=d['costs'][d['costs']>0]
        median=float(np.median(positive)) if len(positive) else 1.
        frequency={}
        for route in bank['routes']:
            for i, vertex in enumerate(route):
                frequency[(route[i-1], vertex)]=frequency.get((route[i-1], vertex), 0)+1
        compact={}
        for pattern in patterns:
            edges=[]
            for origin, destination in pattern['edges']:
                if origin not in mapping or destination not in mapping:
                    raise ValueError('Pattern endpoint is not a visit vertex')
                edges.append(tuple(sorted((mapping[origin], mapping[destination]))))
            pattern['compact']=list(dict.fromkeys(edges))
            pattern['simple']=len(pattern['compact'])==pattern['size']
            pattern['marginal']=sum(frequency[edge] for edge in pattern['edges'])
            pattern['utility']=utility(d, pattern['compact'], median)
        width=max(pattern['size'] for pattern in patterns)
        pool=[i for i, pattern in enumerate(patterns) if pattern['size']==width and pattern['simple']]
        if len(pool)<2:
            rows.append(dict(instance=summary['instance'], geometry=summary['geometry'], status='no_size_contrast',
                             mining_index=bank['mining_index'], patterns=len(patterns), max_size=width, pool=len(pool)))
            continue
        names=['size', 'marginal', 'support', 'utility', 'random', 'mdm_first']
        scores=dict(size=lambda p: p['size'], marginal=lambda p: p['marginal'],
                    support=lambda p: p['support'], utility=lambda p: p['utility'])
        chosen={}
        for name in names:
            if name=='mdm_first':
                chosen[name]=len(patterns)-1
            else:
                pick=choose([patterns[i] for i in pool], scores.get(name, lambda pattern: 0.), summary['instance'], bank['mining_index'], name)
                chosen[name]=pool[pick]
        repair=Repair(d)
        solved={}
        selections={}
        for name, index in chosen.items():
            pattern=patterns[index]
            if pattern['size']!=width or not pattern['simple']:
                selections[name]=dict(index=index, size=pattern['size'], outcome='different_cardinality', cost=None, proven=False)
                continue
            key=tuple(pattern['compact'])
            if key not in solved: solved[key]=solve_forced(repair, key, 2.)
            selections[name]=dict(index=index, size=pattern['size'], support=pattern['support'],
                                  marginal=pattern['marginal'], utility=pattern['utility'], **solved[key])
        def paired(name):
            if selections[name].get('outcome')=='different_cardinality': return 'different_cardinality'
            return compare(selections['utility'], selections[name])
        rows.append(dict(instance=summary['instance'], geometry=summary['geometry'], status='compared',
                         mining_index=bank['mining_index'], elite_routes=len(bank['routes']),
                         patterns=len(patterns), pool=len(pool), max_size=width,
                         utility_vs_marginal=paired('marginal'),
                         utility_vs_support=paired('support'),
                         utility_vs_random=paired('random'),
                         utility_vs_mdm_first=paired('mdm_first'),
                         utility_vs_size=paired('size'),
                         same_pattern_as_marginal=selections['utility']['index']==selections['marginal']['index'],
                         same_pattern_as_support=selections['utility']['index']==selections['support']['index'],
                         same_pattern_as_size=selections['utility']['index']==selections['size']['index'],
                         same_pattern_as_mdm_first=selections['utility'].get('index')==selections['mdm_first'].get('index'),
                         selections=selections))
        print(json.dumps({k: v for k, v in rows[-1].items() if k!='selections'}), flush=True)
    def tally(key):
        values=[row[key] for row in rows if row['status']=='compared']
        return {name: values.count(name) for name in ['left', 'right', 'tie', 'unavailable', 'different_cardinality']}
    result=dict(rule=RULE, code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                geometries=len(rows), compared=sum(row['status']=='compared' for row in rows),
                no_mining=sum(row['status']=='no_mining' for row in rows),
                no_pattern=sum(row['status']=='no_pattern' for row in rows),
                no_size_contrast=sum(row['status']=='no_size_contrast' for row in rows),
                utility_vs_marginal=tally('utility_vs_marginal'),
                utility_vs_support=tally('utility_vs_support'),
                utility_vs_random=tally('utility_vs_random'),
                utility_vs_mdm_first=tally('utility_vs_mdm_first'),
                utility_vs_size=tally('utility_vs_size'),
                rows=rows)
    (HERE/'pattern_selector_diagnostic.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k!='rows'}, indent=2))


if __name__=='__main__':
    main()
