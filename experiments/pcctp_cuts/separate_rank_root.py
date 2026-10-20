"""Root-bound test of the certified rank-partition cut.

The separator and the success count are fixed in the research note
rank_partition_root.md before these bounds are read.
"""
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import csr_matrix, vstack

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from experiments.ml_lns.pilot import Repair, read_instance

THRESHOLD=0.5
MAX_BLOCKS=8
ROUNDS=5


def optimize(repair, constraints):
    began=time.perf_counter()
    result=milp(repair.obj, integrality=np.zeros(repair.width), bounds=Bounds(repair.lb, repair.ub),
                constraints=constraints, options={'time_limit': 20., 'disp': False})
    elapsed=time.perf_counter()-began
    if result.status!=0 or result.x is None:
        return None, elapsed, int(result.status)
    return result.x, elapsed, 0


def partition(d, edges, values):
    parent=list(range(d_n(d)))
    def find(vertex):
        while parent[vertex]!=vertex:
            parent[vertex]=parent[parent[vertex]]
            vertex=parent[vertex]
        return vertex
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra!=rb: parent[rb]=ra
    for index, (a, b) in enumerate(edges):
        if values[index]>=THRESHOLD: union(a, b)
    root=min(d['mandatory'])
    groups={}
    for vertex in range(d_n(d)):
        groups.setdefault(find(vertex), []).append(vertex)
    blocks=[sorted(group) for comp, group in groups.items() if comp!=find(root)]
    blocks.sort(key=lambda block: block[0])
    merged=blocks[MAX_BLOCKS:]
    blocks=blocks[:MAX_BLOCKS]
    base=sorted(groups[find(root)]+[vertex for block in merged for vertex in block])
    label=[-1]*d_n(d)
    for vertex in base: label[vertex]=0
    for index, block in enumerate(blocks, start=1):
        for vertex in block: label[vertex]=index
    return label, blocks, base


def d_n(d):
    return len(d['visits'])


def rank(d, label, blocks, base):
    prizes=d['prizes']
    base_prize=float(sum(prizes[v] for v in base))
    block_prize=[float(sum(prizes[v] for v in block)) for block in blocks]
    forced=[index for index, block in enumerate(blocks) if any(v in d['mandatory'] for v in block)]
    needs=[]
    base_set=set(base)
    block_sets=[set(block) for block in blocks]
    for neighborhood in d['neighborhoods']:
        if neighborhood & base_set: continue
        covers=[index for index, block in enumerate(block_sets) if neighborhood & block]
        if not covers: return None
        needs.append(covers)
    best=None
    width=len(blocks)
    for mask in range(1<<width):
        if any((mask>>index)&1==0 for index in forced): continue
        prize=base_prize+sum(block_prize[index] for index in range(width) if (mask>>index)&1)
        if prize+1e-8<d['quota']: continue
        if any(not any((mask>>index)&1 for index in covers) for covers in needs): continue
        card=bin(mask).count('1')
        if best is None or card<best: best=card
    return best


def cross_sum(edges, label, values):
    total=0.
    for index, (a, b) in enumerate(edges):
        if label[a]!=label[b]: total+=float(values[index])
    return total


def integer_cross(route, label):
    return sum(label[route[i-1]]!=label[route[i]] for i in range(len(route)))


def main():
    banks=json.loads((ROOT/'experiments/original_mdm/elite_banks/collection.json').read_text())
    proven={}
    reopening=json.loads((ROOT/'experiments/original_mdm/reopening_diagnostic.json').read_text())
    for row in reopening['rows']:
        if row['proven']['unrestricted'] and row['costs']['unrestricted'] is not None:
            proven[row['geometry']]=row['costs']['unrestricted']
    rows=[]
    for summary in banks:
        d=read_instance(ROOT/'data'/summary['instance'])
        mapping={int(v): i for i, v in enumerate(d['visits'])}
        native=json.loads((ROOT/'experiments/original_mdm/elite_banks'/f"{summary['instance']}.seed42.json").read_text())
        route=[mapping[v] for v in native['route']]
        repair=Repair(d)
        constraints=repair.constraints
        values, elapsed, status=optimize(repair, constraints)
        if values is None:
            rows.append(dict(geometry=summary['geometry'], status='lp_failed', lp_status=status))
            print(rows[-1], flush=True)
            continue
        initial=float(values@repair.obj)
        cuts=0
        rejected=False
        rounds=[]
        for _ in range(ROUNDS):
            label, blocks, base=partition(d, repair.edges, values)
            rho=rank(d, label, blocks, base)
            if not blocks or rho is None or rho<1:
                rounds.append(dict(blocks=len(blocks), rho=rho, added=False))
                break
            current=cross_sum(repair.edges, label, values)
            if integer_cross(route, label)<rho+1:
                rejected=True
                rounds.append(dict(blocks=len(blocks), rho=rho, integer_cross=integer_cross(route, label), added=False, invalid=True))
                break
            if current>=rho+1-1e-6:
                rounds.append(dict(blocks=len(blocks), rho=rho, fractional_cross=current, added=False))
                break
            indices=np.array([index for index, (a, b) in enumerate(repair.edges) if label[a]!=label[b]])
            row=csr_matrix((np.ones(len(indices)), (np.zeros(len(indices), dtype=int), indices)), shape=(1, repair.width))
            constraints=LinearConstraint(vstack([constraints.A, row]),
                                         np.r_[constraints.lb, rho+1], np.r_[constraints.ub, np.inf])
            values, round_elapsed, status=optimize(repair, constraints)
            elapsed+=round_elapsed
            if values is None:
                rows.append(dict(geometry=summary['geometry'], status='lp_failed_after_cut', cuts=cuts))
                break
            cuts+=1
            rounds.append(dict(blocks=len(blocks), rho=rho, fractional_cross=current, added=True))
        else:
            values=values
        if rows and rows[-1]['geometry']==summary['geometry'] and rows[-1]['status'].startswith('lp_failed'):
            print(rows[-1], flush=True)
            continue
        final=float(values@repair.obj)
        optimum=proven.get(summary['geometry'])
        rows.append(dict(geometry=summary['geometry'], status='invalid_cut' if rejected else 'ok',
                         initial_bound=initial, final_bound=final, improved=final>initial+1e-5,
                         cuts=cuts, seconds=elapsed, proven_optimum=optimum,
                         gap_before=None if optimum is None else (optimum-initial)/max(optimum, 1.),
                         gap_after=None if optimum is None else (optimum-final)/max(optimum, 1.),
                         rounds=rounds))
        print({k: v for k, v in rows[-1].items() if k!='rounds'}, flush=True)
    improved=sum(row.get('improved') is True and row['status']=='ok' for row in rows)
    result=dict(threshold=THRESHOLD, max_blocks=MAX_BLOCKS, rounds=ROUNDS,
                geometries=len(rows), improved=improved, useful=improved>=5,
                invalid=sum(row['status']=='invalid_cut' for row in rows), rows=rows)
    (HERE/'rank_root_diagnostic.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k!='rows'}, indent=2))


if __name__=='__main__':
    main()
