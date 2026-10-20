"""Integer proofs with and without the frozen root rank cuts.

The 30 second limit and the 80 percent speed rule are fixed in
rank_partition_root.md before these proofs are read.
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
from experiments.ml_lns.pilot import Repair, read_instance, validate
from experiments.pcctp_cuts.separate_rank_root import ROUNDS, cross_sum, optimize, partition, rank

LIMIT=30.


def accumulate(repair, route):
    constraints=repair.constraints
    values, elapsed, status=optimize(repair, constraints)
    cuts=0
    if values is None:
        return None, elapsed, cuts
    for _ in range(ROUNDS):
        label, blocks, base=partition(repair.d, repair.edges, values)
        rho=rank(repair.d, label, blocks, base)
        if not blocks or rho is None or rho<1: break
        if sum(label[route[i-1]]!=label[route[i]] for i in range(len(route)))<rho+1:
            raise RuntimeError('Feasible route violates a generated cut')
        if cross_sum(repair.edges, label, values)>=rho+1-1e-6: break
        indices=np.array([index for index, (a, b) in enumerate(repair.edges) if label[a]!=label[b]])
        row=csr_matrix((np.ones(len(indices)), (np.zeros(len(indices), dtype=int), indices)), shape=(1, repair.width))
        constraints=LinearConstraint(vstack([constraints.A, row]), np.r_[constraints.lb, rho+1], np.r_[constraints.ub, np.inf])
        values, round_elapsed, status=optimize(repair, constraints)
        elapsed+=round_elapsed
        if values is None: return None, elapsed, cuts
        cuts+=1
    return constraints, elapsed, cuts


def prove(repair, constraints, seconds):
    began=time.perf_counter()
    result=milp(repair.obj, integrality=repair.integrality, bounds=Bounds(repair.lb, repair.ub),
                constraints=constraints, options={'time_limit': seconds, 'mip_rel_gap': 0., 'disp': False})
    elapsed=time.perf_counter()-began
    proven=result.status==0 and result.x is not None
    cost=None
    if proven:
        vector=np.rint(result.x[:repair.f0])
        if np.max(np.abs(result.x[:repair.f0]-vector))>1e-5: raise ValueError('Nonintegral proof')
        adjacency={v: [] for v in range(repair.n) if vector[repair.m+v]>.5}
        for index, (a, b) in enumerate(repair.edges):
            if vector[index]>.5:
                adjacency[a].append(b); adjacency[b].append(a)
        tour=[min(repair.d['mandatory'])]; previous=None
        while True:
            nxt=next(v for v in adjacency[tour[-1]] if v!=previous)
            if nxt==tour[0]: break
            previous=tour[-1]; tour.append(nxt)
        cost=validate(repair.d, tour)
        if abs(cost-float(result.fun))>1e-4: raise ValueError('Proof objective mismatch')
    return dict(proven=proven, cost=cost, seconds=elapsed, status=int(result.status))


def main():
    banks=json.loads((ROOT/'experiments/original_mdm/elite_banks/collection.json').read_text())
    rows=[]
    for summary in banks:
        d=read_instance(ROOT/'data'/summary['instance'])
        mapping={int(v): i for i, v in enumerate(d['visits'])}
        native=json.loads((ROOT/'experiments/original_mdm/elite_banks'/f"{summary['instance']}.seed42.json").read_text())
        route=[mapping[v] for v in native['route']]
        repair=Repair(d)
        plain=prove(repair, repair.constraints, LIMIT)
        constraints, separation, cuts=accumulate(repair, route)
        if constraints is None:
            cut=dict(proven=False, cost=None, seconds=separation, status='separation_failed')
        else:
            cut=prove(repair, constraints, LIMIT)
            cut['seconds']=separation+cut['seconds']
            cut['cuts']=cuts
            cut['separation_seconds']=separation
        if plain['proven'] and cut['proven'] and abs(plain['cost']-cut['cost'])>1e-4:
            raise RuntimeError(f'Different optima in {summary["geometry"]}')
        both=plain['proven'] and cut['proven']
        faster=both and cut['seconds']<=0.8*plain['seconds']
        slower=both and plain['seconds']<=0.8*cut['seconds']
        rows.append(dict(geometry=summary['geometry'], plain=plain, with_cuts=cut,
                         only_plain=plain['proven'] and not cut['proven'],
                         only_cuts=cut['proven'] and not plain['proven'],
                         faster=faster, slower=slower))
        print(json.dumps(dict(geometry=summary['geometry'], plain_proven=plain['proven'], plain_cost=plain['cost'],
                             plain_seconds=plain['seconds'], cut_proven=cut['proven'], cut_cost=cut['cost'],
                             cut_seconds=cut['seconds'], cuts=cut.get('cuts'))), flush=True)
    only_cuts=sum(row['only_cuts'] for row in rows)
    only_plain=sum(row['only_plain'] for row in rows)
    faster=sum(row['faster'] for row in rows)
    slower=sum(row['slower'] for row in rows)
    result=dict(limit_seconds=LIMIT, geometries=len(rows), only_cuts=only_cuts, only_plain=only_plain,
                faster=faster, slower=slower, continues=only_cuts>only_plain or (only_cuts==only_plain and faster>slower),
                rows=rows)
    (HERE/'rank_certification.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k!='rows'}, indent=2))


if __name__=='__main__':
    main()
