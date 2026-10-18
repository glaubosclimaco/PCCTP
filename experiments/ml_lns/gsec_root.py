"""Frozen root-LP GSEC screening, using native HiGHS and exact candidate checks."""
import time
import json
from pathlib import Path
import sys
import numpy as np
import highspy
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_flow
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from experiments.ml_lns.pilot import Repair,read_instance,validate
from experiments.ml_lns.highs_warm_start import start_vector
SCALE=1000000


def separate(repair, values):
    root=min(repair.d['mandatory']);a=[];b=[];cap=[]
    for i,(u,v) in enumerate(repair.edges):
        capacity=int(np.floor(max(0.,values[i])*SCALE))
        if capacity:a.extend([u,v]);b.extend([v,u]);cap.extend([capacity,capacity])
    graph=csr_matrix((np.array(cap,dtype=np.int64),(a,b)),shape=(repair.n,repair.n))
    candidates=[];seen=set()
    for v in range(repair.n):
        if v==root or values[repair.m+v]<=1e-7:continue
        flow=maximum_flow(graph,root,v)
        residual=(graph-flow.flow).tocsr();reachable={root};queue=[root]
        while queue:
            u=queue.pop()
            for k in range(residual.indptr[u],residual.indptr[u+1]):
                w=int(residual.indices[k])
                if residual.data[k]>0 and w not in reachable:
                    reachable.add(w);queue.append(w)
        side=frozenset(set(range(repair.n))-reachable)
        assert v in side and root not in side
        representative=max(side,key=lambda u:values[repair.m+u])
        key=(side,representative)
        if key in seen:continue
        seen.add(key)
        row={i:1. for i,(u,w) in enumerate(repair.edges) if (u in side)!=(w in side)}
        lhs=sum(values[i] for i in row)
        if lhs<2*values[repair.m+representative]-1e-6:
            row[repair.m+representative]=-2.
            candidates.append((side,representative,row))
    return candidates


def probe(path):
    began=time.perf_counter();profile=json.loads(path.read_text())
    d=read_instance(ROOT/'data'/profile['instance'])
    mapping={int(v):i for i,v in enumerate(d['visits'])};route=[mapping[v] for v in profile['route']]
    incumbent=validate(d,route);repair=Repair(d);initial=start_vector(repair,route)
    matrix=repair.constraints.A.tocsr();lp=highspy.HighsLp()
    lp.num_col_=repair.width;lp.num_row_=matrix.shape[0]
    lp.col_cost_=repair.obj;lp.col_lower_=repair.lb;lp.col_upper_=repair.ub
    lp.row_lower_=repair.constraints.lb;lp.row_upper_=repair.constraints.ub
    lp.a_matrix_.format_=highspy.MatrixFormat.kRowwise
    lp.a_matrix_.start_=matrix.indptr;lp.a_matrix_.index_=matrix.indices;lp.a_matrix_.value_=matrix.data
    h=highspy.Highs();h.setOptionValue('output_flag',False);h.setOptionValue('threads',1)
    h.passModel(lp);bounds=[];cuts=[];final_optimal=False
    for iteration in range(4):
        remaining=2.-(time.perf_counter()-began)
        if remaining<=0:break
        h.setOptionValue('time_limit',h.getRunTime()+remaining);h.run()
        if h.getModelStatus()!=highspy.HighsModelStatus.kOptimal:break
        if time.perf_counter()-began>2.:break
        bound=float(h.getObjectiveValue());assert bound<=incumbent+1e-5
        bounds.append(bound);final_optimal=True
        if iteration==3:break
        candidates=separate(repair,np.array(h.getSolution().col_value))
        if not candidates:break
        if time.perf_counter()-began>2.:break
        for side,v,row in candidates:
            assert sum(initial[i]*a for i,a in row.items())>=-1e-8
            indices=np.array(list(row),dtype=np.int32);values=np.array(list(row.values()))
            h.addRow(0.,highspy.kHighsInf,len(row),indices,values)
            cuts.append(dict(side=sorted(side),representative=v))
        final_optimal=False
    elapsed=time.perf_counter()-began
    reduction=None
    if len(bounds)>=2 and incumbent-bounds[0]>1e-6:
        reduction=(bounds[-1]-bounds[0])/(incumbent-bounds[0])
    passes=bool(final_optimal and elapsed<=2. and reduction is not None and reduction>=.1)
    return dict(instance=profile['instance'],geometry=profile['geometry'],mdm_cost=incumbent,
                bounds=bounds,cut_count=len(cuts),cuts=cuts,seconds=elapsed,
                final_lp_optimal=final_optimal,relative_gap_reduction=reduction,passes=passes)


if __name__=='__main__':
    profiles=ROOT/'experiments/original_mdm/move_profiles_monotonic_elite_bounded'
    rows=[]
    for path in sorted(profiles.glob('*.profiled_baseline.json')):
        p=json.loads(path.read_text())
        if not p['completed'] or p['iterations']!=70:continue
        row=probe(path);rows.append(row)
        print(json.dumps({k:v for k,v in row.items() if k!='cuts'}),flush=True)
    result=dict(rows=rows,passing=sum(r['passes'] for r in rows))
    result['gate_passed']=len(rows)==11 and result['passing']>=6
    Path(__file__).with_name('gsec_root_results.json').write_text(json.dumps(result,indent=2)+'\n')
