"""Native HiGHS SCF diagnostic with validated complete starts and callbacks."""
import time
STARTUP=time.perf_counter()
import json
from pathlib import Path
import sys
import numpy as np
import highspy
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from experiments.ml_lns.pilot import Repair, read_instance, validate


def start_vector(repair, route):
    vector=np.zeros(repair.width)
    index={e:i for i,e in enumerate(repair.edges)}
    root=min(repair.d['mandatory']); pos=route.index(root)
    ordered=route[pos:]+route[:pos]
    for v in route:vector[repair.m+v]=1.
    for k,a in enumerate(ordered):
        b=ordered[(k+1)%len(ordered)]; edge=tuple(sorted((a,b)))
        e=index[edge];vector[e]=1.
        flow=len(ordered)-1-k
        vector[repair.f0+2*e+int(a>b)]=flow
    value=repair.constraints.A@vector
    assert np.all(vector>=repair.lb-1e-8) and np.all(vector<=repair.ub+1e-8)
    assert np.all(value>=repair.constraints.lb-1e-8) and np.all(value<=repair.constraints.ub+1e-8)
    assert abs(repair.obj@vector-validate(repair.d,route))<1e-6
    return vector


def decode(repair, vector):
    if np.max(abs(vector[:repair.f0]-np.rint(vector[:repair.f0])))>1e-5:
        raise ValueError('Fractional incumbent')
    adj={v:[] for v in range(repair.n) if vector[repair.m+v]>.5}
    for i,(a,b) in enumerate(repair.edges):
        if vector[i]>.5:adj[a].append(b);adj[b].append(a)
    if any(len(ns)!=2 for ns in adj.values()):raise ValueError('Degree')
    tour=[min(repair.d['mandatory'])];previous=None
    while True:
        nxt=next(v for v in adj[tour[-1]] if v!=previous)
        if nxt==tour[0]:break
        if nxt in tour:raise ValueError('Disconnected or repeated cycle')
        previous=tour[-1];tour.append(nxt)
    if set(tour)!=set(adj):raise ValueError('Disconnected cycle')
    cost=validate(repair.d,tour)
    if abs(cost-repair.obj@vector)>1e-4:raise ValueError('Objective mismatch')
    return tour,cost


def solve(d, route, warm, seconds=10., began=None, root_cuts=None, allowed_edges=None):
    began=time.perf_counter() if began is None else began
    repair=Repair(d);initial=start_vector(repair,route)
    if allowed_edges is not None:
        allowed_edges={tuple(sorted(map(int,e))) for e in allowed_edges}
        if not allowed_edges<=set(repair.edges):raise ValueError('Invalid pool edge')
        incumbent_edges={tuple(sorted((route[i-1],v))) for i,v in enumerate(route)}
        if not incumbent_edges<=allowed_edges:raise ValueError('Pool excludes incumbent edge')
        for i,e in enumerate(repair.edges):
            if e not in allowed_edges:repair.ub[i]=0.

    if root_cuts:
        from scipy.optimize import LinearConstraint
        from scipy.sparse import csr_matrix,vstack
        additions=[]
        for cut in root_cuts:
            side=set(cut['side']);v=cut['representative']
            assert min(d['mandatory']) not in side and v in side
            row={i:1. for i,(a,b) in enumerate(repair.edges) if (a in side)!=(b in side)}
            row[repair.m+v]=-2.
            assert sum(initial[i]*a for i,a in row.items())>=-1e-8
            additions.append(csr_matrix(([row[i] for i in row],([0]*len(row),list(row))),shape=(1,repair.width)))
        repair.constraints=LinearConstraint(vstack([repair.constraints.A]+additions),
                    np.r_[repair.constraints.lb,np.zeros(len(additions))],
                    np.r_[repair.constraints.ub,np.full(len(additions),np.inf)])
    matrix=repair.constraints.A.tocsr()
    lp=highspy.HighsLp();lp.num_col_=repair.width;lp.num_row_=matrix.shape[0]
    lp.col_cost_=repair.obj;lp.col_lower_=repair.lb;lp.col_upper_=repair.ub
    lp.row_lower_=repair.constraints.lb;lp.row_upper_=repair.constraints.ub
    lp.integrality_=[highspy.HighsVarType.kInteger if x else highspy.HighsVarType.kContinuous for x in repair.integrality]
    lp.a_matrix_.format_=highspy.MatrixFormat.kRowwise
    lp.a_matrix_.start_=matrix.indptr;lp.a_matrix_.index_=matrix.indices;lp.a_matrix_.value_=matrix.data
    h=highspy.Highs();h.setOptionValue('output_flag',False);h.setOptionValue('threads',1)
    h.setOptionValue('random_seed',42);h.setOptionValue('mip_rel_gap',0.)
    assert h.passModel(lp)==highspy.HighsStatus.kOk
    best=list(route);bestcost=validate(d,best);events=[];errors=[];calls=0
    def capture(event):
        nonlocal best,bestcost,calls
        calls+=1
        if time.perf_counter()-began>seconds:return
        try:
            vector=np.array(event.data_out.mip_solution,copy=True)
            tour,cost=decode(repair,vector)
            elapsed=time.perf_counter()-began
            if elapsed<=seconds and cost<bestcost-1e-7:
                best,bestcost=tour,cost;events.append(dict(seconds=elapsed,cost=cost))
        except Exception as error:errors.append(str(error))
    h.cbMipImprovingSolution.subscribe(capture)
    start_status=None
    if warm:start_status=str(h.setSolution(repair.width,np.arange(repair.width,dtype=np.int32),initial))
    remaining=seconds-(time.perf_counter()-began)
    if remaining<=0:raise RuntimeError('Setup exhausted budget')
    h.setOptionValue('time_limit',remaining);h.run()
    # Accept a final result only after decoding within the same wall deadline.
    solution=h.getSolution()
    if solution.value_valid and time.perf_counter()-began<=seconds:
        tour,cost=decode(repair,np.array(solution.col_value))
        if time.perf_counter()-began<=seconds and cost<bestcost-1e-7:
            best,bestcost=tour,cost;events.append(dict(seconds=time.perf_counter()-began,cost=cost))
    if errors:raise RuntimeError('Invalid callback: '+repr(errors))
    elapsed=time.perf_counter()-began;info=h.getInfo()
    optimal=h.getModelStatus()==highspy.HighsModelStatus.kOptimal
    proven=optimal and elapsed<=seconds and abs(info.objective_function_value-bestcost)<1e-5
    return dict(warm_start=warm,start_status=start_status,version=h.version(),
                proof_scope='restricted edge pool' if allowed_edges is not None else 'full PCCTP',
                allowed_edge_count=len(allowed_edges) if allowed_edges is not None else repair.m,
                cost=validate(d,best),initial_cost=validate(d,route),
                improved=bestcost<validate(d,route)-1e-7,seconds=elapsed,
                overrun=max(0.,elapsed-seconds),callback_calls=calls,events=events,
                status=h.modelStatusToString(h.getModelStatus()),proven_optimal=proven,
                nodes=int(info.mip_node_count),dual_bound=float(info.mip_dual_bound),
                route=[int(d['visits'][v]) for v in best])


if __name__=='__main__':
    directory=ROOT/'experiments/original_mdm/move_profiles_monotonic_elite_bounded'
    out=Path(__file__).with_name('highs_warm_start_results.json')
    result=dict(import_seconds=time.perf_counter()-STARTUP,rows=[],excluded_incomplete=[])
    for path in sorted(directory.glob('*.profiled_baseline.json')):
        profile=json.loads(path.read_text())
        if not profile['completed'] or profile['iterations']!=70:
            result['excluded_incomplete'].append(profile['instance']);continue
        row=dict(instance=profile['instance'],geometry=profile['geometry'])
        for warm in (False,True):
            began=time.perf_counter();profile=json.loads(path.read_text())
            d=read_instance(ROOT/'data'/profile['instance'])
            mapping={int(v):i for i,v in enumerate(d['visits'])}
            route=[mapping[int(v)] for v in profile['route']]
            row['warm' if warm else 'cold']=solve(d,route,warm,began=began)
        result['rows'].append(row)
        result['warm_improvements']=sum(r['warm']['improved'] for r in result['rows'])
        result['gate_passed']=len(result['rows'])==11 and result['warm_improvements']>=2
        out.write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps({k:v for k,v in row.items() if k not in ('warm','cold')}|{k:{a:b for a,b in row[k].items() if a!='route'} for k in ('warm','cold')}),flush=True)
