"""Fixed-support DFJ diagnostic; frozen rule in route_polishing_protocol.md."""
import time
STARTUP=time.perf_counter()
import json
from pathlib import Path
from itertools import combinations
import sys
import numpy as np
from scipy.optimize import milp, Bounds, LinearConstraint
from scipy.sparse import csr_matrix

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from experiments.ml_lns.pilot import read_instance, validate
LIMIT=10.


def polish(path, variable_visits=False):
    began=time.perf_counter()
    profile=json.loads(path.read_text())
    d=read_instance(ROOT/'data'/profile['instance'])
    mapping={int(v):i for i,v in enumerate(d['visits'])}
    route=[mapping[int(v)] for v in profile['route']]
    old=validate(d,route)
    assert abs(old-profile['cost'])<1e-6
    chosen=list(range(len(d['visits']))) if variable_visits else sorted(route)
    edges=list(combinations(chosen,2)); m=len(edges)
    width=m+len(chosen) if variable_visits else m
    costs=np.r_[[d['costs'][a,b] for a,b in edges],np.zeros(width-m)]
    lb=np.zeros(width);ub=np.ones(width)
    rows=[]; lower=[]; upper=[]
    for v in chosen:
        row={i:1. for i,e in enumerate(edges) if v in e}
        if variable_visits:row[m+v]=-2.
        rows.append(row)
        lower.append(0. if variable_visits else 2.);upper.append(lower[-1])
    if variable_visits:
        for v in d['mandatory']:lb[m+v]=1.
        for ns in d['neighborhoods']:
            rows.append({m+v:1. for v in ns});lower.append(1.);upper.append(np.inf)
        rows.append({m+v:float(d['prizes'][v]) for v in chosen})
        lower.append(float(d['quota']));upper.append(np.inf)
        for i,(a,b) in enumerate(edges):
            for v in (a,b):
                rows.append({i:1.,m+v:-1.});lower.append(-np.inf);upper.append(0.)
    rows.append(dict(enumerate(costs)));lower.append(-np.inf);upper.append(old+1e-7)
    best=route; cuts=0; rounds=0; proven=False; late=False; status=None
    while True:
        remaining=LIMIT-(time.perf_counter()-began)
        if remaining<=0:break
        rr=[];cc=[];vv=[]
        for r,row in enumerate(rows):
            for c,v in row.items():rr.append(r);cc.append(c);vv.append(v)
        matrix=csr_matrix((vv,(rr,cc)),shape=(len(rows),width))
        remaining=LIMIT-(time.perf_counter()-began)
        if remaining<=0:break
        result=milp(costs,integrality=np.ones(width),bounds=Bounds(lb,ub),
                    constraints=LinearConstraint(matrix,lower,upper),
                    options={'time_limit':remaining,'mip_rel_gap':0.})
        rounds+=1;status=int(result.status)
        if time.perf_counter()-began>LIMIT:
            late=True;break
        if result.x is None:break
        if max(abs(result.x-np.rint(result.x)))>1e-5:raise ValueError('Fractional incumbent')
        selected={v for v in chosen if not variable_visits or result.x[m+v]>.5}
        adj={v:[] for v in selected}
        for i,(a,b) in enumerate(edges):
            if result.x[i]>.5:adj[a].append(b);adj[b].append(a)
        assert all(len(ns)==2 for ns in adj.values())
        unseen=set(selected);components=[]
        while unseen:
            pending=[min(unseen)];seen=set()
            while pending:
                v=pending.pop()
                if v in seen:continue
                seen.add(v);pending.extend(adj[v])
            unseen-=seen;components.append(seen)
        if len(components)==1:
            tour=[min(d['mandatory'])];previous=None
            while True:
                nxt=next(v for v in adj[tour[-1]] if v!=previous)
                if nxt==tour[0]:break
                previous=tour[-1];tour.append(nxt)
            assert set(tour)==selected
            if not variable_visits:assert selected==set(route)
            cost=validate(d,tour)
            assert abs(cost-result.fun)<1e-4
            if time.perf_counter()-began>LIMIT:late=True;break
            if cost<validate(d,best)-1e-7:best=tour
            proven=status==0
            break
        for component in components:
            if variable_visits and min(d['mandatory']) in component:continue
            row={i:1. for i,(a,b) in enumerate(edges) if (a in component)!=(b in component)}
            for v in sorted(component) if variable_visits else [None]:
                cut=dict(row)
                if variable_visits:cut[m+v]=-2.
                rows.append(cut);lower.append(0. if variable_visits else 2.)
                upper.append(np.inf);cuts+=1
    elapsed=time.perf_counter()-began
    return dict(instance=profile['instance'],geometry=profile['geometry'],mdm_cost=old,
                polished_cost=validate(d,best),improved=validate(d,best)<old-1e-7,
                initial_support_size=len(route),variable_visits=variable_visits,
                proven_optimal=proven,
                seconds=elapsed,late_result_discarded=late,status=status,
                rounds=rounds,cuts=cuts,route=[int(d['visits'][v]) for v in best])


if __name__=='__main__':
    directory=ROOT/'experiments/original_mdm/move_profiles_monotonic_elite_bounded'
    variable_visits='--variable-visits' in sys.argv
    output=Path(__file__).with_name('full_dfj_results.json' if variable_visits else 'route_polishing_results.json')
    result=dict(import_seconds=time.perf_counter()-STARTUP,limit_seconds=LIMIT,
                rows=[],excluded_incomplete=[])
    for path in sorted(directory.glob('*.profiled_baseline.json')):
        p=json.loads(path.read_text())
        if not p['completed'] or p['iterations']!=70:
            result['excluded_incomplete'].append(p['instance']);continue
        row=polish(path,variable_visits);result['rows'].append(row)
        result['improvements']=sum(r['improved'] for r in result['rows'])
        result['pilot_gate_passed']=len(result['rows'])==11 and result['improvements']>=2
        output.write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps({k:v for k,v in row.items() if k!='route'}),flush=True)
