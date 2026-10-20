"""Joint insertion repair, exact only within one-new-vertex-per-retained-gap model.

Internal visit indices throughout. Not integrated into frozen validation search.
A returned route must satisfy `accepted` before a caller uses it under a deadline.
"""
import sys,time
from pathlib import Path
import numpy as np
from scipy.optimize import Bounds,LinearConstraint
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'experiments/ml_lns'))
from pilot import validate
from experiments.ml_branching.oracle_probe import sparse_rows
from highs_backend import solve as mip

def repair(d,retained,seconds=1.,slots=None):
 began=time.perf_counter();native=time.clock_gettime(time.CLOCK_MONOTONIC)
 base=list(retained);S=set(base);n=len(d['visits']);assert len(base)==len(S)>=3 and d['mandatory']<=S
 absent=sorted(set(range(n))-S)
 options=[(v,i) for v in absent for i in (range(len(base)) if slots is None else sorted(set(slots.get(v,[]))))]
 assert all(0<=i<len(base) for v,i in options)
 cost0=float(sum(d['costs'][base[i-1],v] for i,v in enumerate(base)))
 if not options:
  try:cost=validate(d,base)
  except ValueError:return dict(status='Infeasible restricted insertion model',route=None,cost=None,accepted=False,proof_scope='restricted insertion neighborhood')
  timely=time.clock_gettime(time.CLOCK_MONOTONIC)-native<=seconds
  return dict(status='Optimal',route=base,cost=cost,accepted=timely,proof_scope='restricted insertion neighborhood',proven_local_optimal=timely,assignments=[],variables=0)
 obj=np.array([d['costs'][base[i-1],v]+d['costs'][v,base[i]]-d['costs'][base[i-1],base[i]] for v,i in options]);rows=[];lo=[];hi=[]
 def add(row,lower=-np.inf,upper=np.inf):rows.append(row);lo.append(lower);hi.append(upper)
 for v in absent:add({j:1. for j,(w,i) in enumerate(options) if w==v},upper=1.)
 for i in range(len(base)):add({j:1. for j,(v,k) in enumerate(options) if k==i},upper=1.)
 for cover in d['neighborhoods']:
  if not cover&S:add({j:1. for j,(v,i) in enumerate(options) if v in cover},lower=1.)
 deficit=max(0.,d['quota']-sum(d['prizes'][list(S)]));add({j:float(d['prizes'][v]) for j,(v,i) in enumerate(options)},lower=deficit)
 remaining=seconds-(time.perf_counter()-began)
 if remaining<=0:return dict(status='Setup exhausted budget',route=None,cost=None,accepted=False,proof_scope='restricted insertion neighborhood')
 result=mip(obj,np.ones(len(options)),Bounds(np.zeros(len(options)),np.ones(len(options))),LinearConstraint(sparse_rows(rows,len(options)),lo,hi),remaining)
 candidate=None;cost=None;assignments=[]
 if result.x is not None:
  assert np.max(np.abs(result.x-np.rint(result.x)))<1e-5
  assignments=[options[j] for j,x in enumerate(result.x) if x>.5];assert len({v for v,i in assignments})==len(assignments)==len({i for v,i in assignments})
  at={i:v for v,i in assignments};candidate=[]
  for i,v in enumerate(base):
   if i in at:candidate.append(at[i])
   candidate.append(v)
  cost=validate(d,candidate);assert abs(cost-cost0-result.fun)<1e-5
 elapsed=time.clock_gettime(time.CLOCK_MONOTONIC)-native;accepted=candidate is not None and elapsed<=seconds
 return dict(status=result.message,route=candidate,cost=cost,accepted=accepted,proof_scope='restricted insertion neighborhood',proven_local_optimal=result.status==0 and accepted,assignments=assignments,variables=len(options),constraints=len(rows),native_elapsed=elapsed)
