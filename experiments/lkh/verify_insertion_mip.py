"""Independent assignment enumeration, including negative costs and restricted infeasibility."""
import json,random
from pathlib import Path
from itertools import product
import numpy as np
from insertion_mip_repair import repair
from pilot import validate
if __name__=='__main__':
 rng=random.Random(20261019);rows=[];negative=0;assignments_checked=0
 for trial in range(20):
  n=6;c=np.zeros((n,n))
  for i in range(n):
   for j in range(i):c[i,j]=c[j,i]=rng.randrange(1,100)
  prizes=np.array([rng.randrange(1,20) for _ in range(n)]);d=dict(visits=np.arange(n),mandatory={0},costs=c,prizes=prizes,quota=rng.randrange(1,int(sum(prizes))+1),neighborhoods=[{1,3},{4,5}]);base=[0,1,2];slots=None if trial%2==0 else {3:[0],4:[1],5:[2]};values=[]
  for positions in product((-1,0,1,2),repeat=3):
   assignments_checked+=1;used=[i for i in positions if i>=0]
   if len(set(used))!=len(used):continue
   if slots is not None and any(i>=0 and i not in slots[v] for v,i in zip((3,4,5),positions)):continue
   at={i:v for v,i in zip((3,4,5),positions) if i>=0};route=[]
   for i,v in enumerate(base):
    if i in at:route.append(at[i])
    route.append(v)
   try:value=validate(d,route)
   except ValueError:continue
   values.append(value)
  r=repair(d,base,seconds=2.,slots=slots)
  if values:assert r['proven_local_optimal'] and r['cost']==min(values)
  else:assert r['cost'] is None and 'Infeasible' in r['status']
  negative+=sum(c[base[i-1],v]+c[v,base[i]]-c[base[i-1],base[i]]<0 for v in (3,4,5) for i in range(3))
  rows.append(dict(trial=trial,optimum=min(values) if values else None,solver_cost=r['cost'],restricted_slots=slots is not None))
 # Original full cycle feasible, but four required additions cannot fit three gaps.
 n=7;d=dict(visits=np.arange(n),mandatory={0},costs=np.ones((n,n))-np.eye(n),prizes=np.ones(n),quota=7,neighborhoods=[{3,4}]);validate(d,list(range(n)));r=repair(d,[0,1,2],2.);assert r['cost'] is None and 'Infeasible' in r['status']
 assert negative>0
 report=dict(status='PASS',tiny_models=21,assignments_checked=assignments_checked,negative_insertion_coefficients=int(negative),scope='finite exact assignment checks; restricted infeasibility explicitly differs from original feasibility',rows=rows)
 Path(__file__).with_name('insertion_mip_verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='rows'}))
