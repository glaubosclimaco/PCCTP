"""GRASP-pool guided visit fix-and-optimize, globally free route edges."""
import sys,json,time,hashlib,random
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'experiments/ml_lns'))
from pilot import read_instance,validate
from highs_warm_start import solve
if __name__=='__main__':
 assert (HERE/'full_cuts60/completion.json').exists()
 assert json.loads((HERE/'visit_fix_verification.json').read_text())['status']=='PASS'
 out=HERE/'consensus_fix60';out.mkdir(exist_ok=True)
 for case,source in enumerate(sorted((HERE/'edge_pool_study').glob('g_*.json'))):
  for arm in (('consensus','random') if case%2==0 else ('random','consensus')):
   target=out/(source.stem+'.'+arm+'.json')
   if target.exists():continue
   began=time.perf_counter();native=time.clock_gettime(time.CLOCK_MONOTONIC);cpu=time.process_time()
   old=json.loads(source.read_text());d=read_instance(ROOT/'data'/old['instance']);mapping={int(v):i for i,v in enumerate(d['visits'])};routes=[[mapping[v] for v in route] for route in old['input_routes']]
   initial=min(routes,key=lambda r:validate(d,r));cost=validate(d,initial);optional=sorted(set(range(len(mapping)))-d['mandatory']);sets=list(map(set,routes));frequency={v:sum(v in s for s in sets)/len(sets) for v in optional};k=min(12,len(optional))
   if arm=='consensus':released=set(sorted(optional,key=lambda v:(-frequency[v]*(1-frequency[v]),v))[:k])
   else:released=set(random.Random(42).sample(optional,k))
   fixes={v:int(v in initial) for v in optional if v not in released}
   try:result=solve(d,initial,True,seconds=60.,began=began,fixed_visits=fixes)
   except RuntimeError as error:
    if str(error)!='Setup exhausted budget':raise
    result=dict(cost=cost,route=[int(d['visits'][v]) for v in initial],status=str(error),proven_optimal=False,proof_scope='fixed visit neighborhood')
   returned=time.clock_gettime(time.CLOCK_MONOTONIC)-native;raw=result;late=returned>60 and result['cost']<cost
   if late:result=dict(result,cost=cost,route=[int(d['visits'][v]) for v in initial],proven_optimal=False)
   route=[mapping[v] for v in result['route']];assert validate(d,route)==result['cost']<=cost and all(int(v in route)==val for v,val in fixes.items())
   row=dict(instance=old['instance'],arm=arm,initial_cost=cost,result=result,raw_result=raw,native_return_seconds=returned,cpu_seconds=time.process_time()-cpu,late_improvement_discarded=late,released=[int(d['visits'][v]) for v in sorted(released)],fixed_visits={str(int(d['visits'][v])):val for v,val in fixes.items()},source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),scope='offline visit fix-and-optimize; all route edges free; bounds local')
   target.write_text(json.dumps(row,indent=2)+'\n');print(json.dumps(dict(instance=old['instance'],arm=arm,cost=result['cost'])),flush=True)
 (out/'completion.json').write_text(json.dumps(dict(cases=9,arms=2))+'\n')
