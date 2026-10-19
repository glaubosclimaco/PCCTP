"""All-vertex release, degree model and external integer connectivity cuts."""
import sys,json,time,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'experiments/ml_lns'))
from pilot import read_instance,validate
from fragment_repair import FragmentRepair
if __name__=='__main__':
 assert (HERE/'rins60_corrected/completion.json').exists()
 assert json.loads((HERE/'full_cuts_verification.json').read_text())['status']=='PASS'
 out=HERE/'full_cuts60';out.mkdir(exist_ok=True)
 for source in sorted((HERE/'edge_pool_study').glob('g_*.json')):
  target=out/source.name
  if target.exists():continue
  began=time.perf_counter();native=time.clock_gettime(time.CLOCK_MONOTONIC);cpu=time.process_time()
  old=json.loads(source.read_text());d=read_instance(ROOT/'data'/old['instance']);mapping={int(v):i for i,v in enumerate(d['visits'])}
  routes=[[mapping[v] for v in route] for route in old['input_routes']];initial=min(routes,key=lambda r:validate(d,r));cost=validate(d,initial)
  model=FragmentRepair(d,'highspy',True,'cuts');remaining=60-(time.perf_counter()-began)
  if remaining<=0:raise RuntimeError('Setup exhausted budget')
  result=model.solve(initial,set(range(len(mapping))),remaining,'incident');returned=time.clock_gettime(time.CLOCK_MONOTONIC)-native;raw=result
  assert result['edge_variables']==len(mapping)*(len(mapping)-1)//2 and result['fixed_edge_constant']==0
  late=returned>60 and result['cost']<cost
  if late:result=dict(result,cost=cost,route=[int(d['visits'][v]) for v in initial],proven_optimal=False)
  validate(d,[mapping[v] for v in result['route']]);assert result['cost']<=cost
  r=dict(instance=old['instance'],initial_cost=cost,result=result,raw_result=raw,native_return_seconds=returned,cpu_seconds=time.process_time()-cpu,late_improvement_discarded=late,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),scope='offline full graph integer outer cuts; input generation excluded; no lazy callback claim')
  target.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(dict(instance=r['instance'],cost=result['cost'],cuts=result['cuts_added'])),flush=True)
 (out/'completion.json').write_text(json.dumps(dict(cases=9))+'\n')
