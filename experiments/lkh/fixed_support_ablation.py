"""Conditional replay ablation from exact saved pre-repair incumbent."""
import sys,json,time,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'experiments/ml_lns'))
from pilot import read_instance,validate
from fixed_support_engine import polish
if __name__=='__main__':
 assert (HERE/'ruin_total60_absolute/completion.json').exists(), 'Complete absolute batch before conditional replay'
 audit=json.loads((HERE/'ruin_total60_absolute_audit.json').read_text());assert audit['status']=='PASS' and audit['batch_complete']
 out=HERE/'fixed_support_ablation';out.mkdir(exist_ok=True)
 sources=sorted((HERE/'ruin_total60_absolute').glob('g_*.hybrid.json'));assert len(sources)==27
 for source in sources:
  target=out/source.name
  if target.exists():continue
  began=time.perf_counter();native=time.clock_gettime(time.CLOCK_MONOTONIC);old=json.loads(source.read_text());allowance=max(0.,60.-old['post_offset'])
  d=read_instance(ROOT/'data'/old['instance']);mapping={int(v):i for i,v in enumerate(d['visits'])};initial=old['pre_route'];cost=old['pre_cost']
  if initial and time.perf_counter()-began<allowance:
   result=polish(d,[mapping[v] for v in initial],began=began,seconds=allowance,seed=old['seed'])
  else:result=dict(route=initial,cost=cost,timeout=False,exit_code=None,seconds=time.perf_counter()-began)
  elapsed=time.clock_gettime(time.CLOCK_MONOTONIC)-native
  if initial:
   assert set(result['route'])==set(initial)
   assert validate(d,[mapping[v] for v in result['route']])==result['cost']<=cost
  raw=result
  if elapsed>allowance:result=dict(result,route=initial,cost=cost)
  row=dict(instance=old['instance'],seed=old['seed'],initial_route=initial,initial_cost=cost,hybrid_cost=old['cost'],control=result,raw_result=raw,phase_allowance=allowance,native_elapsed=elapsed,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),engine_sha256=hashlib.sha256((HERE/'fixed_support_engine.py').read_bytes()).hexdigest(),scope='conditional replay on identical pre-repair start; fixed support; not a new independent end-to-end run')
  target.write_text(json.dumps(row,indent=2)+'\n');print(json.dumps(dict(instance=old['instance'],seed=old['seed'],hybrid=old['cost'],fixed=result['cost'])),flush=True)
 (out/'completion.json').write_text(json.dumps(dict(runs=27))+'\n')
