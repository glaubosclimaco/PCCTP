"""Fresh equal-total-time TRAIN comparison; never seeds from offline route pool."""
import json,time,random,sys,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from experiments.lkh.feasible_fallback import construct
from experiments.lkh.ruin_repair import run as search
from experiments.original_mdm.timed_experiment import run_one
from experiments.ml_lns.pilot import read_instance,validate
OUT=HERE/'ruin_total60'
def run(name,seed,mode):
 target=OUT/f'{name}.seed{seed}.{mode}.json'
 if target.exists():return json.loads(target.read_text())
 began=time.clock_gettime(time.CLOCK_MONOTONIC);cpu=time.process_time();d=read_instance(ROOT/'data'/name)
 f=construct(d);fallback=dict(cost=f['cost'],route=[int(d['visits'][v]) for v in f['route']],seconds=time.clock_gettime(time.CLOCK_MONOTONIC)-began)
 accepted=fallback['seconds']<=60;cost=fallback['cost'] if accepted else None;route=fallback['route'] if accepted else [];origin='fallback' if accepted else None
 native=run_one(name,seed,'original',60. if mode=='mdm' else 10.,OUT/f'{name}.seed{seed}.{mode}.native.json',True,parent_started=began,hard_timeout_seconds=90.)
 if native['feasible_before_deadline'] and (cost is None or native['cost']<cost):cost=native['cost'];route=native['route'];origin='native'
 pre_cost=cost;pre_route=route[:];post=None;offset=time.clock_gettime(time.CLOCK_MONOTONIC)-began
 if mode=='hybrid' and route and offset<60:
  post=search(name,'random',60-offset,initial_override=route,seed=seed)
  if time.clock_gettime(time.CLOCK_MONOTONIC)-began<=60 and post['cost']<cost:cost=post['cost'];route=post['route'];origin='post'
 elapsed=time.clock_gettime(time.CLOCK_MONOTONIC)-began
 if route:
  mapping={int(v):i for i,v in enumerate(d['visits'])};assert validate(d,[mapping[v] for v in route])==cost
 row=dict(instance=name,seed=seed,mode=mode,cost=cost,route=route,feasible=bool(route),origin=origin,pre_cost=pre_cost,pre_route=pre_route,fallback=fallback,native=native,post=post,post_offset=offset,wall_seconds=elapsed,cpu_seconds=time.process_time()-cpu,scope='fresh equal total60s; initialization charged; no offline pool starts')
 target.write_text(json.dumps(row,indent=2)+'\n');return row
if __name__=='__main__':
 assert (HERE/'ruin_repair60/completion.json').exists();OUT.mkdir(exist_ok=True)
 cases=json.loads((ROOT/'experiments/original_mdm/hard_training_baseline/summary.json').read_text())['rows'];roles=json.loads((ROOT/'experiments/ml_lns/decoder_initialized_states/protocol.json').read_text())['roles'];assert all(roles[c['geometry']]=='train' for c in cases)
 protocol=dict(cases=cases,seeds=[42,43,44],budget=60,source_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in ('ruin_total60.py','ruin_repair.py','feasible_fallback.py')})
 p=OUT/'protocol.json'
 if p.exists():assert json.loads(p.read_text())==protocol
 else:p.write_text(json.dumps(protocol,indent=2)+'\n')
 for c in cases:
  for seed in (42,43,44):
   order=['mdm','hybrid'];random.Random(c['instance']+str(seed)).shuffle(order)
   for mode in order:
    r=run(c['instance'],seed,mode);print(json.dumps({k:r[k] for k in ('instance','seed','mode','cost','wall_seconds')}),flush=True)
 (OUT/'completion.json').write_text(json.dumps(dict(runs=54))+'\n')
