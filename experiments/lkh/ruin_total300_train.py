"""Long-horizon TRAIN comparison; fresh starts, frozen 300s/10s-or-60s allocation."""
import json,time,random,sys,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from experiments.lkh.feasible_fallback import construct
from experiments.lkh.ruin_repair_absolute import run as search
from experiments.original_mdm.timed_experiment import run_one
from experiments.ml_lns.pilot import read_instance,validate
OUT=HERE/'ruin_total300_train'
def accept_post_events(cost,route,post,parent_started,budget=60.):
 improved=False
 for event in post['trace']:
  if event['native_timestamp']<=parent_started+budget and event['cost']<cost:
   cost=event['cost'];route=event['route'];improved=True
 return cost,route,improved

def run(name,seed,mode):
 assert mode in ('mdm','hybrid10','hybrid60')
 target=OUT/f'{name}.seed{seed}.{mode}.json'
 if target.exists():return json.loads(target.read_text())
 began=time.clock_gettime(time.CLOCK_MONOTONIC);cpu=time.process_time();d=read_instance(ROOT/'data'/name)
 f=construct(d);fallback=dict(cost=f['cost'],route=[int(d['visits'][v]) for v in f['route']],seconds=time.clock_gettime(time.CLOCK_MONOTONIC)-began)
 accepted=fallback['seconds']<=300;cost=fallback['cost'] if accepted else None;route=fallback['route'] if accepted else [];origin='fallback' if accepted else None
 native=run_one(name,seed,'original',300. if mode=='mdm' else float(mode[6:]),OUT/f'{name}.seed{seed}.{mode}.native.json',True,parent_started=began,hard_timeout_seconds=330.)
 if native['feasible_before_deadline'] and (cost is None or native['cost']<cost):cost=native['cost'];route=native['route'];origin='native'
 pre_cost=cost;pre_route=route[:];post=None;offset=time.clock_gettime(time.CLOCK_MONOTONIC)-began
 if mode!='mdm' and route and offset<300:
  post=search(name,'random',300-offset,initial_override=route,seed=seed)
  cost,route,improved=accept_post_events(cost,route,post,began,budget=300.)
  if improved:origin='post_event'
 elapsed=time.clock_gettime(time.CLOCK_MONOTONIC)-began
 if route:
  mapping={int(v):i for i,v in enumerate(d['visits'])};assert validate(d,[mapping[v] for v in route])==cost
 row=dict(instance=name,seed=seed,mode=mode,cost=cost,route=route,feasible=bool(route),origin=origin,native_started=began,pre_cost=pre_cost,pre_route=pre_route,fallback=fallback,native=native,post=post,post_offset=offset,wall_seconds=elapsed,cpu_seconds=time.process_time()-cpu,scope='fresh equal total300s; initialization charged; no offline pool starts')
 target.write_text(json.dumps(row,indent=2)+'\n');return row
if __name__=='__main__':
 assert (HERE/'joint_insertion_train/completion.json').exists(), 'Wait for component batch, no overlap'
 OUT.mkdir(exist_ok=True)
 cases=json.loads((ROOT/'experiments/original_mdm/hard_training_baseline/summary.json').read_text())['rows'];roles=json.loads((ROOT/'experiments/ml_lns/decoder_initialized_states/protocol.json').read_text())['roles'];assert all(roles[c['geometry']]=='train' for c in cases)
 dependencies=['experiments/original_mdm/timed_experiment.py','experiments/lkh/screen.py','experiments/ml_lns/pilot.py','experiments/ml_branching/oracle_probe.py','research_notes/ruin_long_horizon_protocol.md','experiments/lkh/audit_ruin_total300_train.py']
 protocol=dict(cases=cases,seeds=[70,71,72],budget=300,dependency_hashes={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in dependencies},source_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in ('ruin_total300_train.py','ruin_repair_absolute.py','feasible_fallback.py')})
 binaries=[ROOT/'experiments/original_mdm/build_timed_elite_bounded/original_mdm',ROOT/'experiments/original_mdm/build_timed_elite_bounded/FPMax/fim_maximal',ROOT/'experiments/original_mdm/build_timed_elite_bounded/FPMax/fpmax_hnmp',Path('/private/tmp/pcctp-lkh3-research/LKH')]
 protocol['binary_hashes']={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in binaries}
 protocol['arms']=['mdm','hybrid10','hybrid60']
 p=OUT/'protocol.json'
 if p.exists():assert json.loads(p.read_text())==protocol
 else:p.write_text(json.dumps(protocol,indent=2)+'\n')
 for c in cases:
  for seed in (70,71,72):
   order=['mdm','hybrid10','hybrid60'];random.Random(c['instance']+str(seed)).shuffle(order)
   for mode in order:
    r=run(c['instance'],seed,mode);print(json.dumps({k:r[k] for k in ('instance','seed','mode','cost','wall_seconds')}),flush=True)
 (OUT/'completion.json').write_text(json.dumps(dict(runs=81))+'\n')
