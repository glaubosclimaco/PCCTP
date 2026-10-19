"""Randomized coverage/prize repair; permits temporary infeasibility after removal."""
import sys,time,random,json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
import numpy as np
from experiments.ml_lns.pilot import read_instance,validate
from experiments.lkh.screen import polish

def repair(d,route,rng,end):
 route=route[:];n=len(d['visits']);mean=max(1.,float(np.mean(d['prizes'])))
 while time.perf_counter()<end:
  selected=set(route);uncovered=[s for s in d['neighborhoods'] if not s&selected];deficit=max(0.,d['quota']-sum(d['prizes'][list(selected)]))
  if not uncovered and deficit<=0 and len(route)>=3:return route
  choices=[]
  for v in sorted(set(range(n))-selected):
   gain=sum(v in s for s in uncovered)+min(float(d['prizes'][v]),deficit)/mean
   if len(route)<3:gain+=1
   if gain<=0:continue
   before=np.roll(route,1);delta=d['costs'][before,v]+d['costs'][v,route]-d['costs'][before,route];pos=int(np.argmin(delta))
   score=float(delta[pos])/gain*rng.uniform(.8,1.2);choices.append((score,v,pos))
  if not choices:return None
  _,v,pos=min(choices);route.insert(pos,v)
 return None

def run(name,arm,budget=60.,initial_override=None,seed=42):
 began=time.perf_counter();native=time.clock_gettime(time.CLOCK_MONOTONIC);cpu=time.process_time();end=began+budget
 source=HERE/'edge_pool_study'/(name+'.json');old=json.loads(source.read_text()) if initial_override is None else dict(input_routes=[]);d=read_instance(ROOT/'data'/name);mapping={int(v):i for i,v in enumerate(d['visits'])};routes=[[mapping[v] for v in r] for r in old['input_routes']];best=min(routes,key=lambda r:validate(d,r)) if initial_override is None else [mapping[v] for v in initial_override];initial=validate(d,best);cost=initial;current=best[:];rng=random.Random(seed);trace=[];attempts=0;feasible=0;lkh_calls=0;lkh_timeouts=0;lkh_completed=0;lkh_failed=0
 while time.perf_counter()<end and time.clock_gettime(time.CLOCK_MONOTONIC)-native<budget:
  optional=sorted(set(current)-d['mandatory']);k=min(6,len(optional),max(0,len(current)-3))
  if k==0:break
  if arm=='related':
   seed=rng.choice(optional);removed=set(sorted(optional,key=lambda v:(d['costs'][seed,v],v))[:k])
  else:removed=set(rng.sample(optional,k))
  candidate=repair(d,[v for v in current if v not in removed],rng,end);attempts+=1
  if candidate is None:break
  validate(d,candidate);feasible+=1;remaining=end-time.perf_counter()
  if remaining<=0:break
  result=polish(d,candidate,seconds=min(.2,remaining));lkh_calls+=1
  if result['timeout']:lkh_timeouts+=1
  elif result['exit_code']==0:lkh_completed+=1
  else:lkh_failed+=1
  candidate=[mapping[v] for v in result['route']];value=validate(d,candidate);elapsed=time.clock_gettime(time.CLOCK_MONOTONIC)-native
  if elapsed>budget or time.perf_counter()>end:break
  if value<cost:
   best=candidate[:];cost=value;trace.append(dict(native_seconds=elapsed,native_timestamp=native+elapsed,cost=cost,route=result['route']))
  if value<=cost*1.01:current=candidate
  if attempts%20==0:current=best[:]
 return dict(instance=name,arm=arm,initial_cost=initial,cost=cost,route=[int(d['visits'][v]) for v in best],trace=trace,attempts=attempts,feasible_repairs=feasible,lkh_calls=lkh_calls,lkh_timeouts=lkh_timeouts,lkh_completed=lkh_completed,lkh_failed=lkh_failed,native_started=native,native_elapsed=time.clock_gettime(time.CLOCK_MONOTONIC)-native,cpu_seconds=time.process_time()-cpu,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest() if initial_override is None else None,scope='fresh supplied start, pool unread' if initial_override is not None else 'offline pool start; generation excluded')
if __name__=='__main__':
 assert (HERE/'consensus_fix60/completion.json').exists()
 out=HERE/'ruin_repair60';out.mkdir(exist_ok=True)
 for i,source in enumerate(sorted((HERE/'edge_pool_study').glob('g_*.json'))):
  for arm in (('related','random') if i%2==0 else ('random','related')):
   target=out/(source.stem+'.'+arm+'.json')
   if target.exists():continue
   r=run(source.stem,arm);target.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:r[k] for k in ('instance','arm','cost','attempts')}),flush=True)
 (out/'completion.json').write_text(json.dumps(dict(cases=9,arms=2))+'\n')
