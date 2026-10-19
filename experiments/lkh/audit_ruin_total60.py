"""Audit wrapper outputs and separately bound predeadline postprocessing events."""
import json
from pathlib import Path
from raw_route_audit import load,check
HERE=Path(__file__).resolve().parent
if __name__=='__main__':
 rows=[];checked=0
 folder=HERE/'ruin_total60'
 for p in sorted(folder.glob('g_*.json')):
  if '.native.' in p.name:continue
  r=json.loads(p.read_text());d=load(r['instance'])
  for route,cost in ((r['route'],r['cost']),(r['pre_route'],r['pre_cost']),(r['fallback']['route'],r['fallback']['cost'])):
   if route:check(d,route,cost);checked+=1
  assert r['feasible']==bool(r['route'])
  if r['native']['feasible_before_deadline']:check(d,r['native']['route'],r['native']['cost']);checked+=1
  if r['origin']=='post':assert r['wall_seconds']<=60 and r['cost']<r['pre_cost']
  else:assert r['cost']==r['pre_cost']
  anytime=r['pre_cost'];accepted=[];post=r['post']
  if post:
   check(d,post['route'],post['cost']);checked+=1
   assert post['initial_cost']==r['pre_cost']
   # Parent elapsed at or after post return minus post's measured duration
   # is an upper bound on post start; all use CLOCK_MONOTONIC durations.
   start_upper=r['wall_seconds']-post['native_elapsed']
   assert start_upper>=r['post_offset']-1e-5
   previous=post['initial_cost']
   for e in post['trace']:
    check(d,e['route'],e['cost']);checked+=1
    assert e['cost']<previous;previous=e['cost']
    upper=start_upper+e['native_seconds']
    if upper<=60:
     anytime=min(anytime,e['cost']);accepted.append(dict(cost=e['cost'],parent_time_upper=upper))
   assert post['cost']==previous
  rows.append(dict(instance=r['instance'],seed=r['seed'],mode=r['mode'],wrapper_cost=r['cost'],conservative_anytime_cost=anytime,pre_cost=r['pre_cost'],bounded_events=accepted))
 report=dict(status='PASS' if rows else 'NO_RESULTS',runs=len(rows),routes_checked=checked,batch_complete=(folder/'completion.json').exists(),scope='original scalar routes; wrapper costs retained; separate conservative parent-clock event reconstruction, not replacement of preregistered outcome',rows=rows)
 (HERE/'ruin_total60_audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='rows'}))
