"""Independent scalar route and absolute-event deadline audit."""
import json,hashlib
from pathlib import Path
from raw_route_audit import load,check
HERE=Path(__file__).resolve().parent
if __name__=='__main__':
 folder=HERE/'ruin_total60_absolute';rows=[];checked=0
 for p in sorted(folder.glob('g_*.json')):
  if '.native.' in p.name:continue
  r=json.loads(p.read_text());d=load(r['instance']);began=r['native_started'];deadline=began+60
  fallback=r['fallback'];check(d,fallback['route'],fallback['cost']);checked+=1
  candidates=[]
  if fallback['seconds']<=60:candidates.append((fallback['cost'],fallback['route']))
  native=r['native']
  if native['feasible_before_deadline']:
   check(d,native['route'],native['cost']);checked+=1;candidates.append((native['cost'],native['route']))
  pre=min(candidates,key=lambda x:x[0]) if candidates else (None,[])
  assert (r['pre_cost'],r['pre_route'])==pre
  cost,route=pre;post=r['post'];last_time=began
  if post:
   assert r['mode']=='hybrid' and post['source_sha256'] is None
   assert post['lkh_calls']==post['lkh_timeouts']+post['lkh_completed']+post['lkh_failed']
   assert post['attempts']>=post['feasible_repairs']>=post['lkh_calls']
   assert began+r['post_offset']-1e-5<=post['native_started']<=began+r['wall_seconds']
   assert post['initial_cost']==cost;previous=cost
   for event in post['trace']:
    check(d,event['route'],event['cost']);checked+=1
    timestamp=event['native_timestamp']
    assert abs(timestamp-post['native_started']-event['native_seconds'])<1e-5
    assert last_time<=timestamp<=post['native_started']+post['native_elapsed']+1e-5
    assert event['cost']<previous;previous=event['cost'];last_time=timestamp
    if timestamp<=deadline and event['cost']<cost:cost,route=event['cost'],event['route']
   check(d,post['route'],post['cost']);checked+=1
   assert post['cost']==previous
  assert (r['cost'],r['route'])==(cost,route)
  assert r['feasible']==bool(route)
  if route:check(d,route,cost);checked+=1
  rows.append(dict(instance=r['instance'],seed=r['seed'],mode=r['mode'],cost=cost,feasible=bool(route),result_sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
 complete=(folder/'completion.json').exists()
 if complete:assert len(rows)==54
 report=dict(status='PASS' if rows else 'NO_RESULTS',runs=len(rows),routes_checked=checked,batch_complete=complete,rows=rows,scope='raw scalar feasible costs and absolute monotonic-event selection; no statistical conclusion')
 (HERE/'ruin_total60_absolute_audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='rows'}))
