"""Original scalar incumbent and accepted-time audit of ruin/recreate outputs."""
import json,hashlib
from pathlib import Path
from raw_route_audit import load,check
HERE=Path(__file__).resolve().parent
if __name__=='__main__':
 rows=[];checked=0
 for p in sorted((HERE/'ruin_repair60').glob('g_*.json')):
  r=json.loads(p.read_text());d=load(r['instance']);source=HERE/'edge_pool_study'/(r['instance']+'.json');old=json.loads(source.read_text());assert r['source_sha256']==hashlib.sha256(source.read_bytes()).hexdigest()
  assert r['initial_cost']==old['initial_cost'];assert r['arm'] in ('related','random')
  check(d,r['route'],r['cost']);checked+=1;previous=r['initial_cost'];last=-1.
  for event in r['trace']:
   check(d,event['route'],event['cost']);checked+=1
   assert last<=event['native_seconds']<=60 and event['cost']<previous
   last=event['native_seconds'];previous=event['cost']
  assert r['cost']==previous and r['cost']<=r['initial_cost']
  if r['trace']:assert r['route']==r['trace'][-1]['route']
  else:assert r['route'] in old['input_routes']
  assert r['attempts']>=r['feasible_repairs']>=len(r['trace'])
  baseline=json.loads((HERE/'full_graph60'/(r['instance']+'.full.json')).read_text())
  rows.append(dict(instance=r['instance'],arm=r['arm'],cost=r['cost'],initial_cost=r['initial_cost'],full_cost=baseline['result']['cost'],gain_percent=100*(r['initial_cost']-r['cost'])/r['initial_cost'],header_target=d['header_target']))
 complete=(HERE/'ruin_repair60/completion.json').exists()
 if complete:assert len(rows)==18
 report=dict(status='PASS' if rows else 'NO_RESULTS',runs=len(rows),routes_checked=checked,batch_complete=complete,rows=rows,scope='raw feasible incumbents, decreasing accepted trace and recorded native deadline; no per-repair or CPU certification')
 if complete:
  a=[r for r in rows if r['arm']=='related'];b={r['instance']:r['cost'] for r in rows if r['arm']=='random'};wins=sum(r['gain_percent']>0 for r in a);mean=sum(r['gain_percent'] for r in a)/9;rwins=sum(r['cost']<b[r['instance']] for r in a);fwins=sum(r['cost']<r['full_cost'] for r in a)
  report.update(wins=wins,mean_gain_percent=mean,wins_over_random=rwins,wins_over_full=fwins,gate_passed=wins>=3 and mean>=.2 and rwins>=5 and fwins>=3)
 (HERE/'ruin_repair60_audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='rows'}))
