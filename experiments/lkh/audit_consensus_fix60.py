"""Rebuild release choices and verify original-data feasibility and visit fixings."""
import json,hashlib,random
from pathlib import Path
from raw_route_audit import load,check
HERE=Path(__file__).resolve().parent
if __name__=='__main__':
 rows=[];checked=0
 for p in sorted((HERE/'consensus_fix60').glob('g_*.json')):
  r=json.loads(p.read_text());d=load(r['instance']);source=HERE/'edge_pool_study'/(r['instance']+'.json');old=json.loads(source.read_text());assert r['source_sha256']==hashlib.sha256(source.read_bytes()).hexdigest()
  tours=old['input_routes'];initial=min(tours,key=lambda tour:sum(d['matrix'][tour[i-1]*d['n']+v] for i,v in enumerate(tour)))
  optional=[v for v,t in enumerate(d['types']) if t==0];k=min(12,len(optional));freq={v:sum(v in t for t in tours)/len(tours) for v in optional}
  if r['arm']=='consensus':released=set(sorted(optional,key=lambda v:(-freq[v]*(1-freq[v]),v))[:k])
  else:
   assert r['arm']=='random';released=set(random.Random(42).sample(optional,k))
  assert r['released']==sorted(released)
  fixes={str(v):int(v in initial) for v in optional if v not in released};assert r['fixed_visits']==fixes
  for key in ('result','raw_result'):
   z=r[key];check(d,z['route'],z['cost']);checked+=1
   assert all(int(int(v) in z['route'])==value for v,value in fixes.items())
   assert z['cost']<=old['initial_cost']==r['initial_cost']
  late=r['native_return_seconds']>60 and r['raw_result']['cost']<r['initial_cost'];assert late==r['late_improvement_discarded']
  if late:assert r['result']['cost']==r['initial_cost']
  else:assert r['result']==r['raw_result']
  rows.append(dict(instance=r['instance'],arm=r['arm'],initial_cost=r['initial_cost'],cost=r['result']['cost']))
 complete=(HERE/'consensus_fix60/completion.json').exists()
 if complete:assert len(rows)==18
 report=dict(status='PASS' if rows else 'NO_RESULTS',batch_complete=complete,runs=len(rows),routes_checked=checked,rows=rows,scope='raw scalar feasibility, reconstructed release rules, fixed visits and discard rule')
 if complete:
  cons=[r for r in rows if r['arm']=='consensus'];randoms={r['instance']:r['cost'] for r in rows if r['arm']=='random'}
  gains=[100*(r['initial_cost']-r['cost'])/r['initial_cost'] for r in cons];wins=sum(g>0 for g in gains);rwins=sum(r['cost']<randoms[r['instance']] for r in cons)
  fwins=sum(r['cost']<json.loads((HERE/'full_graph60'/(r['instance']+'.full.json')).read_text())['result']['cost'] for r in cons)
  report.update(wins=wins,mean_gain_percent=sum(gains)/9,wins_over_random=rwins,wins_over_full=fwins,gate_passed=wins>=3 and sum(gains)/9>=.2 and rwins>=5 and fwins>=3)
 (HERE/'consensus_fix60_audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='rows'}))
