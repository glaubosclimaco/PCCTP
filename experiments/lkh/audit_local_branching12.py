"""Independent route and binary Hamming-distance audit."""
import json,hashlib
from pathlib import Path
from raw_route_audit import load,check
from audit_edge_pool import edges
HERE=Path(__file__).resolve().parent
if __name__=='__main__':
 rows=[];checked=0
 for p in sorted((HERE/'local_branching12').glob('g_*.json')):
  r=json.loads(p.read_text());d=load(r['instance']);source=HERE/'edge_pool_study'/(r['instance']+'.json');old=json.loads(source.read_text());assert hashlib.sha256(source.read_bytes()).hexdigest()==r['source_sha256']
  initial=min(old['input_routes'],key=lambda tour:sum(d['matrix'][tour[i-1]*d['n']+v] for i,v in enumerate(tour)))
  assert r['initial_cost']==old['initial_cost']
  for key in ('result','raw_result'):
   z=r[key];check(d,z['route'],z['cost']);checked+=1
   assert len(edges(initial)^edges(z['route']))+len(set(initial)^set(z['route']))<=12
   assert z['proof_scope']=='local branching neighborhood' and z['cost']<=r['initial_cost']
  late=r['native_return_seconds']>60 and r['raw_result']['cost']<r['initial_cost'];assert late==r['late_improvement_discarded']
  if late:assert r['result']['cost']==r['initial_cost']
  else:assert r['result']==r['raw_result']
  base=json.loads((HERE/'full_graph60'/(r['instance']+'.full.json')).read_text())
  rows.append(dict(instance=r['instance'],initial_cost=r['initial_cost'],cost=r['result']['cost'],full_cost=base['result']['cost'],gain_percent=100*(r['initial_cost']-r['result']['cost'])/r['initial_cost']))
 complete=(HERE/'local_branching12/completion.json').exists()
 if complete:assert len(rows)==9
 mean=sum(r['gain_percent'] for r in rows)/len(rows) if rows else None;wins=sum(r['gain_percent']>0 for r in rows);fullwins=sum(r['cost']<r['full_cost'] for r in rows)
 report=dict(status='PASS' if rows else 'NO_RESULTS',cases=len(rows),routes_checked=checked,batch_complete=complete,mean_gain_percent=mean,wins=wins,wins_over_full=fullwins,gate_passed=complete and wins>=3 and mean>=.2 and fullwins>=3,rows=rows,scope='original scalar routes, Hamming radius and discard rule; no global bound')
 (HERE/'local_branching12_audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='rows'}))
