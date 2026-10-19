"""Original-data audit of incumbent/LP agreement fixings and resulting tours."""
import json,hashlib,sys
from pathlib import Path
from itertools import combinations
from raw_route_audit import load,check
from audit_edge_pool import edges
HERE=Path(__file__).resolve().parent
if __name__=='__main__':
 batch=sys.argv[1] if len(sys.argv)>1 else 'rins60'
 assert batch in ('rins60','rins60_corrected')
 rows=[];checked=0
 for p in sorted((HERE/batch).glob('g_*.json')):
  r=json.loads(p.read_text());d=load(r['instance']);source=HERE/'edge_pool_study'/(r['instance']+'.json');old=json.loads(source.read_text());assert r['source_sha256']==hashlib.sha256(source.read_bytes()).hexdigest()
  visits=[v for v,t in enumerate(d['types']) if t!=2];E=list(combinations(visits,2))
  initial=min(old['input_routes'],key=lambda tour:sum(d['matrix'][tour[i-1]*d['n']+v] for i,v in enumerate(tour)))
  def vector(tour):
   selected=set(tour);es=edges(tour);return [int(e in es) for e in E]+[int(v in selected) for v in visits]
  initial_vector=vector(initial);info=r['raw_result']['rins'];assert info['binary_count']==len(initial_vector)
  if info['lp_status']!='Optimal':assert info['fixed']==[]
  assert len({f['index'] for f in info['fixed']})==len(info['fixed'])
  for f in info['fixed']:
   assert f['value']==initial_vector[f['index']] and abs(f['relaxation_value']-f['value'])<=1e-6
  for key in ('result','raw_result'):
   z=r[key];check(d,z['route'],z['cost']);checked+=1;v=vector(z['route'])
   assert all(v[f['index']]==f['value'] for f in info['fixed'])
   assert z['cost']<=r['initial_cost']
  late=r['native_return_seconds']>60 and r['raw_result']['cost']<r['initial_cost'];assert late==r['late_improvement_discarded']
  if late:assert r['result']['cost']==r['initial_cost']
  else:assert r['result']==r['raw_result']
  baseline=json.loads((HERE/'full_graph60'/(r['instance']+'.full.json')).read_text())
  rows.append(dict(instance=r['instance'],cost=r['result']['cost'],gain_percent=100*(r['initial_cost']-r['result']['cost'])/r['initial_cost'],full_cost=baseline['result']['cost'],fixed=len(info['fixed']),lp_status=info['lp_status']))
 complete=(HERE/batch/'completion.json').exists()
 if complete:assert len(rows)==9
 mean=sum(r['gain_percent'] for r in rows)/len(rows) if rows else None;wins=sum(r['gain_percent']>0 for r in rows);fullwins=sum(r['cost']<r['full_cost'] for r in rows)
 report=dict(status='PASS' if rows else 'NO_RESULTS',cases=len(rows),routes_checked=checked,batch_complete=complete,mean_gain_percent=mean,wins=wins,wins_over_full=fullwins,gate_passed=complete and wins>=3 and mean>=.2 and fullwins>=3,rows=rows)
 (HERE/(batch+'_audit.json')).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='rows'}))
