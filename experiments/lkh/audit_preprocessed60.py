"""Check reductions and returned tours using original scalar instance data."""
import json,hashlib
from pathlib import Path
from raw_route_audit import load,check
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
if __name__=='__main__':
 rows=[];checked=0
 for p in sorted((HERE/'preprocessed60').glob('g_*.json')):
  r=json.loads(p.read_text());d=load(r['instance']);n=d['n'];visits=[v for v,t in enumerate(d['types']) if t!=2];mandatory={v for v,t in enumerate(d['types']) if t==1};index={v:i for i,v in enumerate(visits)}
  covers=[{index[v] for v in visits if d['matrix'][v*n+w]<=d['radius']} for w,t in enumerate(d['types']) if t==2]
  original={index[v] for v in mandatory};red=r['preprocessing'];forced=set(red['mandatory']);retained=list(map(set,red['neighborhoods']))
  assert original<=forced
  total=sum(d['prizes'][v] for v in visits)
  for v in forced-original:
   assert total-d['prizes'][visits[v]]<d['quota'] or {v} in covers
  assert all(s in covers for s in retained)
  assert all(s&forced or any(t<=s for t in retained) for s in covers)
  assert red['removed_rows']==len(covers)-len(retained)
  assert r['reduction_source_sha256']==hashlib.sha256((ROOT/'experiments/preprocessing/structural.py').read_bytes()).hexdigest()
  for key in ('result','raw_result'):
   z=r[key];check(d,z['route'],z['cost']);checked+=1
   assert z['cost']<=r['initial_cost']
  late=r['native_return_seconds']>60 and r['raw_result']['cost']<r['initial_cost'];assert late==r['late_improvement_discarded']
  if late:assert r['result']['cost']==r['initial_cost']
  else:assert r['result']==r['raw_result']
  baseline=json.loads((HERE/'full_graph60'/(r['instance']+'.full.json')).read_text());assert baseline['initial_cost']==r['initial_cost']
  rows.append(dict(instance=r['instance'],cost=r['result']['cost'],baseline_cost=baseline['result']['cost'],gain_percent=100*(baseline['result']['cost']-r['result']['cost'])/baseline['result']['cost'],removed_rows=red['removed_rows'],new_mandatory=red['new_mandatory']))
 complete=(HERE/'preprocessed60/completion.json').exists()
 if complete:assert len(rows)==9
 mean=sum(r['gain_percent'] for r in rows)/len(rows) if rows else None;wins=sum(r['gain_percent']>0 for r in rows)
 report=dict(status='PASS' if rows else 'NO_RESULTS',batch_complete=complete,cases=len(rows),routes_checked=checked,mean_gain_percent=mean,wins=wins,gate_passed=complete and wins>=3 and mean>=.2,rows=rows,scope='original scalar feasibility and logical row/forcing certificates; no speed certification')
 (HERE/'preprocessed60_audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='rows'}))
