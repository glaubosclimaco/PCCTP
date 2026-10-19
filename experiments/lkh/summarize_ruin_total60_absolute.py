"""Geometry-level paired summaries and passive engine counters; no pooled seed tests."""
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent
if __name__=='__main__':
 path=HERE/'ruin_total60_absolute_audit.json';audit=json.loads(path.read_text())
 if audit['status']=='NO_RESULTS':print('NO_RESULTS');raise SystemExit(0)
 assert audit['status']=='PASS';groups={}
 for row in audit['rows']:groups.setdefault(row['instance'],{})[(row['seed'],row['mode'])]=row
 geometry=[];pending=[];engine=[]
 for name,records in sorted(groups.items()):
  hybrids=[]
  for (seed,mode),record in records.items():
   if mode!='hybrid':continue
   raw=json.loads((HERE/'ruin_total60_absolute'/f'{name}.seed{seed}.hybrid.json').read_text());post=raw['post']
   if post:hybrids.append(dict(seed=seed,start_support_size=len(raw['pre_route']),attempts=post['attempts'],calls=post['lkh_calls'],completed=post['lkh_completed'],timeouts=post['lkh_timeouts'],no_successful_exit=post['lkh_failed']))
  if hybrids:
   calls=sum(r['calls'] for r in hybrids);engine.append(dict(instance=name,runs=len(hybrids),calls=calls,completed=sum(r['completed'] for r in hybrids),timeouts=sum(r['timeouts'] for r in hybrids),no_successful_exit=sum(r['no_successful_exit'] for r in hybrids),timeout_fraction=sum(r['timeouts'] for r in hybrids)/calls if calls else None,details=hybrids))
  if not all((seed,arm) in records for seed in (42,43,44) for arm in ('mdm','hybrid')):
   pending.append(dict(instance=name,runs=len(records)));continue
  pairs=[]
  for seed in (42,43,44):
   a=records[(seed,'mdm')];b=records[(seed,'hybrid')]
   if not a['feasible'] or not b['feasible']:raise ValueError('Missing feasible pair: do not silently exclude')
   pairs.append(dict(seed=seed,mdm_cost=a['cost'],hybrid_cost=b['cost']))
  mdm=sum(r['mdm_cost'] for r in pairs)/3;hybrid=sum(r['hybrid_cost'] for r in pairs)/3
  geometry.append(dict(instance=name,mdm_mean=mdm,hybrid_mean=hybrid,gain_percent=100*(mdm-hybrid)/mdm,pairs=pairs))
 gains=[r['gain_percent'] for r in geometry];mean=sum(gains)/len(gains) if gains else None;wins=sum(g>1e-9 for g in gains);losses=sum(g< -1e-9 for g in gains)
 report=dict(status='PARTIAL' if len(gains)<9 else 'COMPLETE',batch_complete=audit['batch_complete'],runs=audit['runs'],completed_geometries=len(gains),mean_gain_percent=mean,wins=wins,losses=losses,gate_passed=len(gains)==9 and mean>=.2 and wins>=6 and losses<=1,scope='paired seed means within geometry; arithmetic relative gain, unchanged gate; no confirmatory significance',rows=geometry,pending=pending,engine=engine)
 (HERE/'ruin_total60_absolute_summary.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k not in ('rows','pending','engine')}))
