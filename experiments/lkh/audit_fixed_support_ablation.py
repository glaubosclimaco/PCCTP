"""Independent original-data feasibility, fixed support and same-start pairing."""
import json,hashlib
from pathlib import Path
from raw_route_audit import load,check
HERE=Path(__file__).resolve().parent
if __name__=='__main__':
 folder=HERE/'fixed_support_ablation';rows=[];checked=0
 for p in sorted(folder.glob('g_*.json')):
  r=json.loads(p.read_text());source=HERE/'ruin_total60_absolute'/p.name;old=json.loads(source.read_text());d=load(r['instance'])
  assert r['source_sha256']==hashlib.sha256(source.read_bytes()).hexdigest()
  assert r['engine_sha256']==hashlib.sha256((HERE/'fixed_support_engine.py').read_bytes()).hexdigest()
  assert r['instance']==old['instance'] and r['seed']==old['seed']
  assert r['initial_route']==old['pre_route'] and r['initial_cost']==old['pre_cost'] and r['hybrid_cost']==old['cost']
  assert abs(r['phase_allowance']-max(0.,60-old['post_offset']))<1e-8
  for key in ('control','raw_result'):
   z=r[key]
   if z['route']:
    check(d,z['route'],z['cost']);checked+=1
    assert set(z['route'])==set(r['initial_route']) and z['cost']<=r['initial_cost']
  if r['native_elapsed']>r['phase_allowance']:
   assert r['control']['cost']==r['initial_cost'] and r['control']['route']==r['initial_route']
  else:assert r['control']==r['raw_result']
  rows.append(dict(instance=r['instance'],seed=r['seed'],hybrid_cost=r['hybrid_cost'],fixed_cost=r['control']['cost']))
 complete=(folder/'completion.json').exists()
 if complete:assert len(rows)==27
 groups={}
 for r in rows:groups.setdefault(r['instance'],[]).append(r)
 geometry=[]
 for name,pairs in sorted(groups.items()):
  if {r['seed'] for r in pairs}!={42,43,44}:continue
  assert len(pairs)==3
  a=sum(r['hybrid_cost'] for r in pairs)/3;b=sum(r['fixed_cost'] for r in pairs)/3
  geometry.append(dict(instance=name,hybrid_mean=a,fixed_mean=b,hybrid_gain_percent=100*(b-a)/b))
 wins=sum(r['hybrid_gain_percent']>1e-9 for r in geometry)
 report=dict(status='PASS' if rows else 'NO_RESULTS',runs=len(rows),routes_checked=checked,batch_complete=complete,scope='same-start conditional ablation, no new independent end-to-end sample',rows=rows,geometry=geometry,component_gate_passed=complete and wins>=3,wins=wins)
 (HERE/'fixed_support_ablation_audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k not in ('rows','geometry')}))
