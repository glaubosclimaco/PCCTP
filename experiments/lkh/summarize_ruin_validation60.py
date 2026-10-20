"""Validation-only five-seed paired summaries; unchanged screening rule."""
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent
if __name__=='__main__':
 path=HERE/'ruin_validation60_audit.json'
 if not path.exists():print('NO_RESULTS');raise SystemExit(0)
 audit=json.loads(path.read_text())
 if audit['status']=='NO_RESULTS':print('NO_RESULTS');raise SystemExit(0)
 assert audit['status']=='PASS';protocol=json.loads((HERE/'ruin_validation60/protocol.json').read_text());names={c['instance'] for c in protocol['cases']};groups={}
 for r in audit['rows']:
  assert r['instance'] in names and r['seed'] in range(45,50) and r['mode'] in ('mdm','hybrid')
  key=(r['seed'],r['mode']);group=groups.setdefault(r['instance'],{});assert key not in group;group[key]=r
 rows=[];missing=[];failures=[]
 for name in sorted(names):
  records=groups.get(name,{})
  if not all((s,m) in records for s in range(45,50) for m in ('mdm','hybrid')):
   missing.append(dict(instance=name,runs=len(records)));continue
  if not all(r['feasible'] for r in records.values()):
   failures.append(name);continue
  mdm=sum(records[(s,'mdm')]['cost'] for s in range(45,50))/5;hybrid=sum(records[(s,'hybrid')]['cost'] for s in range(45,50))/5
  rows.append(dict(instance=name,mdm_mean=mdm,hybrid_mean=hybrid,gain_percent=100*(mdm-hybrid)/mdm))
 gains=[r['gain_percent'] for r in rows];mean=sum(gains)/len(gains) if gains else None;wins=sum(g>1e-9 for g in gains);losses=sum(g< -1e-9 for g in gains)
 complete=audit['batch_complete'];gate=complete and len(rows)==7 and not failures and mean>=.2 and wins>=3 and losses<=1
 report=dict(status='COMPLETE' if complete else 'PARTIAL',runs=audit['runs'],completed_geometries=len(rows),mean_gain_percent=mean,wins=wins,losses=losses,gate_passed=gate,rows=rows,missing=missing,feasibility_failures=failures,scope='limited reused-geometry validation screening, not untouched test; no retuning')
 (HERE/'ruin_validation60_summary.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k not in ('rows','missing')}))
