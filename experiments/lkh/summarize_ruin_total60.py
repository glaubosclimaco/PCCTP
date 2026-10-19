"""Paired geometry summaries; original and reconstructed outcomes stay separate."""
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent
if __name__=='__main__':
 audit=json.loads((HERE/'ruin_total60_audit.json').read_text());assert audit['status']=='PASS'
 groups={}
 for r in audit['rows']:groups.setdefault(r['instance'],{})[(r['seed'],r['mode'])]=r
 complete=[];pending=[]
 for name,records in sorted(groups.items()):
  if not all((seed,arm) in records for seed in (42,43,44) for arm in ('mdm','hybrid')):
   pending.append(dict(instance=name,runs=len(records)));continue
  baseline=sum(records[(s,'mdm')]['wrapper_cost'] for s in (42,43,44))/3
  strict=sum(records[(s,'hybrid')]['wrapper_cost'] for s in (42,43,44))/3
  conservative=sum(records[(s,'hybrid')]['conservative_anytime_cost'] for s in (42,43,44))/3
  complete.append(dict(instance=name,mdm_mean=baseline,wrapper_hybrid_mean=strict,conservative_hybrid_mean=conservative,wrapper_gain_percent=100*(baseline-strict)/baseline,conservative_gain_percent=100*(baseline-conservative)/baseline))
 summaries={}
 for outcome in ('wrapper','conservative'):
  gains=[r[outcome+'_gain_percent'] for r in complete];mean=sum(gains)/len(gains) if gains else None;wins=sum(g>1e-9 for g in gains);losses=sum(g< -1e-9 for g in gains)
  summaries[outcome]=dict(completed_geometries=len(gains),mean_gain_percent=mean,wins=wins,losses=losses,gate_passed=(len(gains)==9 and mean>=.2 and wins>=6 and losses<=1))
 report=dict(batch_complete=audit['batch_complete'],runs=audit['runs'],scope='complete seed triplets only, paired within geometry; conservative reconstruction is separate sensitivity outcome, no significance or held-out claim',summaries=summaries,rows=complete,pending=pending)
 (HERE/'ruin_total60_summary.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k not in ('rows','pending')}))
