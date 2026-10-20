"""Both prespecified long-horizon TRAIN comparisons; geometry means, no inference."""
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent
if __name__=='__main__':
 path=HERE/'ruin_total300_train_audit.json'
 if not path.exists():print('NO_RESULTS');raise SystemExit(0)
 audit=json.loads(path.read_text());assert audit['status']=='PASS'
 groups={}
 for r in audit['rows']:groups.setdefault(r['instance'],{})[r['seed'],r['mode']]=r
 rows=[];failures=[]
 for name,records in sorted(groups.items()):
  if len(records)!=9:continue
  if not all(r['feasible'] for r in records.values()):failures.append(name);continue
  means={m:sum(records[s,m]['cost'] for s in (70,71,72))/3 for m in ('mdm','hybrid10','hybrid60')}
  rows.append(dict(instance=name,means=means,
   gains={m:100*(means['mdm']-means[m])/means['mdm'] for m in ('hybrid10','hybrid60')},
   hybrid60_vs_10_percent=100*(means['hybrid10']-means['hybrid60'])/means['hybrid10']))
 comparisons={}
 for m in ('hybrid10','hybrid60'):
  gains=[r['gains'][m] for r in rows];mean=sum(gains)/len(gains) if gains else None
  wins=sum(g>1e-9 for g in gains);losses=sum(g< -1e-9 for g in gains)
  comparisons[m]=dict(mean_gain_percent=mean,wins=wins,losses=losses,
   gate_passed=audit['batch_complete'] and len(rows)==9 and not failures and mean>=.2 and wins>=6 and losses<=1)
 report=dict(status='COMPLETE' if audit['batch_complete'] else 'PARTIAL',runs=audit['runs'],
  completed_geometries=len(rows),comparisons=comparisons,rows=rows,feasibility_failures=failures,
  scope='exploratory TRAIN at300seconds; both prespecified arms retained; no test or inferential conclusion')
 (HERE/'ruin_total300_train_summary.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({k:v for k,v in report.items() if k!='rows'}))
