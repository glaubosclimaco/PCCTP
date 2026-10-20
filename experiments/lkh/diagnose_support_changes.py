"""TRAIN-only descriptive feasibility witnesses for net support changes."""
import json,hashlib
from pathlib import Path
from raw_route_audit import load,check
HERE=Path(__file__).resolve().parent
if __name__=='__main__':
 rows=[]
 for p in sorted((HERE/'ruin_total60_absolute').glob('g_*.hybrid.json')):
  r=json.loads(p.read_text());d=load(r['instance']);n=d['n'];old=set(r['pre_route']);new=set(r['route']);removed=sorted(old-new);added=sorted(new-old)
  check(d,r['pre_route'],r['pre_cost']);check(d,r['route'],r['cost'])
  cover=[{v for v,t in enumerate(d['types']) if t!=2 and d['matrix'][v*n+w]<=d['radius']} for w,t in enumerate(d['types']) if t==2]
  mandatory={v for v,t in enumerate(d['types']) if t==1};prize_before=sum(d['prizes'][v] for v in old);prize_after=sum(d['prizes'][v] for v in new)
  details=[]
  for v in removed:
   assert v not in mandatory
   coverage_blocked=any(s&old=={v} for s in cover);prize_blocked=prize_before-d['prizes'][v]<d['quota']
   details.append(dict(vertex=v,coverage_blocked=coverage_blocked,prize_blocked=prize_blocked))
  reduced=old-set(removed);collective_cover=any(not s&reduced for s in cover);collective_prize=sum(d['prizes'][v] for v in reduced)<d['quota']
  rows.append(dict(instance=r['instance'],seed=r['seed'],initial_cost=r['pre_cost'],final_cost=r['cost'],removed=removed,added=added,removed_details=details,coverage_rows_not_satisfied_by_mandatory=sum(not s&mandatory for s in cover),prize_slack_before=prize_before-d['quota'],prize_slack_after=prize_after-d['quota'],net_removal_without_additions_violates_coverage=collective_cover,net_removal_without_additions_violates_prize=collective_prize,source_sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
 assert len(rows)==27
 report=dict(scope='TRAIN net-support feasibility witnesses only; does not reconstruct actual destruction path or establish causality',runs=27,support_changed_runs=sum(bool(r['removed'] or r['added']) for r in rows),runs_with_net_coverage_repair_need=sum(r['net_removal_without_additions_violates_coverage'] for r in rows),runs_with_net_prize_repair_need=sum(r['net_removal_without_additions_violates_prize'] for r in rows),runs_with_any_net_repair_need=sum(r['net_removal_without_additions_violates_coverage'] or r['net_removal_without_additions_violates_prize'] for r in rows),rows=rows)
 (HERE/'support_change_diagnosis.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='rows'}))
