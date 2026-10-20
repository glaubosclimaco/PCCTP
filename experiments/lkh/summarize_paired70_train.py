"""Prespecified TRAIN screen; incomplete geometries never enter the decision."""
import hashlib,json
from pathlib import Path
HERE=Path(__file__).resolve().parent

def summarize(audit,protocol):
    assert audit['status']=='PASS'
    rows=audit['rows'];geometries=[]
    for case in protocol['cases']:
        matched=[r for r in rows if r['instance']==case['instance']]
        complete=(len(matched)==10 and {r['seed'] for r in matched}==set(range(81,91))
                  and all(r['status']=='PAIR_COMPLETED' for r in matched))
        record=dict(instance=case['instance'],geometry=case['geometry'],pairs=len(matched),complete=complete)
        if complete:
            baseline=sum(r['reference_cost'] for r in matched)/10
            record.update(reference_mean=baseline,reference_mean_seconds=sum(r['duration'] for r in matched)/10,
                candidates={mode:dict(mean=sum(r['candidate_costs'][mode] for r in matched)/10,
                    gain_percent=100*(baseline-sum(r['candidate_costs'][mode] for r in matched)/10)/baseline,
                    seed_wins=sum(r['candidate_costs'][mode]<r['reference_cost'] for r in matched),
                    seed_losses=sum(r['candidate_costs'][mode]>r['reference_cost'] for r in matched))
                    for mode in ('hybrid10','hybrid60')})
        geometries.append(record)
    completed=[g for g in geometries if g['complete']]
    eligible=audit['eligible_complete_study'] and len(completed)==9
    comparisons={}
    for mode in ('hybrid10','hybrid60'):
        gains=[g['candidates'][mode]['gain_percent'] for g in completed]
        mean=sum(gains)/len(gains) if gains else None
        wins=sum(x>1e-12 for x in gains);losses=sum(x< -1e-12 for x in gains)
        comparisons[mode]=dict(mean_gain_percent=mean,wins=wins,losses=losses,
            ties=len(gains)-wins-losses,gate_passed=eligible and mean>=.2 and wins>=6 and losses<=1)
    return dict(status='COMPLETE' if audit['batch_complete'] else 'PARTIAL',pairs=audit['pairs'],
        planned_pairs=90,completed_geometries=len(completed),eligible_complete_study=eligible,
        reference_incomplete=audit['reference_incomplete'],candidate_failure_pairs=audit['candidate_failure_pairs'],
        comparisons=comparisons,geometries=geometries,
        scope='descriptive TRAIN screen, equal reference-defined per-pair wall budget; seeds not independent geometries; no test/107-instance or publication verdict')

if __name__=='__main__':
    audit_path=HERE/'paired70_train_audit.json';protocol_path=HERE/'paired70_train/protocol.json'
    if not audit_path.exists() or not protocol_path.exists():print('NO_RESULTS');raise SystemExit(0)
    report=summarize(json.loads(audit_path.read_text()),json.loads(protocol_path.read_text()))
    report['audit_sha256']=hashlib.sha256(audit_path.read_bytes()).hexdigest()
    report['protocol_sha256']=hashlib.sha256(protocol_path.read_bytes()).hexdigest()
    (HERE/'paired70_train_summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='geometries'}))
