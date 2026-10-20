"""Prespecified exploratory TRAIN gate; no inference or reserved evaluation."""
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
if __name__=='__main__':
    path=HERE/'joint_insertion_train_audit.json'
    if not path.exists():
        print('NO_RESULTS');raise SystemExit(0)
    audit=json.loads(path.read_text())
    assert audit['status']=='PASS'
    groups={}
    for r in audit['rows']:
        group=groups.setdefault(r['instance'],{})
        key=(r['seed'],r['sample'],r['arm'])
        assert key not in group
        group[key]=r
    rows=[]
    arms=('greedy','mip_nearest','mip_random')
    for name,group in sorted(groups.items()):
        if len(group)!=36:continue
        means={a:sum(group[s,k,a]['cost'] for s in (42,43,44) for k in range(4))/12
               for a in arms}
        rows.append(dict(instance=name,means=means,
            gain_vs_greedy=100*(means['greedy']-means['mip_nearest'])/means['greedy'],
            gain_vs_random=100*(means['mip_random']-means['mip_nearest'])/means['mip_random']))
    gains=[r['gain_vs_greedy'] for r in rows]
    mean=sum(gains)/len(gains) if gains else None
    wins=sum(g>1e-9 for g in gains)
    random_wins=sum(r['gain_vs_random']>1e-9 for r in rows)
    complete=audit['batch_complete']
    gate=complete and len(rows)==9 and mean>=.2 and wins>=3 and random_wins>=3
    report=dict(status='COMPLETE' if complete else 'PARTIAL',runs=audit['runs'],
        completed_geometries=len(rows),mean_gain_vs_greedy_percent=mean,
        wins_vs_greedy=wins,wins_vs_random=random_wins,gate_passed=gate,rows=rows,
        scope='conditional exploratory TRAIN; state generation excluded; no end-to-end claim')
    (HERE/'joint_insertion_train_summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='rows'}))
