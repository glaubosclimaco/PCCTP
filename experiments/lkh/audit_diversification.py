"""Full-batch deadline and comparison audit for exploratory diversification."""
import json,statistics,hashlib
from pathlib import Path
HERE=Path(__file__).resolve().parent
if __name__=='__main__':
    folder=HERE/'diversification';assert (folder/'completion.json').exists(),'Batch still running'
    scalar=json.loads((HERE/'diversification_raw_audit.json').read_text())
    assert scalar['status']=='PASS' and scalar['batch_complete'] and scalar['files']==27
    previous=json.loads((HERE/'support_common_fallback_audit.json').read_text())
    rows=[];hashes={}
    for case in previous['rows']:
        name=case['instance'];arms={}
        for mode in ('fixed_lkh','support_2opt','support_lkh'):
            p=folder/f'{name}.{mode}.json';r=json.loads(p.read_text());hashes[p.name]=hashlib.sha256(p.read_bytes()).hexdigest()
            assert scalar['result_hashes'][p.name]==hashes[p.name]
            assert r['trace'] and r['cost']==r['trace'][-1]['cost'] and r['route']==r['trace'][-1]['route']
            last=-1;cost=float('inf');support=set(r['trace'][0]['route'])
            for e in r['trace']:
                assert last<=e['seconds']<=15 and e['cost']<cost
                if mode=='fixed_lkh':assert set(e['route'])==support
                last=e['seconds'];cost=e['cost']
            arms[mode]=dict(cost=r['cost'],iterations=r['iterations'],wall_seconds=r['wall_seconds'])
        baseline=statistics.mean(p['arms']['mdm']['cost'] for p in case['pairs'])
        rows.append(dict(instance=name,arms=arms,cached_mdm_mean=baseline,gain_percent=100*(baseline-arms['support_lkh']['cost'])/baseline))
    mean=statistics.mean(r['gain_percent'] for r in rows);wins=sum(r['gain_percent']>0 for r in rows)
    report=dict(status='PASS',scope='TRAIN viability with historical cached baseline; not contemporaneous paired confirmation',rows=rows,mean_gain_percent=mean,wins_vs_cached_mdm=wins,gate_passed=mean>=.5 and wins>=6,result_hashes=hashes)
    (HERE/'diversification_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('rows','result_hashes')}))
