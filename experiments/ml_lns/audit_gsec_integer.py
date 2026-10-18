"""Saved-data audit for root strengthening and its paired integer pilot."""
import json,hashlib
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from experiments.ml_lns.pilot import read_instance,validate


if __name__=='__main__':
    here=Path(__file__).resolve().parent
    root=json.loads((here/'gsec_root_results.json').read_text())
    assert len(root['rows'])==11
    for row in root['rows']:
        b=row['bounds'];assert all(y>=x-1e-6 for x,y in zip(b,b[1:]))
        assert all(v<=row['mdm_cost']+1e-5 for v in b)
        assert len(row['cuts'])==row['cut_count']
    assert sum(r['passes'] for r in root['rows'])==root['passing']
    path=here/'gsec_integer_results.json';result=json.loads(path.read_text())
    assert len(result['rows'])==11,'Incomplete batch'
    seen=set();proofs={'plain':0,'cuts':0};gains={'plain':0,'cuts':0}
    for k,row in enumerate(result['rows']):
        assert row['geometry'] not in seen;seen.add(row['geometry'])
        assert row['order']==(['plain','cuts'] if k%2==0 else ['cuts','plain'])
        d=read_instance(ROOT/'data'/row['instance']);mapping={int(v):i for i,v in enumerate(d['visits'])}
        for mode in ('plain','cuts'):
            r=row[mode];cost=validate(d,[mapping[v] for v in r['route']])
            assert abs(cost-r['cost'])<1e-6 and cost<=r['initial_cost']+1e-6
            assert r['warm_start'] and r['start_status']=='HighsStatus.kOk' and r['version']=='1.13.1'
            assert r['dual_bound']<=cost+1e-4
            assert r['separation_seconds']<=r['seconds']
            if mode=='plain':assert r['root_cut_count']==0 and r['separation_seconds']==0
            last_cost=r['initial_cost'];last_time=0.
            for event in r['events']:
                assert last_time<=event['seconds']<=10. and event['cost']<last_cost-1e-7
                last_time=event['seconds'];last_cost=event['cost']
            assert abs(last_cost-cost)<1e-6
            assert r['improved']==(cost<r['initial_cost']-1e-7)
            if r['proven_optimal']:
                assert r['status']=='Optimal' and r['seconds']<=10. and abs(cost-r['dual_bound'])<1e-4
            proofs[mode]+=int(r['proven_optimal']);gains[mode]+=int(r['improved'])
    extra=proofs['cuts']-proofs['plain']
    no_loss=all(r['cuts']['cost']<=r['plain']['cost']+1e-6 for r in result['rows'])
    gate=(gains['cuts']>=2 and gains['cuts']>gains['plain']) or (extra>=2 and no_loss)
    assert result['cut_gains']==gains['cuts'] and result['plain_gains']==gains['plain']
    assert result['net_extra_proofs']==extra and result['gate_passed']==gate
    audit=dict(status='PASS',geometries=11,root_passing=root['passing'],proofs=proofs,gains=gains,
               net_extra_proofs=extra,gate_passed=gate,result_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    (here/'gsec_integer_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps(audit))
