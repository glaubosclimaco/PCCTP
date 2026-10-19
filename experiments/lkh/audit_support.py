"""Saved-route audit for fixed-state support diagnostic."""
import json,sys,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from experiments.ml_lns.pilot import read_instance,validate
HERE=Path(__file__).resolve().parent
if __name__=='__main__':
    folder=HERE/'support_diagnostic';assert (folder/'completion.json').exists()
    source=hashlib.sha256((HERE/'support_search.py').read_bytes()).hexdigest()
    cases=json.loads((ROOT/'experiments/original_mdm/hard_training_baseline/summary.json').read_text())['rows'];rows=[];hashes={}
    for c in cases:
        name=c['instance'];d=read_instance(ROOT/'data'/name);mapping={int(v):i for i,v in enumerate(d['visits'])};arms={}
        initial=json.loads((HERE/'total_time'/f'{name}.seed42.mdm.json').read_text())
        for mode in ('fixed_lkh','support_2opt','support_lkh'):
            p=folder/f'{name}.{mode}.json';r=json.loads(p.read_text());hashes[p.name]=hashlib.sha256(p.read_bytes()).hexdigest()
            assert r['source_hash']==source and r['initial_cost']==initial['cost']
            assert validate(d,[mapping[v] for v in r['route']])==r['cost']<=r['initial_cost']
            previous=r['initial_cost'];last=-1
            for i,e in enumerate(r['trace']):
                assert last<=e['seconds']<=5 and e['cost']<=previous
                if i:assert e['cost']<previous
                assert validate(d,[mapping[v] for v in e['route']])==e['cost']
                previous=e['cost'];last=e['seconds']
            assert r['trace'][-1]['cost']==r['cost'] and r['trace'][-1]['route']==r['route']
            if mode=='fixed_lkh':assert set(r['route'])==set(initial['route'])
            arms[mode]=r['cost']
        rows.append(dict(instance=name,costs=arms,published_mdm_best=c['published']['best'],support_beats_fixed=arms['support_lkh']<arms['fixed_lkh'],support_beats_2opt=arms['support_lkh']<arms['support_2opt']))
    result=dict(status='PASS',scope='fixed-state additional-time TRAIN diagnostic; not end-to-end or global record verification',rows=rows,result_hashes=hashes,wins_vs_fixed=sum(r['support_beats_fixed'] for r in rows),wins_vs_2opt=sum(r['support_beats_2opt'] for r in rows),below_supplied_mdm_table=sum(r['costs']['support_lkh']<r['published_mdm_best'] for r in rows))
    result['gate_passed']=result['wins_vs_fixed']>=3
    (HERE/'support_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('rows','result_hashes')}))
