"""Exploratory nonlinear ranking, with held-out geometries and strong simple controls."""
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from nonlinear_ranker import fit,predict

HERE=Path(__file__).resolve().parent

def main():
    began=time.perf_counter(); states=[]; manifest={}
    for path in sorted((HERE/'strong_oracle').glob('*.json')):
        row=json.loads(path.read_text())
        if row['role'] not in ('development','train'): raise ValueError('Reserved labels')
        states.append(row); manifest[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
    geometries=sorted({s['geometry'] for s in states})
    if len(geometries)<3: raise ValueError('Insufficient geometries')
    generic=states[0]['generic_feature_indices']
    indices=dict(generic=generic,raw=list(range(12)),joint=generic+list(range(12,22)))
    rows=[]
    for geometry in geometries:
        train=[s for s in states if s['geometry']!=geometry]
        models={f'{kind}_{name}':fit(train,columns,kind) for kind in ('linear','rff') for name,columns in indices.items()}
        for test in [s for s in states if s['geometry']==geometry]:
            gains=np.array([r['normalized_gain'] for r in test['records']])
            chosen=dict(uniform=float(gains.mean()),sampled_oracle=float(gains.max()))
            for name,index in [('deletion_savings',6),('joint_savings',20)]:
                chosen[name]=float(gains[np.argmax([r['features'][index] for r in test['records']])])
            for name,model in models.items(): chosen[name]=float(gains[np.argmax(predict(model,test['records']))])
            rows.append(dict(geometry=geometry,instance=test['instance'],phase=test['phase'],gain=chosen))
    summary={}
    for method in rows[0]['gain']:
        values=[np.mean([r['gain'][method] for r in rows if r['geometry']==g]) for g in geometries]
        summary[method]=dict(mean_geometry_gain=float(np.mean(values)),positive_geometries=int(sum(v>1e-10 for v in values)))
    result=dict(scope='Exploratory development-only LOGO; no confirmatory or end-to-end claim',
                architecture=dict(width=64,alpha=1.,seed=2026,bandwidth='sqrt(input dimension)',
                                  scaling='train only, standardized then clipped to [-10,10]',
                                  pair_weighting='equal per informative geometry'),
                geometries=len(geometries),states=len(states),source_hashes=manifest,
                elapsed_seconds=time.perf_counter()-began,summary=summary,rows=rows)
    (HERE/'nonlinear_development_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('rows','source_hashes')},indent=2))

if __name__=='__main__': main()
