"""Development-only grouped ranking audit. Never opens validation/test labels.

States from a geometry stay together. Summary weights geometries equally;
candidate rows and snapshots are not independent experimental units.
"""
import json
from pathlib import Path
import numpy as np
from evaluate_ranker import fit_predict

HERE=Path(__file__).resolve().parent

def main():
    states=[]
    for path in sorted((HERE/'strong_oracle').glob('*.json')):
        row=json.loads(path.read_text())
        if row['role'] not in ('development','train'):
            raise ValueError('Reserved labels must not enter this exploratory audit')
        if not row['records']: continue
        row['source']=path.name
        states.append(row)
    groups=sorted({s['geometry'] for s in states})
    if len(groups)<3: raise ValueError('Insufficient geometries')
    generic=states[0]['generic_feature_indices']
    feature_sets=dict(generic=generic,raw_coverage_prize=list(range(12)),
                      joint_only=generic+list(range(12,22)),all=list(range(22)))
    rows=[]
    for geometry in groups:
        train=[s for s in states if s['geometry']!=geometry]
        for test in [s for s in states if s['geometry']==geometry]:
            gains=np.array([r['normalized_gain'] for r in test['records']])
            selected=dict(uniform=float(gains.mean()),sampled_oracle=float(gains.max()))
            for label,index in [('deletion_savings',6),('joint_savings',20),('joint_flexibility',12)]:
                selected[label]=float(gains[np.argmax([r['features'][index] for r in test['records']])])
            for method,indices in feature_sets.items():
                scores=fit_predict(train,test,indices)
                selected[method]=float(gains[int(np.argmax(scores))])
            rows.append(dict(geometry=geometry,phase=test['phase'],source=test['source'],gain=selected))
    summary={}
    for method in rows[0]['gain']:
        values=[np.mean([r['gain'][method] for r in rows if r['geometry']==g]) for g in groups]
        summary[method]=dict(geometry_mean_gain=float(np.mean(values)),
                             positive_geometries=int(sum(v>1e-10 for v in values)))
    result=dict(scope='Development-only leave-one-geometry-out; not confirmatory or end-to-end',
                geometries=len(groups),states=len(states),
                queries=sum(len(s['records']) for s in states),
                informative_states=sum(any(r['normalized_gain']>1e-10 for r in s['records']) for s in states),
                ridge_alpha=1.,summary=summary,rows=rows)
    (HERE/'strong_development_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))

if __name__=='__main__': main()
