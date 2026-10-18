"""Exploratory leave-one-training-geometry-out ranking for the new DP task."""
import hashlib
import json
from pathlib import Path
import numpy as np
from nonlinear_ranker import fit,predict

HERE=Path(__file__).resolve().parent

def main():
    states=[];hashes={}
    for path in sorted((HERE/'dp_training_labels').glob('g_*.json')):
        state=json.loads(path.read_text());assert state['role']=='train'
        states.append(state);hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
    assert len(states)>=10 and len({s['geometry'] for s in states})==len(states)
    rows=[]
    for test in states:
        train=[s for s in states if s['geometry']!=test['geometry']]
        gains=np.array([r['normalized_gain'] for r in test['records']]);x=np.array([r['features'] for r in test['records']])
        eligible=x[:,8]<=1.
        choices=dict(uniform=float(gains.mean()),sampled_oracle=float(gains.max()),
                     minimum_components=float(gains[np.argmin(x[:,8])]),
                     capped_uniform=float(gains[eligible].mean()) if eligible.any() else 0.,
                     capped_joint_saving=float(gains[np.argmax(np.where(eligible,x[:,20],-np.inf))]) if eligible.any() else 0.)
        for kind in ('linear','rff'):
            for name,columns in [('structural',list(range(12))),('joint',list(range(22)))]:
                model=fit(train,columns,kind)
                assert test['geometry'] not in model['training_geometries']
                scores=predict(model,test['records'])
                choices[kind+'_'+name]=float(gains[np.argmax(np.where(eligible,scores,-np.inf))]) if eligible.any() else 0.
        rows.append(dict(instance=test['instance'],geometry=test['geometry'],gain=choices))
    summary={m:dict(mean_geometry_gain=float(np.mean([r['gain'][m] for r in rows])),
                    positive_geometries=sum(r['gain'][m]>0 for r in rows)) for m in rows[0]['gain']}
    result=dict(scope='Exploratory training-only LOGO; not reserved-test or end-to-end evidence',
                source_hashes=hashes,states=len(states),rows=rows,summary=summary,
                controls='All learned choices and capped rules know the same 18-component limit; uniform includes limits',
                architecture=dict(width=64,alpha=1.,seed=2026,selection='first maximum, no tuning'))
    (HERE/'dp_ranking_training_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('rows','source_hashes')},indent=2))

if __name__=='__main__':main()
