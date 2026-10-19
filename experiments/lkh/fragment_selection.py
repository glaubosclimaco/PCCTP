"""Fixed-ridge geometry-held-out exploratory selection, before integration."""
import sys,json,statistics
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'experiments/ml_lns'))
import numpy as np
from pilot import read_instance,features
HERE=Path(__file__).resolve().parent
if __name__=='__main__':
    audit=json.loads((HERE/'hard_fragment_audit.json').read_text());assert audit['status']=='PASS' and audit['gate_passed']
    groups={}
    for p in sorted((HERE/'hard_fragment').glob('*.neighbor*.json')):
        r=json.loads(p.read_text());name=r['instance'];d=read_instance(ROOT/'data'/name);mapping={int(v):i for i,v in enumerate(d['visits'])};route=[mapping[v] for v in r['initial_route']];released={mapping[v] for v in r['released']}
        x=features(d,route,released)+[int(r['kind']==k) for k in ('random','coverage','prize')]
        y=100*(r['initial_cost']-r['result']['cost'])/r['initial_cost']
        groups.setdefault(name,[]).append(dict(index=r['index'],kind=r['kind'],x=x,y=y))
    rows=[]
    for name,test in groups.items():
        train=[r for g,rs in groups.items() if g!=name for r in rs];x=np.array([r['x'] for r in train]);y=np.array([r['y'] for r in train]);mean=x.mean(axis=0);scale=x.std(axis=0);scale[scale<1e-10]=1
        z=(x-mean)/scale;center=y.mean();w=np.linalg.solve(z.T@z+np.eye(z.shape[1]),z.T@(y-center));pred=((np.array([r['x'] for r in test])-mean)/scale)@w+center
        chosen=max(range(len(test)),key=lambda i:(float(pred[i]),-test[i]['index']));random_gain=statistics.mean(r['y'] for r in test if r['kind']=='random')
        rows.append(dict(instance=name,selected_index=test[chosen]['index'],selected_kind=test[chosen]['kind'],selected_gain=test[chosen]['y'],random_expected_gain=random_gain,uniform_expected_gain=statistics.mean(r['y'] for r in test),delta_vs_random=test[chosen]['y']-random_gain,predictions=pred.tolist()))
    delta=statistics.mean(r['delta_vs_random'] for r in rows);wins=sum(r['delta_vs_random']>1e-10 for r in rows)
    report=dict(scope='exploratory TRAIN leave-one-geometry-out; no hyperparameter tuning, not confirmatory',alpha=1,rows=rows,mean_delta_vs_random_pp=delta,wins_vs_random=wins,gate_passed=delta>=.1 and wins>=6)
    (HERE/'fragment_selection_results.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='rows'}))
