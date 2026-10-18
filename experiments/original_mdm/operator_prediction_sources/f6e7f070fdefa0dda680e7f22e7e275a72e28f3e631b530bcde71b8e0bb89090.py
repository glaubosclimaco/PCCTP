"""Training-geometry LOGO prediction; observational accuracy is not policy value."""
import hashlib
import argparse
import json
from pathlib import Path
import sys
import time
import numpy as np
from scipy.optimize import minimize
from scipy.special import expit

ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from experiments.ml_lns.pilot import read_instance,validate

NAMES=['log_visits','log_route_size','selected_fraction','mandatory_fraction',
       'quota_fraction','prize_ratio','tight_coverage_fraction','mean_edge_ratio',
       'edge_variation','optional_deletion_saving']

def state_features(d,route):
    selected=set(route);n=len(d['visits']);mean=max(1.,float(np.mean(d['costs'])))
    edges=np.array([d['costs'][route[i-1],v] for i,v in enumerate(route)])
    saving=sum(max(0.,d['costs'][route[i-1],v]+d['costs'][v,route[(i+1)%len(route)]]-
                   d['costs'][route[i-1],route[(i+1)%len(route)]])
               for i,v in enumerate(route) if v not in d['mandatory'])
    return [np.log1p(n),np.log1p(len(route)),len(route)/n,len(d['mandatory'])/n,
            d['quota']/max(1.,float(sum(d['prizes']))),
            float(sum(d['prizes'][route]))/max(1.,d['quota']),
            sum(len(ns&selected)==1 for ns in d['neighborhoods'])/max(1,len(d['neighborhoods'])),
            float(edges.mean())/mean,float(edges.std())/max(1.,float(edges.mean())),
            saving/max(1.,float(edges.sum()))]

def design(x,moves,mean,scale,context):
    hot=np.eye(15)[moves]
    if not context:return hot
    z=np.clip((x-mean)/scale,-5.,5.)
    if context=='size':z=z[:,:2]
    return np.column_stack((hot,(hot[:,:,None]*z[:,None,:]).reshape(len(x),-1)))

def fit(x,moves,y,logtime,weight,context):
    mean=np.average(x,axis=0,weights=weight)
    scale=np.sqrt(np.average((x-mean)**2,axis=0,weights=weight));scale[scale<1e-9]=1.
    a=design(x,moves,mean,scale,context);penalty=np.full(a.shape[1],.01);penalty[:15]=.0001
    def objective(b):
        linear=a@b
        loss=np.dot(weight,np.logaddexp(0.,linear)-y*linear)+.5*np.dot(penalty*b,b)
        gradient=a.T@(weight*(expit(linear)-y))+penalty*b
        return loss,gradient
    result=minimize(objective,np.zeros(a.shape[1]),jac=True,method='L-BFGS-B',options=dict(maxiter=300,ftol=1e-10))
    if not result.success:raise RuntimeError(result.message)
    timing=np.linalg.solve(a.T@(weight[:,None]*a)+np.diag(penalty),a.T@(weight*logtime))
    return dict(mean=mean,scale=scale,probability=result.x,timing=timing,context=context)

def metrics(y,p,t,pred_t):
    p=np.clip(p,1e-8,1.-1e-8)
    return dict(brier=float(np.mean((p-y)**2)),log_loss=float(np.mean(-y*np.log(p)-(1-y)*np.log1p(-p))),
                log_time_MAE=float(np.mean(abs(t-pred_t))))

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--all-profiles',action='store_true')
    parser.add_argument('--partial',action='store_true',help='Explicit snapshot of an unfinished expanded collection')
    args=parser.parse_args()
    collection_complete=(HERE/'move_profiles/completion.093590aaa738.json').exists()
    if args.all_profiles and not args.partial and not collection_complete:
        parser.error('Expanded collection unfinished; use --partial for a labeled exploratory snapshot')
    began=time.perf_counter();blocks=[];manifest={};audit=[]
    pattern='g_*.json' if args.all_profiles else 'g_*.4e1c42c7aa48.json'
    seen=set()
    for path in sorted((HERE/'move_profiles').glob(pattern)):
        s=json.loads(path.read_text())
        if 'records' not in s:continue
        key=(s['instance'],s['seed'])
        if key in seen:raise ValueError('Duplicate profile for instance/seed')
        seen.add(key);d=read_instance(ROOT/'data'/s['instance'])
        source=ROOT/'experiments/ml_lns/decoder_initialized_states'/f"{s['instance']}.seed{s['seed']}.json"
        assert hashlib.sha256(source.read_bytes()).hexdigest()==s['source_sha256']
        assert json.loads(source.read_text())['role']==s['role']
        mapping={int(v):i for i,v in enumerate(d['visits'])};records=s['records'];x=[];moves=[];y=[];t=[]
        cache={};feature_times=[]
        for r in records:
            route=[mapping[v] for v in r['before_route']]
            assert validate(d,route)==r['before_cost']
            assert validate(d,[mapping[v] for v in r['after_route']])==r['after_cost']
            route_key=tuple(route)
            if route_key not in cache:
                feature_began=time.perf_counter();cache[route_key]=state_features(d,route)
                feature_times.append(time.perf_counter()-feature_began)
            x.append(cache[route_key]);moves.append(r['move'])
            y.append(float(r['after_cost']<r['before_cost']));t.append(np.log(max(1e-6,r['operator_seconds'])))
        audit.append(dict(instance=s['instance'],role=s['role'],calls=len(records),sample_capped=s['sample_capped'],
                          profiled_run_completed=s['profiled_run_completed'],unique_input_states=len(cache),
                          median_feature_seconds=float(np.median(feature_times))))
        manifest[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        if s['role']=='train':blocks.append(dict(geometry=s['geometry'],x=np.array(x),moves=np.array(moves),y=np.array(y),t=np.array(t)))
    assert len(blocks)>=4
    rows=[]
    geometries=sorted({b['geometry'] for b in blocks})
    for geometry in geometries:
        test_blocks=[b for b in blocks if b['geometry']==geometry]
        test={name:np.concatenate([b[name] for b in test_blocks]) for name in ('x','moves','y','t')}
        train=[b for b in blocks if b['geometry']!=geometry]
        x=np.concatenate([b['x'] for b in train]);moves=np.concatenate([b['moves'] for b in train])
        y=np.concatenate([b['y'] for b in train]);t=np.concatenate([b['t'] for b in train])
        geometry_counts={g:sum(len(b['y']) for b in train if b['geometry']==g) for g in geometries if g!=geometry}
        weight=np.concatenate([np.full(len(b['y']),1./geometry_counts[b['geometry']]/len(geometry_counts)) for b in train])
        scores={}
        for name,context in [('operator_only',False),('operator_size','size'),('operator_state',True)]:
            model=fit(x,moves,y,t,weight,context)
            a=design(test['x'],test['moves'],model['mean'],model['scale'],context)
            scores[name]=metrics(test['y'],expit(a@model['probability']),test['t'],a@model['timing'])
        rows.append(dict(geometry=geometry,observations=len(test['y']),positive=int(test['y'].sum()),metrics=scores))
    summary={name:{metric:float(np.mean([r['metrics'][name][metric] for r in rows]))
                   for metric in rows[0]['metrics'][name]} for name in rows[0]['metrics']}
    result=dict(scope='Exploratory LOGO on TRAIN geometries; observed executed actions only, no policy-value claim',
                training_geometries=len(geometries),
                expanded_collection_complete=collection_complete,partial_snapshot=bool(args.all_profiles and args.partial),
                code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),source_hashes=manifest,
                audit=audit,feature_names=NAMES,rows=rows,summary=summary,elapsed_seconds=time.perf_counter()-began,
                configuration=dict(regularization_intercepts=.0001,regularization_slopes=.01,
                                   geometry_weight='equal',state_scaling='train only; clip5',time_floor_seconds=1e-6))
    fingerprint=hashlib.sha256(json.dumps(dict(sources=manifest,code=result['code_sha256']),sort_keys=True).encode()).hexdigest()[:12]
    target=HERE/(f'operator_prediction.{fingerprint}'+('.partial_snapshot' if args.partial else '')+'.json' if args.all_profiles else 'operator_prediction_development.json')
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('source_hashes','rows','audit')},indent=2))

if __name__=='__main__':main()
