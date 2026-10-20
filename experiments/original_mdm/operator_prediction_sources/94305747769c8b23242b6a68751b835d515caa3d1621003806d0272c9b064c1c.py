"""Freeze TRAIN-only models and verify native inference against Python.

No scheduling policy is implemented here. The measured native call excludes
route feature extraction and cannot be reported as the full policy overhead.
"""
import ctypes
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import time
import numpy as np
from scipy.special import expit
from operator_prediction import state_features,fit,design,NAMES
from experiments.ml_lns.pilot import read_instance,validate

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profile-dir',type=Path,default=HERE/'move_profiles')
    parser.add_argument('--require-monotonic-elite',action='store_true')
    args=parser.parse_args();profile_dir=args.profile_dir.resolve()
    if args.require_monotonic_elite:
        assert list(profile_dir.glob('completion.*.json')),'Complete the profile collection before exporting models'
    versions=HERE/'operator_prediction_sources';versions.mkdir(exist_ok=True)
    for name in ['export_operator_models.py','operator_prediction.py','operator_inference.cpp']:
        path=HERE/name;content=path.read_bytes()
        (versions/(hashlib.sha256(content).hexdigest()+path.suffix)).write_bytes(content)
    blocks=[];sources={};seen=set()
    for path in sorted(profile_dir.glob('g_*.json')):
        s=json.loads(path.read_text())
        if 'records' not in s or s['role']!='train':continue
        if args.require_monotonic_elite:
            assert s.get('clock_source')=='clock_gettime(CLOCK_MONOTONIC)' and s.get('elite_count_bounded') is True
        key=(s['instance'],s['seed']);assert key not in seen;seen.add(key)
        source=ROOT/'experiments/ml_lns/decoder_initialized_states'/f"{s['instance']}.seed{s['seed']}.json"
        assert hashlib.sha256(source.read_bytes()).hexdigest()==s['source_sha256']
        assert json.loads(source.read_text())['role']=='train'
        d=read_instance(ROOT/'data'/s['instance']);mapping={int(v):i for i,v in enumerate(d['visits'])}
        x=[];moves=[];y=[];t=[];cache={}
        for r in s['records']:
            route=tuple(mapping[v] for v in r['before_route'])
            assert validate(d,list(route))==r['before_cost']
            assert validate(d,[mapping[v] for v in r['after_route']])==r['after_cost']
            if route not in cache:cache[route]=state_features(d,list(route))
            x.append(cache[route]);moves.append(r['move']);y.append(float(r['after_cost']<r['before_cost']))
            t.append(np.log(max(1e-6,r['operator_seconds'])))
        blocks.append(dict(geometry=s['geometry'],x=np.array(x),moves=np.array(moves),y=np.array(y),t=np.array(t)))
        sources[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
    geometries=sorted({b['geometry'] for b in blocks});assert len(geometries)==20
    x=np.concatenate([b['x'] for b in blocks]);moves=np.concatenate([b['moves'] for b in blocks])
    y=np.concatenate([b['y'] for b in blocks]);t=np.concatenate([b['t'] for b in blocks])
    counts={g:sum(len(b['y']) for b in blocks if b['geometry']==g) for g in geometries}
    weight=np.concatenate([np.full(len(b['y']),1./counts[b['geometry']]/len(geometries)) for b in blocks])
    models={name:fit(x,moves,y,t,weight,context) for name,context in
            [('operator_only',False),('operator_size','size'),
             ('operator_generic','generic'),('operator_state',True)]}
    build=HERE/'operator_inference_build';build.mkdir(exist_ok=True)
    source=HERE/'operator_inference.cpp';library=build/'operator_inference.dylib'
    subprocess.run([shutil.which('g++-15') or 'g++','-std=c++11','-O3','-dynamiclib',str(source),'-o',str(library)],check=True)
    native=ctypes.CDLL(str(library)).pcctp_operator_scores
    pointer=ctypes.POINTER(ctypes.c_double)
    native.argtypes=[pointer,pointer,pointer,pointer,pointer,ctypes.c_int,pointer,pointer];native.restype=None
    def ptr(a):return a.ctypes.data_as(pointer)
    rng=np.random.default_rng(2026);test=x[rng.choice(len(x),100,replace=False)]
    test=np.vstack((test,x.mean(axis=0)+20*x.std(axis=0),x.mean(axis=0)-20*x.std(axis=0)))
    checks=0;maxerror=0.;timings={}
    for mode,(name,model) in enumerate(models.items()):
        arrays=[np.ascontiguousarray(model[k],dtype=float) for k in ('mean','scale','probability','timing')]
        p=np.empty(15);lt=np.empty(15)
        for features in test:
            features=np.ascontiguousarray(features,dtype=float)
            native(ptr(features),*[ptr(a) for a in arrays],mode,ptr(p),ptr(lt))
            a=design(np.tile(features,(15,1)),np.arange(15),model['mean'],model['scale'],model['context'])
            error=max(float(np.max(abs(p-expit(a@model['probability'])))),float(np.max(abs(lt-a@model['timing']))))
            assert error<1e-10;maxerror=max(maxerror,error);checks+=15
        began=time.perf_counter()
        for repeat in range(1000):native(ptr(test[0]),*[ptr(a) for a in arrays],mode,ptr(p),ptr(lt))
        timings[name]=(time.perf_counter()-began)/1000
    hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in
            ['export_operator_models.py','operator_prediction.py','operator_inference.cpp']}
    payload=dict(scope=__doc__,training_geometries=geometries,observations=len(x),source_hashes=sources,
                 profile_directory=str(profile_dir),requires_monotonic_elite=args.require_monotonic_elite,
                 code_hashes=hashes,feature_names=NAMES,models={name:{k:v.tolist() if isinstance(v,np.ndarray) else v
                              for k,v in model.items()} for name,model in models.items()},
                 native_verification=dict(status='PASS',action_scores_checked=checks,max_absolute_error=maxerror,
                                          python_ctypes_call_seconds=timings,includes_route_feature_extraction=False))
    fingerprint=hashlib.sha256(json.dumps(payload,sort_keys=True).encode()).hexdigest()[:12]
    out=HERE/'operator_models';out.mkdir(exist_ok=True)
    (out/f'models.{fingerprint}.json').write_text(json.dumps(payload,indent=2)+'\n')
    print(json.dumps(dict(file=str(out/f'models.{fingerprint}.json'),training_geometries=len(geometries),
                         observations=len(x),verification=payload['native_verification']),indent=2))

if __name__=='__main__':main()
