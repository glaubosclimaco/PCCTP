"""Fixed random Fourier features with a learned pairwise ridge output head.

Development alternative to linear ranking, with no extra library dependency.
Scaling and output weights use training geometries only. Architecture, seed,
bandwidth and regularization are fixed, not optimized against held-out labels.
Geometry weighting prevents many snapshots from multiplying a geometry's weight.
"""
import numpy as np

def transform(x,model):
    z=np.clip((np.asarray(x)[:,model['indices']]-model['mean'])/model['scale'],-10.,10.)
    if model['kind']=='linear': return z
    projected=np.sqrt(2./model['width'])*np.cos(z@model['frequencies']+model['phases'])
    return np.column_stack((z,projected))

def fit(states,indices,kind='rff',width=64,alpha=1.,seed=2026):
    if kind not in ('linear','rff'): raise ValueError('Unknown ranker')
    if not states: raise ValueError('No training states')
    x=np.array([r['features'] for s in states for r in s['records']])[:,indices]
    scale=x.std(axis=0); scale[scale<1e-10]=1.
    rng=np.random.default_rng(seed)
    model=dict(kind=kind,indices=indices,mean=x.mean(axis=0),scale=scale,width=width,
               frequencies=rng.normal(size=(len(indices),width))/np.sqrt(len(indices)),
               phases=rng.uniform(0,2*np.pi,size=width),alpha=alpha,seed=seed)
    grouped={}
    for state in states:
        records=state['records']
        z=transform([r['features'] for r in records],model)
        xx=[]; yy=[]
        for i,r in enumerate(records):
            for j,other in enumerate(records[:i]):
                delta=r['normalized_gain']-other['normalized_gain']
                if abs(delta)<1e-12: continue
                xx.append(z[i]-z[j]); yy.append(delta)
        if xx: grouped.setdefault(state['geometry'],[]).append((np.array(xx),np.array(yy)))
    dimension=len(indices)+(width if kind=='rff' else 0)
    gram=alpha*np.eye(dimension); rhs=np.zeros(dimension)
    for geometry,blocks in grouped.items():
        pairs=sum(len(y) for _,y in blocks)
        for z,y in blocks:
            gram+=z.T@z/pairs
            rhs+=z.T@y/pairs
    model['weights']=np.linalg.solve(gram,rhs)
    model['training_geometries']=sorted({s['geometry'] for s in states})
    model['informative_training_geometries']=sorted(grouped)
    return model

def predict(model,records):
    return transform([r['features'] for r in records],model)@model['weights']
