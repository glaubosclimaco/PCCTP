"""Exact fixed-path contraction of the same edge-fixing repair neighborhood.

Retained paths are mandatory supernodes. Isolated optional vertices remain
selectable. Endpoint degree equations preserve orientation without orientation
variables; SCF on the contracted multigraph excludes disconnected tours.
This is a model reduction, not a claim of a new neighborhood or theorem.
"""
from itertools import combinations
import time
import numpy as np
from scipy.optimize import Bounds,LinearConstraint,milp
from pilot import free_edge_set,route_cost,sparse_rows,validate

class FragmentRepair:
    def __init__(self,d):
        self.d=d; self.n=len(d['visits'])
        self.edges=list(combinations(range(self.n),2))

    def solve(self,route,released,seconds,operator='incident',extra_edges=None,seed=42):
        began=time.perf_counter(); d=self.d; selected=set(route)
        current={tuple(sorted((route[i-1],v))) for i,v in enumerate(route)}
        free=free_edge_set(d,route,released,operator,extra_edges,seed)
        fixed=current-free
        adjacency={v:[] for v in range(self.n)}
        for a,b in fixed: adjacency[a].append(b); adjacency[b].append(a)
        components=[]; owner={}
        for v in range(self.n):
            if v in owner: continue
            stack=[v]; component=set()
            while stack:
                u=stack.pop()
                if u in component: continue
                component.add(u); stack.extend(adjacency[u])
            index=len(components)
            for u in component: owner[u]=index
            components.append(component)
        # A retained cycle cannot connect any additional vertex and is unique.
        closed=any(len(c)>=3 and all(len(adjacency[v])==2 for v in c) for c in components)
        nodes=[c for c in components if c & (selected|released)]
        owner={v:i for i,c in enumerate(nodes) for v in c}
        q=len(nodes)
        usable=[e for e in sorted(free) if all(v in owner and len(adjacency[v])<2 for v in e)]
        m=len(usable); f0=m+q
        arcs=[(e,owner[a],owner[b]) for e,(a,b) in enumerate(usable) if owner[a]!=owner[b]]
        width=f0+2*len(arcs)
        cost0=sum(d['costs'][a,b] for a,b in fixed)
        lower=np.zeros(width); upper=np.r_[np.ones(f0),np.full(2*len(arcs),max(0,q-1))]
        for i,c in enumerate(nodes):
            v=min(c)
            if len(c)>1 or c & d['mandatory'] or (v in selected and v not in released):
                lower[m+i]=1.
        rows=[]; lo=[]; hi=[]
        def add(row,a=-np.inf,b=np.inf): rows.append(row); lo.append(a); hi.append(b)
        for v in owner:
            if len(adjacency[v])==2: continue
            row={e:1. for e,edge in enumerate(usable) if v in edge}
            if len(adjacency[v])==1: add(row,1.,1.)
            else:
                row[m+owner[v]]=-2.; add(row,0.,0.)
        for e,(a,b) in enumerate(usable):
            for i in set((owner[a],owner[b])): add({e:1.,m+i:-1.},b=0.)
        root=owner[min(d['mandatory'])]
        for i in range(q):
            row={}
            for j,(e,a,b) in enumerate(arcs):
                if i==a: row[f0+2*j]=-1.; row[f0+2*j+1]=1.
                if i==b: row[f0+2*j]=1.; row[f0+2*j+1]=-1.
            if i==root:
                for j in range(q):
                    if j!=root: row[m+j]=1.
            else: row[m+i]=-1.
            add(row,0.,0.)
        for j,(e,a,b) in enumerate(arcs):
            for index in (f0+2*j,f0+2*j+1): add({index:1.,e:-(q-1.)},b=0.)
        add({m+i:float(sum(d['prizes'][list(c)])) for i,c in enumerate(nodes)},a=d['quota'])
        add({m+i:float(len(c)) for i,c in enumerate(nodes)},a=3.)
        for ns in d['neighborhoods']:
            add({m+i:1. for i,c in enumerate(nodes) if c&ns},a=1.)
        obj=np.r_[[d['costs'][a,b] for a,b in usable],np.zeros(width-m)]
        solve_began=time.perf_counter()
        if closed:
            status=0; returned=True; best=route; message='Retained incumbent cycle'
        else:
            result=milp(obj,integrality=np.r_[np.ones(f0),np.zeros(width-f0)],
                        bounds=Bounds(lower,upper),constraints=LinearConstraint(sparse_rows(rows,width),lo,hi),
                        options=dict(time_limit=seconds,mip_rel_gap=0.,disp=False))
            status=int(result.status); returned=result.x is not None; message=result.message
            if status not in (0,1): raise RuntimeError(f'Feasible fragment repair failed: {message}')
            best=route
            if returned:
                if np.max(np.abs(result.x[:f0]-np.rint(result.x[:f0])))>1e-5:
                    raise ValueError('Nonintegral fragment solution')
                chosen=fixed|{edge for e,edge in enumerate(usable) if result.x[e]>.5}
                adj={}
                for a,b in chosen: adj.setdefault(a,[]).append(b); adj.setdefault(b,[]).append(a)
                if any(len(neighbors)!=2 for neighbors in adj.values()): raise ValueError('Degree failure')
                decoded=[min(d['mandatory'])]; previous=None
                while True:
                    nxt=next(v for v in adj[decoded[-1]] if v!=previous)
                    if nxt==decoded[0]: break
                    if nxt in decoded: raise ValueError('Repeated cycle vertex')
                    previous=decoded[-1]; decoded.append(nxt)
                if set(decoded)!=set(adj): raise ValueError('Disconnected fragment tour')
                value=validate(d,decoded)
                if abs(value-(result.fun+cost0))>1e-5: raise ValueError('Objective mismatch')
                for v in range(self.n):
                    if v not in released and ((v in decoded)!=(v in selected)):
                        raise ValueError('Fixed visit changed')
                if any(e not in free and ((e in chosen)!=(e in current)) for e in self.edges):
                    raise ValueError('Fixed edge changed')
                if value<route_cost(d,route): best=decoded
        return dict(operator=operator,status=status,returned_solution=returned,proven_optimal=status==0,
                    solve_seconds=time.perf_counter()-solve_began,total_seconds=time.perf_counter()-began,
                    free_x=len(free),free_y=len(released),cost=validate(d,best),
                    removed_optional=[int(d['visits'][v]) for v in sorted(selected-set(best))],
                    added_optional=[int(d['visits'][v]) for v in sorted(set(best)-selected)],
                    route=[int(d['visits'][v]) for v in best],
                    contracted_nodes=q,edge_variables=m,total_variables=width,
                    fixed_edge_constant=float(cost0))
