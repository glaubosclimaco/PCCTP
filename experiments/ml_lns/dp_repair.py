"""Exact fixed-path repair with subset dynamic programming and feasible selections.

q<=18 by default for bounded memory; size/time limits preserve the incumbent.
Both path orientations are considered; undirected symmetry fixes root orientation.
All costs must be nonnegative. No assertion that classic Held--Karp is novel.
"""
import ctypes
from itertools import combinations
from pathlib import Path
import shutil
import subprocess
import time
import numpy as np
from pilot import free_edge_set,route_cost,validate

HERE=Path(__file__).resolve().parent
LIB=None

def library():
    global LIB
    if LIB is None:
        folder=HERE/'dp_build'; folder.mkdir(exist_ok=True)
        source=HERE/'path_dp.cpp'; target=folder/'path_dp.dylib'
        if not target.exists() or target.stat().st_mtime<source.stat().st_mtime:
            subprocess.run([shutil.which('g++-15') or 'g++','-std=c++11','-O3','-dynamiclib',
                            str(source),'-o',str(target)],check=True,capture_output=True,text=True)
        LIB=ctypes.CDLL(str(target))
        LIB.pcctp_path_dp.argtypes=[ctypes.c_int,ctypes.c_int,ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_ubyte),ctypes.c_double,ctypes.c_double,
            ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_uint64)]
        LIB.pcctp_path_dp.restype=ctypes.c_int
    return LIB

class DPRepair:
    def __init__(self,d,max_components=18):
        if not 1<=max_components<=20: raise ValueError('Native DP supports 1..20 components')
        if np.any(d['costs']<0): raise ValueError('DP pruning requires nonnegative costs')
        if not np.allclose(d['costs'],d['costs'].T,rtol=0.,atol=1e-9):
            raise ValueError('Root orientation symmetry requires undirected symmetric costs')
        self.d=d; self.max_components=max_components
        library() # Build/load explicitly counted by a caller constructing the model online.

    def solve(self,route,released,seconds,operator='route_context',extra_edges=None,seed=42):
        began=time.perf_counter(); d=self.d; released=set(released)
        current={tuple(sorted((route[i-1],v))) for i,v in enumerate(route)}
        free=free_edge_set(d,route,released,operator,extra_edges,seed); fixed=current-free
        graph={v:[] for v in range(len(d['visits']))}
        for a,b in fixed:graph[a].append(b);graph[b].append(a)
        components=[]; seen=set()
        for v in sorted(set(route)|released):
            if v in seen:continue
            stack=[v]; component=set()
            while stack:
                u=stack.pop()
                if u in component:continue
                component.add(u);stack.extend(graph[u])
            seen|=component
            ends=sorted(u for u in component if len(graph[u])<2)
            if not ends:
                return self.result(route,route,released,free,0,began,len(components)+1,0,'fixed_cycle')
            path=[ends[0]]; previous=None
            while True:
                candidates=[w for w in graph[path[-1]] if w!=previous]
                if not candidates:break
                previous,pathnext=path[-1],candidates[0];path.append(pathnext)
            assert set(path)==component
            components.append(path)
        q=len(components)
        if q>self.max_components:
            return self.result(route,route,released,free,1,began,q,0,'component_limit')
        mandatory=0; optional=[]; owner={v:i for i,c in enumerate(components) for v in c}
        for i,c in enumerate(components):
            if len(c)>1 or c[0] not in released or c[0] in d['mandatory']:mandatory|=1<<i
            else:optional.append((i,c[0]))
        root=owner[min(d['mandatory'])]
        feasible=np.zeros(1<<q,dtype=np.uint8)
        outside={v for i,c in enumerate(components) if mandatory&(1<<i) for v in c}
        assignments=((np.arange(1<<len(optional))[:,None]>>np.arange(len(optional)))&1).astype(np.int8)
        valid=assignments.sum(axis=1)+len(outside)>=3
        valid &= assignments@np.array([d['prizes'][v] for i,v in optional],dtype=float)+sum(d['prizes'][list(outside)])>=d['quota']-1e-6
        residual=[ns for ns in d['neighborhoods'] if not ns&outside]
        if residual:
            incidence=np.array([[v in ns for i,v in optional] for ns in residual],dtype=np.int16)
            valid &= (assignments@incidence.T>=1).all(axis=1)
        masks=mandatory+assignments@np.array([1<<i for i,v in optional],dtype=np.int64)
        feasible[masks[valid]]=1
        orientation=np.array([2 if len(c)>1 else 1 for c in components],dtype=np.int32)
        transitions=np.full((2*q,2*q),np.inf,dtype=float)
        for i,path in enumerate(components):
            for o in range(orientation[i]):
                exit=path[-1] if o==0 else path[0]
                for j,other in enumerate(components):
                    for p in range(orientation[j]):
                        entry=other[0] if p==0 else other[-1]
                        if exit!=entry and tuple(sorted((exit,entry))) in free:
                            transitions[2*i+o,2*j+p]=d['costs'][exit,entry]
        constant=sum(d['costs'][a,b] for a,b in fixed)
        order=np.empty(q,dtype=np.int32); length=ctypes.c_int(); value=ctypes.c_double(); processed=ctypes.c_uint64()
        remaining=max(.001,seconds-(time.perf_counter()-began))
        status=library().pcctp_path_dp(q,root,transitions.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            orientation.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),feasible.ctypes.data_as(ctypes.POINTER(ctypes.c_ubyte)),
            route_cost(d,route)-constant,remaining,order.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            ctypes.byref(length),ctypes.byref(value),ctypes.byref(processed))
        if status not in (0,1):raise RuntimeError('Native DP failure '+str(status))
        best=route
        if length.value:
            decoded=[]
            for state in order[:length.value]:
                path=components[int(state)//2]
                decoded.extend(path if int(state)%2==0 else path[::-1])
            cost=validate(d,decoded)
            assert abs(cost-constant-value.value)<1e-6
            edges={tuple(sorted((decoded[i-1],v))) for i,v in enumerate(decoded)}
            assert (edges^current)<=free
            assert all(v in released or ((v in decoded)==(v in route)) for v in range(len(d['visits'])))
            if cost<route_cost(d,route):best=decoded
        return self.result(route,best,released,free,status,began,q,processed.value,'optimal' if status==0 else 'time_limit')

    def result(self,original,best,released,free,status,began,q,processed,reason):
        return dict(cost=validate(self.d,best),route=[int(self.d['visits'][v]) for v in best],
                    status=status,proven_optimal=status==0,reason=reason,contracted_nodes=q,
                    dp_processed_states=processed,total_seconds=time.perf_counter()-began,
                    free_x=len(free),free_y=len(released-self.d['mandatory']))
