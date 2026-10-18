"""Verify portable native route descriptors against the frozen Python definition."""
import ctypes
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import time
import numpy as np
from operator_prediction import state_features
from experiments.ml_lns.pilot import read_instance

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]

def main():
    build=HERE/'operator_inference_build';build.mkdir(exist_ok=True)
    library=build/'state.dylib'
    subprocess.run([shutil.which('g++-15') or 'g++','-O3','-std=c++11','-dynamiclib',
                    str(HERE/'operator_state.cpp'),'-o',str(library)],check=True)
    fn=ctypes.CDLL(str(library)).pcctp_state_features
    double=ctypes.POINTER(ctypes.c_double);integer=ctypes.POINTER(ctypes.c_int);byte=ctypes.POINTER(ctypes.c_ubyte)
    fn.argtypes=[ctypes.c_int,ctypes.c_int,integer,byte,double,ctypes.c_double,double,
                 ctypes.c_int,byte,ctypes.c_double,ctypes.c_double,double];fn.restype=None
    ptr=lambda a,kind:a.ctypes.data_as(kind)
    rows=[];maximum=0.;checks=0
    for path in sorted((HERE/'move_profiles').glob('g_*.json')):
        s=json.loads(path.read_text())
        if 'records' not in s:continue
        d=read_instance(ROOT/'data'/s['instance']);mapping={int(v):i for i,v in enumerate(d['visits'])};n=len(mapping)
        mandatory=np.array([v in d['mandatory'] for v in range(n)],dtype=np.uint8)
        prize=np.ascontiguousarray(d['prizes'],dtype=float);cost=np.ascontiguousarray(d['costs'],dtype=float)
        cover=np.ascontiguousarray([[v in ns for v in range(n)] for ns in d['neighborhoods']],dtype=np.uint8)
        out=np.empty(10);count=0;elapsed=[];seen=set()
        for record in s['records']:
            key=tuple(record['before_route'])
            if key in seen:continue
            seen.add(key);route=np.array([mapping[v] for v in key],dtype=np.int32)
            args=(n,len(route),ptr(route,integer),ptr(mandatory,byte),ptr(prize,double),d['quota'],
                  ptr(cost,double),len(d['neighborhoods']),ptr(cover,byte),float(cost.mean()),float(prize.sum()),ptr(out,double))
            began=time.perf_counter();fn(*args);elapsed.append(time.perf_counter()-began)
            error=float(np.max(abs(out-state_features(d,route.tolist()))))
            assert error<1e-9;(checks:=checks+1);maximum=max(maximum,error);count+=1
            if count==20:break
        rows.append(dict(instance=s['instance'],states=count,median_call_seconds=float(np.median(elapsed))))
    result=dict(status='PASS',states=checks,max_absolute_error=maximum,rows=rows,
                scope='Native descriptor parity; timing excludes one-off graph preparation, which must be included online',
                code_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in
                             ['verify_operator_state.py','operator_state.cpp','operator_prediction.py']})
    (HERE/'operator_state_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))

if __name__=='__main__':main()
