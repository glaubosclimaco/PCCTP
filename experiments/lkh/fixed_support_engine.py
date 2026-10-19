"""Fixed-support LKH ablation adapter; generous trial cap, charged wall deadline."""
import time,json,tempfile,subprocess,hashlib
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from experiments.ml_lns.pilot import read_instance,validate
BIN=Path('/private/tmp/pcctp-lkh3-research/LKH')


def polish(d,route,began=None,seconds=1.,seed=42):
    began=time.perf_counter() if began is None else began
    old=validate(d,route);chosen=sorted(route);index={v:i+1 for i,v in enumerate(chosen)}
    with tempfile.TemporaryDirectory(dir='/private/tmp',prefix='pcctp-lkh-') as folder:
        folder=Path(folder)
        matrix=['NAME: PCCTP_SUPPORT','TYPE: TSP',f'DIMENSION: {len(chosen)}','EDGE_WEIGHT_TYPE: EXPLICIT','EDGE_WEIGHT_FORMAT: FULL_MATRIX','EDGE_WEIGHT_SECTION']
        for a in chosen:
            values=[d['costs'][a,b] for b in chosen]
            assert all(v==int(v) for v in values)
            matrix.append(' '.join(str(int(v)) for v in values))
        (folder/'instance.tsp').write_text('\n'.join(matrix+['EOF'])+'\n')
        (folder/'initial.tour').write_text('\n'.join(['NAME: INITIAL','TYPE: TOUR',f'DIMENSION: {len(chosen)}','TOUR_SECTION']+[str(index[v]) for v in route]+['-1','EOF'])+'\n')
        (folder/'run.par').write_text('PROBLEM_FILE = instance.tsp\nINITIAL_TOUR_FILE = initial.tour\nTOUR_FILE = result.tour\nRUNS = 1\nMAX_TRIALS = 1000000\nSEED = '+str(seed)+'\nTIME_LIMIT = '+str(.8*seconds)+'\nTRACE_LEVEL = 0\n')
        timeout=False;code=None;log='';best=route
        remaining=seconds-(time.perf_counter()-began)
        if remaining>0:
            try:
                p=subprocess.run([str(BIN),'run.par'],cwd=folder,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=remaining)
                code=p.returncode;log=p.stdout
            except subprocess.TimeoutExpired:timeout=True
        result=folder/'result.tour'
        if code==0 and result.exists() and time.perf_counter()-began<=seconds:
            tokens=result.read_text().split();start=tokens.index('TOUR_SECTION')+1
            order=[]
            for token in tokens[start:]:
                if token=='-1':break
                order.append(int(token)-1)
            assert len(order)==len(chosen) and set(order)==set(range(len(chosen)))
            candidate=[chosen[i] for i in order];cost=validate(d,candidate)
            assert set(candidate)==set(route)
            if time.perf_counter()-began<=seconds and cost<old:best=candidate
        elapsed=time.perf_counter()-began
    return dict(initial_cost=old,cost=validate(d,best),improved=validate(d,best)<old,
                seconds=elapsed,timeout=timeout,exit_code=code,log=log,
                route=[int(d['visits'][v]) for v in best])


