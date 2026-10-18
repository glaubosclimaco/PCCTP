"""One-second fixed-support LKH screening, not an equal-total-time comparison."""
import time,json,tempfile,subprocess,hashlib
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from experiments.ml_lns.pilot import read_instance,validate
BIN=Path('/private/tmp/pcctp-lkh3-research/LKH')


def polish(d,route,began=None,seconds=1.):
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
        (folder/'run.par').write_text('PROBLEM_FILE = instance.tsp\nINITIAL_TOUR_FILE = initial.tour\nTOUR_FILE = result.tour\nRUNS = 1\nMAX_TRIALS = 100\nSEED = 42\nTIME_LIMIT = '+str(.8*seconds)+'\nTRACE_LEVEL = 0\n')
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


if __name__=='__main__':
    base=ROOT/'experiments/original_mdm/hard_training_baseline'
    cases=json.loads((base/'summary.json').read_text())['rows'];assert len(cases)==9
    result=dict(binary_sha256=hashlib.sha256(BIN.read_bytes()).hexdigest(),rows=[],additional_time_screen=True)
    for case in cases:
        began=time.perf_counter();d=read_instance(ROOT/'data'/case['instance'])
        incumbent=json.loads((base/f"{case['instance']}.seed42.original.json").read_text())
        mapping={int(v):i for i,v in enumerate(d['visits'])};route=[mapping[v] for v in incumbent['route']]
        row=polish(d,route,began);row.update(instance=case['instance'],geometry=case['geometry'])
        result['rows'].append(row);result['improvements']=sum(r['improved'] for r in result['rows'])
        result['gate_passed']=len(result['rows'])==9 and result['improvements']>=3
        Path(__file__).with_name('screen_results.json').write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps({k:v for k,v in row.items() if k not in ('route','log')}),flush=True)
