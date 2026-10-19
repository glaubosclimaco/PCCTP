"""Operational wall/CPU diagnostic; no algorithm or performance claim."""
import sys,time,json,resource,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from experiments.original_mdm.timed_experiment import run_one
HERE=Path(__file__).resolve().parent

def calibration():
    wall=time.monotonic();cpu=time.process_time();iterations=0;x=b'PCCTP calibration'
    while time.process_time()-cpu<.2:
        for _ in range(1000):x=hashlib.sha256(x).digest()
        iterations+=1000
    return dict(wall_seconds=time.monotonic()-wall,cpu_seconds=time.process_time()-cpu,iterations=iterations)

if __name__=='__main__':
    out=HERE/'timing_diagnostic';out.mkdir(exist_ok=True)
    names=['g_a280_VT56_T168_W56_50','g_rd400_VT80_T240_W80_25','g_ts225_VT45_T135_W45_50']
    for name in names:
        p=out/f'{name}.json'
        if p.exists():continue
        before=calibration();children=resource.getrusage(resource.RUSAGE_CHILDREN);began=time.monotonic();parent_cpu=time.process_time()
        native=run_one(name,42,'original',15.,out/f'{name}.native.json',True,parent_started=began)
        wall=time.monotonic()-began;parent=time.process_time()-parent_cpu;after=resource.getrusage(resource.RUSAGE_CHILDREN)
        r=dict(instance=name,scope='operational diagnostic only',wall_seconds=wall,parent_cpu_seconds=parent,child_user_cpu_seconds=after.ru_utime-children.ru_utime,child_system_cpu_seconds=after.ru_stime-children.ru_stime,calibration_before=before,calibration_after=calibration(),native=native)
        p.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items() if k!='native'}),flush=True)
    (out/'completion.json').write_text(json.dumps(dict(cases=3))+'\n')
