"""Deadline-charged subprocess adapter for the verified symmetric 2-opt kernel."""
import subprocess,time
from experiments.ml_lns.pilot import validate
BIN='/private/tmp/pcctp-two-opt'

def two_opt(d,route,end):
    end=min(end,time.perf_counter()+.2)
    old=validate(d,route);n=len(route)
    values=[d['costs'][v,w] for v in route for w in route]
    assert all(v==int(v) for v in values)
    remaining=end-time.perf_counter()
    if remaining<=0:return route
    payload=f'{n} {remaining*.8}\n'+' '.join(str(int(v)) for v in values)+'\n'+' '.join(str(i) for i in range(n))+'\n'
    remaining=end-time.perf_counter()
    if remaining<=0:return route
    try:p=subprocess.run([BIN],input=payload,text=True,capture_output=True,timeout=remaining)
    except subprocess.TimeoutExpired:return route
    if p.returncode!=0:raise RuntimeError(p.stderr)
    t=p.stdout.split();candidate=[route[int(i)] for i in t[2:]]
    value=validate(d,candidate);assert value==int(t[0]) and set(candidate)==set(route)
    return candidate if value<old and time.perf_counter()<=end else route
