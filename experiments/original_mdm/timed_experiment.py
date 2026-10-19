"""Development matched total-time pilot with strict pre-deadline incumbents.

Parent launch time is debited in the native clock. Cooperative overruns and hard
timeouts are recorded; improvements reported after the deadline are excluded.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import random
import resource
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from experiments.ml_lns.pilot import read_instance,validate

def run_one(name,seed,policy,budget,target,elite_count_bounded=False,parent_started=None,hard_timeout_seconds=None):
    build=HERE/('build_timed_elite_bounded' if elite_count_bounded else 'build_timed')
    working=HERE/'runs'/f'{name}.seed{seed}.timed.{policy}'
    if elite_count_bounded:working=working.with_name(working.name+'.elite_bounded')
    fp=working/'FPMax';fp.mkdir(parents=True,exist_ok=True)
    for binary in ['fim_maximal','fpmax_hnmp']:
        link=fp/binary
        if not link.exists():link.symlink_to((build/'FPMax'/binary).resolve())
    raw=target.with_suffix('.raw');checkpoint=Path(str(raw)+'.checkpoint')
    raw.unlink(missing_ok=True);checkpoint.unlink(missing_ok=True)
    env=os.environ.copy();env['PCCTP_OPERATOR_POLICY']=policy;env['PCCTP_BUDGET']=str(budget)
    env.pop('PCCTP_SCHEDULE_VERIFY',None)
    timed_out=False;code=None
    usage_before=resource.getrusage(resource.RUSAGE_CHILDREN)
    parent_cpu_before=time.process_time()
    with target.with_suffix('.log').open('w') as log:
        began=time.clock_gettime(time.CLOCK_MONOTONIC) if parent_started is None else parent_started
        env['PCCTP_STARTED_MONOTONIC']=format(began,'.17g')
        try:
            process=subprocess.run([str((build/'original_mdm').resolve()),str((ROOT/'data'/name).resolve()),
                                    str(seed),'2','5',str(raw.resolve())],cwd=working,env=env,
                                   stdout=log,stderr=log,timeout=hard_timeout_seconds if hard_timeout_seconds is not None else budget+15 if budget>0 else 60)
            code=process.returncode
        except subprocess.TimeoutExpired:timed_out=True
        wall=time.clock_gettime(time.CLOCK_MONOTONIC)-began
    usage_after=resource.getrusage(resource.RUSAGE_CHILDREN)
    parent_cpu_seconds=time.process_time()-parent_cpu_before
    payload=None;payload_source=None;payload_read_errors=[]
    for source in (raw,checkpoint):
        if not source.exists():continue
        try:payload=json.loads(source.read_text());payload_source=source.name;break
        except (json.JSONDecodeError,UnicodeDecodeError) as error:
            payload_read_errors.append(dict(file=source.name,error=str(error)))
    d=read_instance(ROOT/'data'/name);mapping={int(v):i for i,v in enumerate(d['visits'])}
    trace=[]
    if payload:
        for s in [payload,*payload['snapshots'],*payload.get('anytime',[])]:
            assert abs(validate(d,[mapping[v] for v in s['route']])-s['cost'])<1e-6
        trace=payload['anytime'];last=float('inf');elapsed=-1.
        for s in trace:
            assert elapsed<=s['seconds']<=budget+1e-9 and s['cost']<last
            elapsed=s['seconds'];last=s['cost']
        assert bool(trace)==payload['feasible_before_deadline']
        assert payload['deadline_cost']==(trace[-1]['cost'] if trace else None)
        if trace:assert payload['deadline_route']==trace[-1]['route']
    row=dict(instance=name,geometry=d['geometry'],seed=seed,operator_policy=policy,
             method='timed_MDM_'+policy,budget_seconds=budget,process_wall_seconds=wall,
             parent_cpu_seconds=parent_cpu_seconds,
             child_user_cpu_seconds=usage_after.ru_utime-usage_before.ru_utime,
             child_system_cpu_seconds=usage_after.ru_stime-usage_before.ru_stime,
             elite_count_bounded=elite_count_bounded,
             overrun_seconds=max(0.,wall-budget),hard_timeout=timed_out,exit_code=code,
             native_returned=bool(payload and payload['completed']),feasible_before_deadline=bool(trace),
             payload_source=payload_source,payload_read_errors=payload_read_errors,
             cost=trace[-1]['cost'] if trace else None,route=trace[-1]['route'] if trace else [],anytime=trace,
             native_cost_after_return=payload['cost'] if payload else None,
             native_elapsed_seconds=payload['seconds'] if payload else None,
             clock_source='clock_gettime(CLOCK_MONOTONIC)',
             input_sha256=hashlib.sha256((ROOT/'data'/name).read_bytes()).hexdigest(),
             startup_debit_seconds=payload.get('startup_debit_seconds') if payload else None)
    temp=target.with_suffix('.partial');temp.write_text(json.dumps(row,indent=2)+'\n');temp.replace(target)
    return row

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--budget',type=float,default=1.)
    p.add_argument('--limit',type=int,default=9);p.add_argument('--seeds',type=int,nargs='+',default=[42,43,44])
    p.add_argument('--bound-elite-size',action='store_true')
    args=p.parse_args();assert args.budget>0
    out=HERE/'timed_development';out.mkdir(exist_ok=True)
    cases=[];seen=set()
    for path in sorted((ROOT/'experiments/ml_lns/decoder_initialized_states').glob('g_*.json')):
        b=json.loads(path.read_text())
        if b['role']!='development' or b['geometry'] in seen:continue
        seen.add(b['geometry']);cases.append(b)
        if len(cases)==args.limit:break
    methods=['original','static','size','generic','joint','adaptive']
    build=HERE/('build_timed_elite_bounded' if args.bound_elite_size else 'build_timed')
    protocol=dict(scope=__doc__,budget=args.budget,seeds=args.seeds,methods=methods,elite_count_bounded=args.bound_elite_size,
                  cases=[dict(instance=b['instance'],geometry=b['geometry'],role='development') for b in cases],
                  code_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in
                               ['timed_experiment.py','build.py','research_hook.cpp','research_hook.hpp','scheduler.hpp']},
                  build_manifest=json.loads((build/'manifest.json').read_text()),
                  binary_sha256=hashlib.sha256((build/'original_mdm').read_bytes()).hexdigest())
    fingerprint=hashlib.sha256(json.dumps(protocol,sort_keys=True).encode()).hexdigest()[:12]
    (out/f'protocol.{fingerprint}.json').write_text(json.dumps(protocol,indent=2)+'\n')
    # No-budget original must still reproduce the supplied fixed-iteration solver.
    control=run_one(cases[0]['instance'],42,'original',0.,out/f'fidelity.{fingerprint}.json',args.bound_elite_size)
    assert control['native_returned']
    if not args.bound_elite_size:assert control['native_cost_after_return']==cases[0]['cost']
    for b in cases:
        for seed in args.seeds:
            order=methods[:];random.Random(b['instance']+str(seed)).shuffle(order);first=[]
            for policy in order:
                target=out/f"{b['instance']}.seed{seed}.{policy}.{fingerprint}.json"
                row=json.loads(target.read_text()) if target.exists() else run_one(b['instance'],seed,policy,args.budget,target,args.bound_elite_size)
                if row['anytime']:first.append((row['anytime'][0]['cost'],row['anytime'][0]['route']))
                print(json.dumps({k:v for k,v in row.items() if k not in ('route','anytime')}),flush=True)
            if first:assert all(r==first[0] for r in first),'Unequal initial construction across policies'
    (out/f'completion.{fingerprint}.json').write_text(json.dumps(dict(cases=len(cases),seeds=args.seeds,methods=methods))+'\n')

if __name__=='__main__':main()
