"""Full70 reference calibration only; no candidate or reserved instances."""
import hashlib,json,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from experiments.original_mdm.timed_experiment import run_one
from experiments.lkh.raw_route_audit import load,check
from experiments.lkh.awake300_train import ClockGuard,sleep_events_since

OUT=HERE/'published70_calibration'

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

if __name__=='__main__':
    assert (HERE/'ruin_total300_awake_train/completion.json').exists(), 'Wait for awake repetition queue'
    assert json.loads((HERE/'ruin_total300_awake_train/environment_audit.json').read_text())['status']=='PASS'
    previous=json.loads((HERE/'ruin_total300_awake_train_audit.json').read_text())
    assert previous['status']=='PASS' and previous['batch_complete']
    cases=json.loads((ROOT/'experiments/original_mdm/hard_training_baseline/summary.json').read_text())['rows']
    split=json.loads((ROOT/'experiments/ml_lns/decoder_initialized_states/protocol.json').read_text())
    assert len(cases)==9 and all(split['roles'][c['geometry']]=='train' for c in cases)
    priority={'g_gil262':0,'g_rd400':1}
    cases=sorted(cases,key=lambda c:(priority.get(c['geometry'],2),c['instance']))
    build=ROOT/'experiments/original_mdm/build_timed_elite_bounded'
    sources=[Path(__file__),HERE/'audit_published70_calibration.py',ROOT/'experiments/original_mdm/timed_experiment.py',
             HERE/'raw_route_audit.py',HERE/'awake300_train.py',ROOT/'research_notes/published70_calibration_protocol.md']
    binaries=[build/'original_mdm',build/'FPMax/fim_maximal',build/'FPMax/fpmax_hnmp']
    protocol=dict(scope='TRAIN reference calibration only, not ten-seed final replication',
        cases=[dict(instance=c['instance'],geometry=c['geometry']) for c in cases],
        seed=80,iterations=70,positive_budget_disabled=True,safety_timeout_seconds=43200,
        source_hashes={str(p.relative_to(ROOT)):digest(p) for p in sources},
        binary_hashes={str(p):digest(p) for p in binaries})
    OUT.mkdir(exist_ok=True);p=OUT/'protocol.json'
    assert not (OUT/'environment_failure.json').exists(),'Preserve and review compromised reference attempt before resuming'
    if p.exists():assert json.loads(p.read_text())==protocol
    else:p.write_text(json.dumps(protocol,indent=2)+'\n')
    for case in cases:
        name=case['instance'];target=OUT/f'{name}.seed80.json'
        if target.exists():
            prior=json.loads(target.read_text());envpath=OUT/prior['environment_file']
            assert digest(envpath)==prior['environment_sha256']
            assert json.loads(envpath.read_text())['status']=='PASS'
            continue
        assert not (OUT/f'{name}.seed80.native.json').exists(),'Review and preserve incomplete reference attempt before reuse'
        began=time.clock_gettime(time.CLOCK_MONOTONIC)
        with ClockGuard() as guard:
            native=run_one(name,80,'original',0.,OUT/f'{name}.seed80.native.json',True,
                           parent_started=began,hard_timeout_seconds=43200.)
        environmental=dict(status='PASS',instance=name,started_unix=guard.started,finished_unix=time.time(),
            max_monitor_gap_seconds=guard.max_gap,
            max_clock_difference_seconds=guard.clock_difference,sleep_events=sleep_events_since(guard.started))
        if guard.max_gap>5 or guard.clock_difference>3 or environmental['sleep_events']:
            environmental['status']='ENVIRONMENT_FAILURE'
            (OUT/'environment_failure.json').write_text(json.dumps(environmental,indent=2)+'\n')
            raise AssertionError('Reference environment compromised; raw attempt preserved')
        envpath=OUT/f'environment.{name}.seed80.json'
        envpath.write_text(json.dumps(environmental,indent=2)+'\n')
        raw=OUT/f'{name}.seed80.native.raw';exports=[];errors=[]
        for source in (raw,Path(str(raw)+'.checkpoint')):
            if not source.exists():continue
            try:payload=json.loads(source.read_text())
            except (UnicodeDecodeError,json.JSONDecodeError) as error:
                errors.append(dict(file=source.name,error=str(error)));continue
            exports.append((source,payload))
        d=load(name);checked=0
        for source,payload in exports:
            check(d,payload['route'],payload['cost']);checked+=1
            for snap in payload['snapshots']:
                check(d,snap['route'],snap['cost']);checked+=1
        chosen=next(((p,r) for p,r in exports if p.name==native['payload_source']),None)
        payload=chosen[1] if chosen else None
        complete=bool(payload and payload['completed'] and payload['iterations']==70
                      and not native['hard_timeout'] and native['exit_code']==0)
        if complete:assert max(s['iteration'] for s in payload['snapshots'])==70
        record=dict(instance=name,seed=80,status='COMPLETED70' if complete else 'CENSORED_OR_FAILED',
            environment_file=envpath.name,environment_sha256=digest(envpath),
            cost=payload['cost'] if complete else None,route=payload['route'] if complete else [],
            observed_cost=payload['cost'] if payload else None,
            observed_iterations=payload['iterations'] if payload else None,
            mining_calls=payload['mining_calls'] if payload else None,
            native=native,wrapper_wall_seconds=time.clock_gettime(time.CLOCK_MONOTONIC)-began,
            routes_checked=checked,input_sha256=d['input_hash'],payload_read_errors=errors,
            export_hashes={p.name:digest(p) for p,r in exports})
        tmp=target.with_suffix('.tmp');tmp.write_text(json.dumps(record,indent=2)+'\n');tmp.replace(target)
        print(json.dumps({k:v for k,v in record.items() if k not in ('native','route','export_hashes')}),flush=True)
    (OUT/'completion.json').write_text(json.dumps(dict(scheduled_runs=9))+'\n')
