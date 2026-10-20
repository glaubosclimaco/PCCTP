"""Completed native70 reference; a censored cost never becomes a comparison baseline."""
import hashlib,json,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from experiments.original_mdm.timed_experiment import run_one
from experiments.lkh.raw_route_audit import load,check


def run(name,seed,out):
    assert seed in range(81,91)
    out=Path(out);out.mkdir(exist_ok=True)
    target=out/f'{name}.seed{seed}.reference.json'
    if target.exists():
        record=json.loads(target.read_text())
        assert record['instance']==name and record['seed']==seed
        return record
    began=time.clock_gettime(time.CLOCK_MONOTONIC)
    native=run_one(name,seed,'original',0.,out/f'{name}.seed{seed}.reference.native.json',
        True,parent_started=began,hard_timeout_seconds=43200.)
    # Measure the native search process, not post-exit independent auditing.
    duration=native['process_wall_seconds']
    source=out/native['payload_source'] if native['payload_source'] else None
    payload=json.loads(source.read_text()) if source else None
    data=load(name);checked=0
    if payload:
        check(data,payload['route'],payload['cost']);checked+=1
        for event in payload['snapshots']:
            check(data,event['route'],event['cost']);checked+=1
    complete=bool(payload and payload['completed'] and payload['iterations']==70
                  and not native['hard_timeout'] and native['exit_code']==0)
    if complete:assert max(e['iteration'] for e in payload['snapshots'])==70
    record=dict(instance=name,seed=seed,status='COMPLETED70' if complete else 'REFERENCE_INCOMPLETE',
        cost=payload['cost'] if complete else None,route=payload['route'] if complete else [],
        reference_duration_seconds=duration if complete else None,
        observed_cost=payload['cost'] if payload else None,
        observed_iterations=payload['iterations'] if payload else None,
        mining_calls=payload['mining_calls'] if payload else None,
        native=native,wrapper_wall_seconds=time.clock_gettime(time.CLOCK_MONOTONIC)-began,
        routes_checked=checked,input_sha256=data['input_hash'],
        payload_sha256=hashlib.sha256(source.read_bytes()).hexdigest() if source else None)
    temp=target.with_suffix('.tmp');temp.write_text(json.dumps(record,indent=2)+'\n');temp.replace(target)
    return record
