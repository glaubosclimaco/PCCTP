"""Audit strict-deadline results, preserve native clock records, expose overruns."""
import hashlib
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from experiments.ml_lns.pilot import read_instance,validate

def main():
    out=HERE/'timed_development';native_dir=out/'native_records';native_dir.mkdir(exist_ok=True)
    rows=[];initial={};cache={}
    invalid_path=HERE/'invalid_timing_protocols.json'
    invalid=json.loads(invalid_path.read_text()) if invalid_path.exists() else {}
    for path in sorted(out.glob('g_*.json')):
        r=json.loads(path.read_text());fingerprint=path.stem.split('.')[-1]
        if fingerprint in invalid:continue
        protocol=json.loads((out/f'protocol.{fingerprint}.json').read_text())
        assert any(c['instance']==r['instance'] and c['role']=='development' for c in protocol['cases'])
        assert r['operator_policy'] in protocol['methods'] and r['budget_seconds']==protocol['budget']
        assert r['method']=='timed_MDM_'+r['operator_policy']
        if r['instance'] not in cache:
            d=read_instance(ROOT/'data'/r['instance']);cache[r['instance']]=(d,{int(v):i for i,v in enumerate(d['visits'])})
        d,mapping=cache[r['instance']]
        assert hashlib.sha256((ROOT/'data'/r['instance']).read_bytes()).hexdigest()==r['input_sha256']
        previous=float('inf');elapsed=-1.
        for s in r['anytime']:
            assert abs(validate(d,[mapping[v] for v in s['route']])-s['cost'])<1e-6
            assert elapsed<=s['seconds']<=r['budget_seconds']+1e-9 and s['cost']<previous
            previous=s['cost'];elapsed=s['seconds']
        assert r['feasible_before_deadline']==bool(r['anytime'])
        assert r['cost']==(r['anytime'][-1]['cost'] if r['anytime'] else None)
        assert r['route']==(r['anytime'][-1]['route'] if r['anytime'] else [])
        if r['anytime']:
            first=r['anytime'][0];key=(r['instance'],r['seed']);value=(first['cost'],first['route'])
            if key in initial:assert initial[key]==value
            initial[key]=value
        raw=path.with_suffix('.raw');checkpoint=Path(str(raw)+'.checkpoint');saved=native_dir/path.name
        source=raw if raw.exists() else checkpoint if checkpoint.exists() else saved if saved.exists() else None
        native_elapsed=None
        if source:
            payload=json.loads(source.read_text())
            assert payload['deadline_cost']==r['cost'] and payload['anytime']==r['anytime']
            native_elapsed=payload['seconds']
            assert native_elapsed<=r['process_wall_seconds']+.01
            saved.write_text(json.dumps(payload,indent=2)+'\n')
        rows.append(dict(file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                         instance=r['instance'],geometry=r['geometry'],seed=r['seed'],policy=r['operator_policy'],
                         cost=r['cost'],feasible=r['feasible_before_deadline'],hard_timeout=r['hard_timeout'],
                         exit_code=r['exit_code'],parent_wall_seconds=r['process_wall_seconds'],
                         native_elapsed_seconds=native_elapsed,
                         native_overrun_seconds=max(0.,native_elapsed-r['budget_seconds']) if native_elapsed is not None else None,
                         native_record_sha256=hashlib.sha256(saved.read_bytes()).hexdigest() if saved.exists() else None))
    result=dict(status='PASS',files=len(rows),rows=rows,
                excluded_invalid_protocols=invalid,
                note='Snapshot of completed files; parent wall includes subprocess timeout-polling latency; native elapsed includes debited launch time. Overruns are not counted as search improvement.')
    (HERE/'timed_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='rows'}))

if __name__=='__main__':main()
