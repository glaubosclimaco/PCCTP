"""Independent saved-route audit of fixed-iteration baseline convergence."""
import json,sys,hashlib
from pathlib import Path
from raw_route_audit import load,check
HERE=Path(__file__).resolve().parent
if __name__=='__main__':
    partial='--partial' in sys.argv;assert sys.argv[1:] in ([],['--partial'])
    folder=HERE/'fixed_work';complete=(folder/'completion.json').exists()
    assert partial or complete,'Batch still running'
    protocol=json.loads((folder/'protocol.json').read_text());rows=[];hashes={};checked=0
    for case in protocol['cases']:
        p=folder/f"{case['instance']}.json"
        if not p.exists():continue
        r=json.loads(p.read_text());native=r['native'];hashes[p.name]=hashlib.sha256(p.read_bytes()).hexdigest()
        assert native['budget_seconds']==0 and native['elite_count_bounded']
        d=load(r['instance']);assert native['input_sha256']==d['input_hash']
        if r['status']=='NO_VALID_EXPORT':rows.append(dict(instance=r['instance'],status=r['status']));continue
        check(d,r['native_route'],r['native_cost']);checked+=1
        lkh=r['lkh'];check(d,lkh['route'],lkh['cost']);checked+=1
        assert lkh['initial_cost']==r['native_cost'] and set(lkh['route'])==set(r['native_route']) and lkh['cost']<=r['native_cost']
        payload=None
        for source in (folder/f"{r['instance']}.native.raw",folder/f"{r['instance']}.native.raw.checkpoint"):
            if source.exists():
                try:payload=json.loads(source.read_text());break
                except json.JSONDecodeError:pass
        assert payload is not None
        for snap in payload['snapshots']:check(d,snap['route'],snap['cost']);checked+=1
        assert r['last_iteration']==max((s['iteration'] for s in payload['snapshots']),default=0)
        if r['status']=='COMPLETED70':assert r['last_iteration']==70 and payload['completed'] and not native['hard_timeout']
        else:assert r['status']=='CENSORED'
        cpu=sum(native[k] for k in ['parent_cpu_seconds','child_user_cpu_seconds','child_system_cpu_seconds'])
        rows.append(dict(instance=r['instance'],status=r['status'],last_iteration=r['last_iteration'],mdm_cost=r['native_cost'],lkh_cost=lkh['cost'],gain_percent=100*(r['native_cost']-lkh['cost'])/r['native_cost'],native_wall_seconds=native['process_wall_seconds'],cpu_seconds=cpu,published_best=case['published']['best']))
    if complete:assert len(rows)==9
    result=dict(status='PASS' if rows else 'NO_COMPLETED_CASES',batch_complete=complete,scope='fixed-work diagnostic; censored cases retained; LKH receives additional time',cases_available=len(rows),routes_checked=checked,rows=rows,result_hashes=hashes)
    (HERE/'fixed_work_audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ('rows','result_hashes')}))
