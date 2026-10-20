"""Independent final-iteration, scalar route and provenance audit of reference runs."""
import hashlib,json
from pathlib import Path
from raw_route_audit import load,check
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

if __name__=='__main__':
    folder=HERE/'published70_calibration';p=folder/'protocol.json'
    if not p.exists():print('NO_RESULTS');raise SystemExit(0)
    protocol=json.loads(p.read_text())
    for n,h in protocol['source_hashes'].items():assert digest(ROOT/n)==h
    for n,h in protocol['binary_hashes'].items():assert digest(Path(n))==h
    rows=[];checked=0;keys=set()
    for p in sorted(folder.glob('g_*.seed80.json')):
        r=json.loads(p.read_text());assert r['seed']==80 and r['instance'] not in keys
        assert r['instance'] in {c['instance'] for c in protocol['cases']}
        keys.add(r['instance']);d=load(r['instance']);assert d['input_hash']==r['input_sha256']
        native=r['native'];assert native['budget_seconds']==0 and native['elite_count_bounded']
        payloads={}
        for name,hashed in r['export_hashes'].items():
            source=folder/name;assert digest(source)==hashed;payload=json.loads(source.read_text());payloads[name]=payload
            check(d,payload['route'],payload['cost']);checked+=1
            assert payload['budget_seconds']==0 and payload['seed']==80
            for snap in payload['snapshots']:check(d,snap['route'],snap['cost']);checked+=1
        chosen=payloads.get(native['payload_source'])
        complete=bool(chosen and chosen['completed'] and chosen['iterations']==70
                      and not native['hard_timeout'] and native['exit_code']==0)
        assert (r['status']=='COMPLETED70')==complete
        if complete:
            assert max(s['iteration'] for s in chosen['snapshots'])==70
            assert (r['route'],r['cost'])==(chosen['route'],chosen['cost'])
            check(d,r['route'],r['cost']);checked+=1
        else:assert r['cost'] is None and r['route']==[]
        rows.append(dict(instance=r['instance'],status=r['status'],cost=r['cost'],
            iterations=r['observed_iterations'],mining_calls=r['mining_calls'],
            native_process_wall_seconds=native['process_wall_seconds'],
            wrapper_wall_seconds=r['wrapper_wall_seconds'],result_sha256=digest(p)))
    finished=(folder/'completion.json').exists()
    if finished:assert len(rows)==9 and keys=={c['instance'] for c in protocol['cases']}
    report=dict(status='PASS' if rows else 'NO_RESULTS',batch_complete=finished,
        runs=len(rows),completed70=sum(r['status']=='COMPLETED70' for r in rows),
        censored_or_failed=sum(r['status']!='COMPLETED70' for r in rows),routes_checked=checked,
        all_references_completed=finished and all(r['status']=='COMPLETED70' for r in rows),
        rows=rows,scope='one-seed TRAIN reference calibration, not final ten-seed comparison')
    (HERE/'published70_calibration_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='rows'}))
