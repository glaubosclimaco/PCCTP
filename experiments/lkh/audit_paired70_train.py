"""Independent full-batch audit of reference-defined paired candidates."""
import hashlib,json,math,random
from pathlib import Path
from raw_route_audit import load
from paired70_reference_audit import audit_reference
from paired70_candidate_audit import audit_candidate
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def audit_batch(folder,root=ROOT,here=HERE):
    folder=Path(folder);path=folder/'protocol.json'
    if not path.exists():return dict(status='NO_RESULTS')
    protocol=json.loads(path.read_text());protocol_hash=digest(path)
    assert protocol['seeds']==list(range(81,91)) and protocol['scheduled_pairs']==90
    for name,hashed in protocol['source_hashes'].items():assert digest(root/name)==hashed
    for name,hashed in protocol['binary_hashes'].items():assert digest(Path(name))==hashed
    assert digest(here/'published70_calibration_audit.json')==protocol['calibration_audit_sha256']
    roles=json.loads((root/'experiments/ml_lns/decoder_initialized_states/protocol.json').read_text())['roles']
    assert len(protocol['cases'])==9 and all(roles[c['geometry']]=='train' for c in protocol['cases'])
    expected={(c['instance'],s) for c in protocol['cases'] for s in range(81,91)}
    assert {(p['instance'],p['seed']) for p in protocol['schedule']}==expected
    assert len(protocol['schedule'])==90
    schedule=[dict(instance=c['instance'],seed=s) for c in protocol['cases'] for s in range(81,91)]
    random.Random(20261020).shuffle(schedule)
    assert protocol['schedule']==schedule
    rows=[];seen=set();checked=0
    for path in sorted(folder.glob('g_*.pair.json')):
        r=json.loads(path.read_text());key=(r['instance'],r['seed'])
        assert key in expected and key not in seen;seen.add(key)
        assert r['protocol_sha256']==protocol_hash
        if protocol.get('environment_guarded'):
            envpath=folder/r['environment_file'];assert digest(envpath)==r['environment_sha256']
            environment=json.loads(envpath.read_text())
            assert environment['status']=='PASS' and not environment['sleep_events']
            for name,limit in (('max_monitor_gap_seconds',5.),('max_clock_difference_seconds',3.)):
                assert math.isfinite(environment[name]) and 0<=environment[name]<=limit
            assert environment['started_unix']<=environment['finished_unix']
        source=folder/r['reference_file'];assert digest(source)==r['reference_sha256']
        reference=json.loads(source.read_text());data=load(r['instance'])
        assert reference['instance']==r['instance'] and reference['seed']==r['seed']
        result=audit_reference(reference,data,source.parent,r['seed']);checked+=result['routes_checked']
        costs={};feasible={};result_hashes={str(source.relative_to(folder)):digest(source)}
        if result['completed70']:
            budget=result['duration'];assert r['budget_seconds']==budget
            assert set(r['candidate_records'])=={'hybrid10','hybrid60'}
            assert r['unavailable_candidate_modes']==[]
            order=['hybrid10','hybrid60'];random.Random(r['instance']+str(r['seed'])+'paired70').shuffle(order)
            assert r['candidate_order']==order
            for mode,entry in r['candidate_records'].items():
                p=folder/entry['file'];assert digest(p)==entry['sha256']
                candidate=json.loads(p.read_text())
                assert (candidate['instance'],candidate['seed'],candidate['mode'])==(r['instance'],r['seed'],mode)
                verified=audit_candidate(candidate,data,budget);checked+=verified['routes_checked']
                assert (entry['cost'],entry['feasible'])==(candidate['cost'],candidate['feasible'])
                native_path=p.with_name(p.stem+'.native.json')
                assert json.loads(native_path.read_text())==candidate['native']
                costs[mode]=candidate['cost'];feasible[mode]=candidate['feasible']
                result_hashes[str(p.relative_to(folder))]=digest(p)
            assert r['status']==('PAIR_COMPLETED' if all(feasible.values()) else 'CANDIDATE_FAILURE')
        else:
            assert r['status']=='REFERENCE_INCOMPLETE' and r['candidate_records']=={}
            assert set(r['unavailable_candidate_modes'])=={'hybrid10','hybrid60'}
            assert 'budget_seconds' not in r and 'candidate_order' not in r
        rows.append(dict(instance=r['instance'],seed=r['seed'],status=r['status'],
            reference_cost=result['cost'],duration=result['duration'],candidate_costs=costs,
            candidate_feasibility=feasible,pair_sha256=digest(path),result_hashes=result_hashes))
    finished=(folder/'completion.json').exists()
    if finished:
        assert len(rows)==90 and seen==expected
        completion=json.loads((folder/'completion.json').read_text())
        assert completion['reference_records']==90
        assert completion['scheduled_pairs']==90 and completion['planned_executions']==270
        assert completion['candidate_records']==sum(len(r['candidate_costs']) for r in rows)
    report=dict(status='PASS' if rows else 'NO_RESULTS',batch_complete=finished,
        pairs=len(rows),planned_pairs=90,planned_executions=270,
        reference_incomplete=sum(r['status']=='REFERENCE_INCOMPLETE' for r in rows),
        candidate_failure_pairs=sum(r['status']=='CANDIDATE_FAILURE' for r in rows),
        eligible_complete_study=finished and all(r['status']=='PAIR_COMPLETED' for r in rows),
        routes_checked=checked,rows=rows,scope='TRAIN, completed70 baseline plus reference-defined candidate budget; not full107-instance replication')
    return report

if __name__=='__main__':
    report=audit_batch(HERE/'paired70_train')
    if report['status']!='NO_RESULTS':
        (HERE/'paired70_train_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='rows'}))
