"""Ten-seed paired TRAIN study; launch blocked until full-reference calibration passes."""
import fcntl,hashlib,json,os,random,sys,time
from contextlib import contextmanager
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from experiments.lkh.paired70_reference import run as reference
from experiments.lkh.paired70_hybrid import run as candidate
from experiments.lkh.paired70_candidate_audit import audit_candidate
from experiments.lkh.paired70_reference_audit import audit_reference
from experiments.lkh.raw_route_audit import load
from experiments.lkh.awake300_train import ClockGuard,sleep_events_since
OUT=HERE/'paired70_rebuilt_train'

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

@contextmanager
def study_lock(path):
    """Hold a kernel lock for the entire study, including resume checks.

    Keep the inode after release: deleting a lock file can let another process
    lock a replacement inode while an existing worker still owns the old one.
    """
    with Path(path).open('a+') as handle:
        try:fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError('Paired70 study already running; no optimizer launched') from None
        try:
            handle.seek(0);handle.truncate()
            json.dump(dict(pid=os.getpid(),started_unix=time.time()),handle)
            handle.flush()
            yield
        finally:fcntl.flock(handle,fcntl.LOCK_UN)

def main():
    with study_lock(HERE/'paired70_train.queue.lock'):
        collect()

def collect():
    guard=HERE/'published70_calibration_audit.json'
    assert guard.exists(), 'Wait for completed and audited seed80 calibration'
    calibration=json.loads(guard.read_text())
    assert calibration['status']=='PASS' and calibration['all_references_completed']
    for name,source in (('verify_paired70_lock.json','paired70_train.py'),
                        ('verify_paired70_deadlines.json','paired70_hybrid.py'),
                        ('verify_paired70_candidate_audit.json','paired70_candidate_audit.py'),
                        ('verify_paired70_reference_audit.json','paired70_reference_audit.py')):
        verification=json.loads((HERE/name).read_text())
        assert verification['status']=='PASS' and verification['source_sha256']==digest(HERE/source)
    integration=json.loads((HERE/'verify_paired70_adapters.json').read_text())
    assert integration['status']=='PASS'
    for filename,hashed in integration['source_hashes'].items():assert digest(HERE/filename)==hashed
    verification=json.loads((HERE/'verify_paired70_batch_audit.json').read_text())
    assert verification['status']=='PASS'
    for filename,hashed in verification['source_hashes'].items():assert digest(HERE/filename)==hashed
    orchestration=json.loads((HERE/'verify_paired70_orchestration.json').read_text())
    assert orchestration['status']=='PASS'
    for filename,hashed in orchestration['source_hashes'].items():assert digest(HERE/filename)==hashed
    for cert in ('verify_paired70_rebuilt_orchestration.json', 'verify_paired70_rebuilt_lock.json'):
        verified=json.loads((HERE/cert).read_text())
        assert verified['status']=='PASS'
        bound=verified.get('source_sha256') or verified['source_hashes']['paired70_rebuilt_train.py']
        assert bound==digest(Path(__file__))
    old=json.loads((ROOT/'experiments/original_mdm/hard_training_baseline/summary.json').read_text())['rows']
    roles=json.loads((ROOT/'experiments/ml_lns/decoder_initialized_states/protocol.json').read_text())['roles']
    cases=[dict(instance=r['instance'],geometry=r['geometry']) for r in old]
    assert len(cases)==9 and all(roles[c['geometry']]=='train' for c in cases)
    schedule=[dict(instance=c['instance'],seed=s) for c in cases for s in range(81,91)]
    random.Random(20261020).shuffle(schedule)
    paths=[Path(__file__),HERE/'paired70_reference.py',HERE/'paired70_hybrid.py',
        HERE/'paired70_candidate_audit.py',HERE/'paired70_reference_audit.py',
        HERE/'audit_paired70_train.py',HERE/'summarize_paired70_train.py',HERE/'raw_route_audit.py',HERE/'awake300_train.py',
        HERE/'ruin_repair_absolute.py',HERE/'feasible_fallback.py',HERE/'screen.py',
        ROOT/'experiments/original_mdm/timed_experiment.py',ROOT/'experiments/ml_lns/pilot.py',
        ROOT/'experiments/ml_branching/oracle_probe.py',ROOT/'research_notes/paired70_ten_seed_design.md']
    paths += [HERE/'audit_paired70_rebuilt.py', HERE/'verify_paired70_rebuilt_orchestration.json', HERE/'verify_paired70_rebuilt_lock.json', HERE/'paired70_train.py', ROOT/'research_notes/paired70_rebuilt_design.md', HERE/'rebuilt_backend/provenance.json']
    rebuilt = json.loads((HERE/'rebuilt_backend/provenance.json').read_text())
    expected_rebuilt = '1908f80deeafc8952acb21d36ba47ed618a71f2c4527c9654338cd24316be949'
    assert rebuilt['actual_sha256'] == expected_rebuilt
    assert digest(HERE/'rebuilt_backend/LKH') == expected_rebuilt
    assert digest(Path('/private/tmp/pcctp-lkh3-research/LKH')) == expected_rebuilt
    build=ROOT/'experiments/original_mdm/build_timed_elite_bounded'
    binaries=[build/'original_mdm',build/'FPMax/fim_maximal',build/'FPMax/fpmax_hnmp',Path('/private/tmp/pcctp-lkh3-research/LKH')]
    calibration_protocol=json.loads((HERE/'published70_calibration/protocol.json').read_text())
    for filename,hashed in calibration_protocol['binary_hashes'].items():
        assert digest(Path(filename))==hashed, 'Reference binary changed since calibration'
    protocol=dict(scope='Separate rebuilt-LKH restart of all90 TRAIN pairs; old five pairs excluded; not independent confirmation',
        reconstruction_provenance_sha256=digest(HERE/'rebuilt_backend/provenance.json'),
        superseded_incomplete_study='paired70_train',
        cases=cases,seeds=list(range(81,91)),schedule=schedule,scheduled_pairs=90,
        planned_executions_if_all_references_complete=270,
        reference_iterations=70,reference_safety_timeout=43200,environment_guarded=True,
        candidates=['hybrid10','hybrid60'],budget_rule='per-pair completed reference process wall duration',
        calibration_audit_sha256=digest(guard),
        source_hashes={str(p.relative_to(ROOT)):digest(p) for p in paths},
        binary_hashes={str(p):digest(p) for p in binaries})
    OUT.mkdir(exist_ok=True);refs=OUT/'references';candidates=OUT/'candidates'
    refs.mkdir(exist_ok=True);candidates.mkdir(exist_ok=True)
    p=OUT/'protocol.json'
    if p.exists():assert json.loads(p.read_text())==protocol
    else:p.write_text(json.dumps(protocol,indent=2)+'\n')
    protocol_hash=digest(p)
    execute_pairs(schedule,OUT,protocol_hash)


def execute_pairs(schedule,out,protocol_hash):
    out=Path(out);refs=out/'references';candidates=out/'candidates'
    refs.mkdir(parents=True,exist_ok=True);candidates.mkdir(exist_ok=True)
    for pair in schedule:
        name,seed=pair['instance'],pair['seed'];target=out/f'{name}.seed{seed}.pair.json'
        if target.exists():
            cached=json.loads(target.read_text());assert cached['protocol_sha256']==protocol_hash
            envpath=out/cached['environment_file']
            assert digest(envpath)==cached['environment_sha256']
            assert json.loads(envpath.read_text())['status']=='PASS'
            continue
        envpath=out/f'environment.{name}.seed{seed}.json'
        assert not envpath.exists(),'Preserve and review an incomplete or compromised pair before resuming'
        assert not (refs/f'{name}.seed{seed}.reference.json').exists(),'Unfinished pair requires review before reuse'
        with ClockGuard() as guard:
            r=reference(name,seed,refs);refpath=refs/f'{name}.seed{seed}.reference.json'
            data=load(name)
            audit_reference(r,data,refs,seed)
            record=dict(instance=name,seed=seed,protocol_sha256=protocol_hash,
                reference_file=str(refpath.relative_to(out)),reference_sha256=digest(refpath),
                status='REFERENCE_INCOMPLETE',candidate_records={},unavailable_candidate_modes=[])
            if r['status']=='COMPLETED70':
                budget=r['reference_duration_seconds']
                order=['hybrid10','hybrid60'];random.Random(name+str(seed)+'paired70').shuffle(order)
                for mode in order:
                    result=candidate(name,seed,mode,budget,candidates)
                    audit_candidate(result,data,budget)
                    source=candidates/f'{name}.seed{seed}.{mode}.json'
                    record['candidate_records'][mode]=dict(file=str(source.relative_to(out)),
                        sha256=digest(source),cost=result['cost'],feasible=result['feasible'])
                record.update(status='PAIR_COMPLETED' if all(c['feasible'] for c in record['candidate_records'].values()) else 'CANDIDATE_FAILURE',
                              budget_seconds=budget,candidate_order=order)
            else:record['unavailable_candidate_modes']=['hybrid10','hybrid60']
        environmental=dict(status='PASS',started_unix=guard.started,finished_unix=time.time(),
            max_monitor_gap_seconds=guard.max_gap,max_clock_difference_seconds=guard.clock_difference,
            sleep_events=sleep_events_since(guard.started))
        if guard.max_gap>5 or guard.clock_difference>3 or environmental['sleep_events']:
            environmental['status']='ENVIRONMENT_FAILURE'
        envpath.write_text(json.dumps(environmental,indent=2)+'\n')
        assert environmental['status']=='PASS','Compromised pair preserved; no subsequent optimizer call'
        record.update(environment_file=envpath.name,environment_sha256=digest(envpath))
        temp=target.with_suffix('.tmp');temp.write_text(json.dumps(record,indent=2)+'\n');temp.replace(target)
        print(json.dumps(record),flush=True)
    rows=[json.loads(p.read_text()) for p in out.glob('g_*.pair.json')]
    assert len(rows)==90
    (out/'completion.json').write_text(json.dumps(dict(scheduled_pairs=90,
        planned_executions=270,reference_records=90,
        candidate_records=sum(len(r['candidate_records']) for r in rows)))+'\n')


if __name__=='__main__':main()
