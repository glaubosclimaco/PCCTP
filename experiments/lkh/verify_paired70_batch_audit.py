"""Synthetic90-pair audit integration, with no solver or benchmark data."""
import copy,hashlib,json,random,tempfile
from pathlib import Path
import audit_paired70_train as batch
from summarize_paired70_train import summarize

def write(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2)+'\n')
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

if __name__=='__main__':
    data=dict(n=4,types=[1,0,0,0],prizes=[1]*4,quota=1,radius=0,
        matrix=[0 if i==j else 1 for i in range(4) for j in range(4)],input_hash='synthetic')
    old_load=batch.load;batch.load=lambda _:data
    rejected=[]
    try:
        with tempfile.TemporaryDirectory(prefix='pcctp-batch-audit-',dir='/private/tmp') as tmp:
            root=Path(tmp);here=root/'observers';folder=here/'paired70_train'
            cases=[dict(instance='g_synthetic'+str(i),geometry='synthetic'+str(i)) for i in range(9)]
            write(root/'experiments/ml_lns/decoder_initialized_states/protocol.json',
                  dict(roles={c['geometry']:'train' for c in cases}))
            write(here/'published70_calibration_audit.json',dict(status='PASS'))
            schedule=[dict(instance=c['instance'],seed=s) for c in cases for s in range(81,91)]
            random.Random(20261020).shuffle(schedule)
            protocol=dict(cases=cases,seeds=list(range(81,91)),scheduled_pairs=90,schedule=schedule,
                source_hashes={},binary_hashes={},environment_guarded=True,
                calibration_audit_sha256=digest(here/'published70_calibration_audit.json'))
            write(folder/'protocol.json',protocol);protocol_hash=digest(folder/'protocol.json')
            for pair in schedule:
                name,seed=pair['instance'],pair['seed']
                raw=folder/'references'/f'{name}.seed{seed}.raw'
                payload=dict(seed=seed,budget_seconds=0,completed=True,iterations=70,
                    route=[0,1,2,3],cost=4,mining_calls=1,snapshots=[
                        dict(iteration=70,seconds=1.8,route=[0,1,2,3],cost=4)])
                write(raw,payload)
                reference=dict(instance=name,seed=seed,input_sha256='synthetic',status='COMPLETED70',
                    native=dict(budget_seconds=0,elite_count_bounded=True,payload_source=raw.name,
                        hard_timeout=False,exit_code=0,process_wall_seconds=2.),
                    observed_cost=4,observed_iterations=70,mining_calls=1,cost=4,route=[0,1,2,3],
                    reference_duration_seconds=2.,wrapper_wall_seconds=2.1,payload_sha256=digest(raw))
                refpath=raw.with_suffix('.reference.json');write(refpath,reference)
                order=['hybrid10','hybrid60'];random.Random(name+str(seed)+'paired70').shuffle(order)
                record=dict(instance=name,seed=seed,protocol_sha256=protocol_hash,
                    reference_file=str(refpath.relative_to(folder)),reference_sha256=digest(refpath),
                    budget_seconds=2.,candidate_order=order,candidate_records={},
                    unavailable_candidate_modes=[],status='PAIR_COMPLETED')
                for mode in order:
                    native=dict(budget_seconds=2.,feasible_before_deadline=True,cost=3,route=[0,1,2],
                                anytime=[dict(seconds=.5,cost=3,route=[0,1,2])])
                    candidate=dict(instance=name,seed=seed,mode=mode,budget_seconds=2.,
                        initialization_budget_seconds=2.,native_started=100.,wall_seconds=2.1,
                        fallback=dict(seconds=.1,cost=4,route=[0,1,2,3]),native=native,
                        pre_cost=3,pre_route=[0,1,2],post_offset=2.1,post=None,
                        cost=3,route=[0,1,2],origin='native',feasible=True)
                    path=folder/'candidates'/f'{name}.seed{seed}.{mode}.json';write(path,candidate)
                    write(path.with_name(path.stem+'.native.json'),native)
                    record['candidate_records'][mode]=dict(file=str(path.relative_to(folder)),
                        sha256=digest(path),cost=3,feasible=True)
                envpath=folder/f'environment.{name}.seed{seed}.json'
                write(envpath,dict(status='PASS',sleep_events=[],started_unix=100.,finished_unix=107.,
                    max_monitor_gap_seconds=1.1,max_clock_difference_seconds=.01))
                record.update(environment_file=envpath.name,environment_sha256=digest(envpath))
                write(folder/f'{name}.seed{seed}.pair.json',record)
            completion=dict(scheduled_pairs=90,planned_executions=270,reference_records=90,candidate_records=180)
            write(folder/'completion.json',completion)
            report=batch.audit_batch(folder,root,here)
            assert report['eligible_complete_study'] and report['pairs']==90
            summary=summarize(report,protocol)
            assert all(c['gate_passed'] and c['wins']==9 and c['mean_gain_percent']==25.
                       for c in summary['comparisons'].values())
            target=next(folder.glob('g_*.pair.json'));original=json.loads(target.read_text())
            for label,change in (
                ('wrong_pair_duration',lambda r:r.update(budget_seconds=3.)),
                ('wrong_reference_hash',lambda r:r.update(reference_sha256='wrong')),
                ('wrong_environment_hash',lambda r:r.update(environment_sha256='wrong')),
                ('missing_candidate',lambda r:r['candidate_records'].pop('hybrid10')),
                ('false_candidate_cost',lambda r:r['candidate_records']['hybrid10'].update(cost=2)),
                ('false_completion_status',lambda r:r.update(status='REFERENCE_INCOMPLETE'))):
                mutated=copy.deepcopy(original);change(mutated);write(target,mutated)
                try:batch.audit_batch(folder,root,here)
                except AssertionError:rejected.append(label)
                else:raise AssertionError('Corruption accepted: '+label)
            write(target,original)
            refpath=folder/original['reference_file'];ref=json.loads(refpath.read_text())
            raw=refpath.parent/ref['native']['payload_source'];payload=json.loads(raw.read_text())
            payload.update(completed=False,iterations=35);payload['snapshots'][0]['iteration']=35
            write(raw,payload)
            ref.update(status='REFERENCE_INCOMPLETE',cost=None,route=[],reference_duration_seconds=None,
                       observed_iterations=35,payload_sha256=digest(raw));write(refpath,ref)
            original.update(status='REFERENCE_INCOMPLETE',reference_sha256=digest(refpath),
                candidate_records={},unavailable_candidate_modes=['hybrid10','hybrid60'])
            del original['budget_seconds'];del original['candidate_order'];write(target,original)
            completion['candidate_records']=178;write(folder/'completion.json',completion)
            censored=batch.audit_batch(folder,root,here)
            assert censored['reference_incomplete']==1 and not censored['eligible_complete_study']
            summary=summarize(censored,protocol)
            assert summary['completed_geometries']==8
            assert all(not c['gate_passed'] for c in summary['comparisons'].values())
    finally:batch.load=old_load
    report=dict(status='PASS',synthetic_pairs=90,synthetic_candidate_records=180,
        valid_cases=['complete90pairs','censored_reference_keeps_two_unavailable_slots'],
        rejected_corruptions=rejected,source_hashes={n:digest(Path(__file__).with_name(n)) for n in
            ('audit_paired70_train.py','paired70_reference_audit.py','paired70_candidate_audit.py',
             'summarize_paired70_train.py')},
        scope='synthetic audit integration only; no optimizer or PCCTP benchmark execution')
    write(Path(__file__).with_suffix('.json'),report);print(json.dumps(report))
