"""Exercise the actual90-pair scheduler with mocked optimizers and real audits."""
import hashlib,json,random,tempfile
from pathlib import Path
import paired70_train as runner

def write(path,value):path.write_text(json.dumps(value)+'\n')
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

if __name__=='__main__':
    data=dict(n=4,types=[1,0,0,0],prizes=[1]*4,quota=1,radius=0,
        matrix=[0 if i==j else 1 for i in range(4) for j in range(4)],input_hash='synthetic')
    calls=[];budgets={}
    def reference(name,seed,folder):
        calls.append(('reference',name,seed));complete=not (name=='g_synthetic0' and seed==81)
        duration=2.+(seed-81)/100
        iterations=70 if complete else 35
        raw=folder/f'{name}.seed{seed}.raw'
        payload=dict(seed=seed,budget_seconds=0,completed=complete,iterations=iterations,
            route=[0,1,2,3],cost=4,mining_calls=1,
            snapshots=[dict(iteration=iterations,seconds=1.,route=[0,1,2,3],cost=4)])
        write(raw,payload)
        record=dict(instance=name,seed=seed,input_sha256='synthetic',
            status='COMPLETED70' if complete else 'REFERENCE_INCOMPLETE',
            native=dict(budget_seconds=0,elite_count_bounded=True,payload_source=raw.name,
                hard_timeout=False,exit_code=0,process_wall_seconds=duration),
            observed_cost=4,observed_iterations=iterations,mining_calls=1,
            cost=4 if complete else None,route=[0,1,2,3] if complete else [],
            reference_duration_seconds=duration if complete else None,
            wrapper_wall_seconds=duration+.1,payload_sha256=digest(raw))
        write(folder/f'{name}.seed{seed}.reference.json',record)
        budgets[(name,seed)]=duration if complete else None
        return record
    def candidate(name,seed,mode,budget,folder):
        # API receives only name, seed, mode and T: no reference route/cost.
        assert budgets[(name,seed)]==budget
        calls.append((mode,name,seed))
        native=dict(budget_seconds=budget,feasible_before_deadline=True,cost=3,route=[0,1,2],
            anytime=[dict(seconds=.5,cost=3,route=[0,1,2])])
        record=dict(instance=name,seed=seed,mode=mode,budget_seconds=budget,
            initialization_budget_seconds=budget,native_started=100.,wall_seconds=budget+.1,
            fallback=dict(seconds=.1,cost=4,route=[0,1,2,3]),native=native,
            pre_cost=3,pre_route=[0,1,2],post_offset=budget+.1,post=None,
            cost=3,route=[0,1,2],origin='native',feasible=True)
        path=folder/f'{name}.seed{seed}.{mode}.json';write(path,record)
        write(path.with_name(path.stem+'.native.json'),native)
        return record
    saved=runner.reference,runner.candidate,runner.load,runner.HERE
    try:
        runner.reference=reference;runner.candidate=candidate;runner.load=lambda _:data
        with tempfile.TemporaryDirectory(prefix='pcctp-orchestration-',dir='/private/tmp') as tmp:
            folder=Path(tmp);runner.HERE=folder
            try:runner.main()
            except AssertionError as error:assert 'calibration' in str(error)
            else:raise AssertionError('Missing calibration did not block launch')
            assert calls==[]
            schedule=[dict(instance='g_synthetic'+str(i),seed=seed)
                      for i in range(9) for seed in range(81,91)]
            random.Random(20261020).shuffle(schedule)
            # Silence per-pair JSON in the fixture only.
            runner.print=lambda *args,**kwargs:None
            runner.execute_pairs(schedule,folder,'synthetic-protocol')
            assert len(calls)==268
            cursor=0
            for pair in schedule:
                name,seed=pair['instance'],pair['seed']
                assert calls[cursor]==('reference',name,seed);cursor+=1
                if budgets[(name,seed)] is None:continue
                order=['hybrid10','hybrid60'];random.Random(name+str(seed)+'paired70').shuffle(order)
                assert calls[cursor:cursor+2]==[(mode,name,seed) for mode in order];cursor+=2
            completion=json.loads((folder/'completion.json').read_text())
            assert completion['reference_records']==90 and completion['candidate_records']==178
            calls.clear();runner.execute_pairs(schedule,folder,'synthetic-protocol')
            assert calls==[], 'Resume duplicated completed pairs'
            target=next(folder.glob('g_*.pair.json'));record=json.loads(target.read_text())
            record['protocol_sha256']='different';write(target,record)
            try:runner.execute_pairs(schedule,folder,'synthetic-protocol')
            except AssertionError:pass
            else:raise AssertionError('Changed protocol was accepted')
            assert calls==[]
    finally:
        runner.reference,runner.candidate,runner.load,runner.HERE=saved
        if hasattr(runner,'print'):del runner.print
    report=dict(status='PASS',synthetic_pairs=90,synthetic_optimizer_calls=268,
        checks=['missing calibration blocks all calls','reference before candidates',
            'both candidates get exact reference T','prespecified candidate order',
            'censored reference skips both candidates','resume makes no optimizer calls',
            'changed protocol rejected'],source_hashes={n:digest(Path(__file__).with_name(n)) for n in
            ('paired70_train.py','paired70_reference_audit.py','paired70_candidate_audit.py')},
        scope='mocked optimizers, synthetic data and real scheduler/audits; not a PCCTP experiment')
    write(Path(__file__).with_suffix('.json'),report);print(json.dumps(report))
