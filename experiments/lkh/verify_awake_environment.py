"""Finite synthetic environment-auditor checks, without an optimizer."""
import copy,hashlib,json,tempfile
from pathlib import Path
from audit_awake_environment import audit

def write(path,value):path.write_text(json.dumps(value)+'\n')
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

if __name__=='__main__':
    rejected=[]
    with tempfile.TemporaryDirectory(prefix='pcctp-environment-',dir='/private/tmp') as tmp:
        folder=Path(tmp);arms=['mdm','hybrid10','hybrid60']
        write(folder/'protocol.json',dict(cases=[dict(instance='g_synthetic'+str(i)) for i in range(9)],arms=arms))
        for i in range(9):
            for seed in (70,71,72):
                name='g_synthetic'+str(i);hashes={}
                for mode in arms:
                    p=folder/f'{name}.seed{seed}.{mode}.json';write(p,dict(synthetic=True));hashes[mode]=digest(p)
                r=dict(status='PASS',instance=name,seed=seed,order=arms,result_hashes=hashes,
                    started_unix=100.,finished_unix=1001.,max_monitor_gap_seconds=1.1,
                    max_clock_difference_seconds=.001,sleep_events=[])
                write(folder/f'environment.{name}.seed{seed}.json',r)
        assert audit(folder)['validated_blocks']==27
        target=folder/'environment.g_synthetic0.seed70.json';original=json.loads(target.read_text())
        for label,change in (
            ('monitor_gap',lambda r:r.update(max_monitor_gap_seconds=6.)),
            ('clock_jump',lambda r:r.update(max_clock_difference_seconds=4.)),
            ('pmset_sleep',lambda r:r.update(sleep_events=['Entering Sleep'])),
            ('nonfinite_clock',lambda r:r.update(max_monitor_gap_seconds=float('nan'))),
            ('changed_result',lambda r:r['result_hashes'].update(mdm='wrong'))):
            r=copy.deepcopy(original);change(r);write(target,r)
            try:audit(folder)
            except AssertionError:rejected.append(label)
            else:raise AssertionError('Corruption accepted: '+label)
        write(target,original)
        write(folder/'completion.json',dict(runs=81))
        write(folder/'environment_audit.json',dict(status='PASS',blocks=27,
            records={p.name:digest(p) for p in folder.glob('environment.g_*.json')}))
        assert audit(folder)['eligible_complete_environment']
    report=dict(status='PASS',synthetic_blocks=27,rejected_corruptions=rejected,
        source_sha256=digest(Path(__file__).with_name('audit_awake_environment.py')),
        scope='synthetic environment audit only; no solver, sleeping computer or benchmark execution')
    write(Path(__file__).with_suffix('.json'),report);print(json.dumps(report))
