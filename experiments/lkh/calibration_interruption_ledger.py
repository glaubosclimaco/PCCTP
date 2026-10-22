"""Passive inventory of excluded calibration attempts; never promotes outcomes."""
import hashlib,json
from pathlib import Path
HERE=Path(__file__).resolve().parent

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def build():
    rows=[]
    for folder in sorted((HERE/'published70_calibration_interrupted').iterdir()):
        review=folder/'review.json';envpath=folder/'environment_failure.json'
        if not review.exists() or not envpath.exists():continue
        record=json.loads(review.read_text());env=json.loads(envpath.read_text())
        assert record['status']=='EXCLUDED_ENVIRONMENT_FAILURE'
        assert all(digest(folder/name)==sha for name,sha in record['files'].items())
        assert env['status']=='ENVIRONMENT_FAILURE'
        native=list(folder.glob('*.native.json'));assert len(native)==1
        run=json.loads(native[0].read_text());assert run['instance']==env['instance']
        rows.append(dict(archive=folder.name,instance=env['instance'],seed=run['seed'],
            exclusion='SLEEP_RECORDED' if env['sleep_events'] else 'MONITOR_GAP_WITHOUT_SLEEP_RECORD',
            sleep_events=len(env['sleep_events']),max_monitor_gap_seconds=env['max_monitor_gap_seconds'],
            max_clock_difference_seconds=env['max_clock_difference_seconds'],
            excluded_process_wall_seconds=run['process_wall_seconds'],
            excluded_child_cpu_seconds=run['child_user_cpu_seconds']+run['child_system_cpu_seconds'],
            review_sha256=digest(review),environment_sha256=digest(envpath)))
    return dict(status='PRESERVATION_INVENTORY',excluded_attempts=len(rows),
        sleep_recorded_attempts=sum(r['sleep_events']>0 for r in rows),
        no_sleep_recorded_attempts=sum(r['sleep_events']==0 for r in rows),
        excluded_process_wall_seconds=sum(r['excluded_process_wall_seconds'] for r in rows),
        rows=rows,scope='Excluded operational attempts only, not algorithm performance, runtime calibration, or primary cost results. A monitor gap without a Sleep record does not identify its cause. No threshold changes or retrospective acceptance.')

if __name__=='__main__':
    result=build();result['source_sha256']=digest(Path(__file__))
    (HERE/'calibration_interruption_ledger.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='rows'}))
