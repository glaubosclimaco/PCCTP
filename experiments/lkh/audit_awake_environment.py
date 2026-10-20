"""Independent environment-record and result-hash verification for awake300."""
import hashlib,json,math
from pathlib import Path
HERE=Path(__file__).resolve().parent

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def audit(folder):
    folder=Path(folder);protocol=json.loads((folder/'protocol.json').read_text())
    expected={(c['instance'],s) for c in protocol['cases'] for s in (70,71,72)}
    arms=set(protocol['arms']);seen=set();rows=[]
    for path in sorted(folder.glob('environment.g_*.json')):
        r=json.loads(path.read_text());key=r['instance'],r['seed']
        assert key in expected and key not in seen;seen.add(key)
        assert r['status']=='PASS' and set(r['order'])==arms and len(r['order'])==3
        assert set(r['result_hashes'])==arms and not r['sleep_events']
        assert math.isfinite(r['max_monitor_gap_seconds']) and 0<=r['max_monitor_gap_seconds']<=5
        assert math.isfinite(r['max_clock_difference_seconds']) and 0<=r['max_clock_difference_seconds']<=3
        assert math.isfinite(r['started_unix']) and math.isfinite(r['finished_unix'])
        assert r['started_unix']<=r['finished_unix']
        for mode,hashed in r['result_hashes'].items():
            result=folder/f"{r['instance']}.seed{r['seed']}.{mode}.json"
            assert digest(result)==hashed
        rows.append(dict(instance=r['instance'],seed=r['seed'],environment_sha256=digest(path),
            max_monitor_gap_seconds=r['max_monitor_gap_seconds'],
            max_clock_difference_seconds=r['max_clock_difference_seconds']))
    complete=(folder/'completion.json').exists()
    if complete:
        assert seen==expected and len(rows)==27
        summary=json.loads((folder/'environment_audit.json').read_text())
        assert summary['status']=='PASS' and summary['blocks']==27
        assert summary['records']=={p.name:digest(p) for p in folder.glob('environment.g_*.json')}
    return dict(status='PASS',batch_complete=complete,validated_blocks=len(rows),expected_blocks=27,
        eligible_complete_environment=complete and seen==expected,rows=rows,
        scope='independent recorded clock/pmset/result-hash checks; not proof of equal CPU or thermal conditions')

if __name__=='__main__':
    result=audit(HERE/'ruin_total300_awake_train')
    (HERE/'ruin_total300_awake_environment_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='rows'}))
