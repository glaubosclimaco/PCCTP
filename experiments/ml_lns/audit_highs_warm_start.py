"""Independent saved-route audit of the native HiGHS paired screening run."""
import hashlib
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from experiments.ml_lns.pilot import read_instance,validate


def audit():
    here=Path(__file__).resolve().parent
    source=here/'highs_warm_start_results.json'
    results=json.loads(source.read_text())
    assert len(results['rows'])==11,'Batch not complete; do not audit checkpoint as final'
    seen=set();rows=[]
    for row in results['rows']:
        assert row['geometry'] not in seen;seen.add(row['geometry'])
        path=ROOT/'data'/row['instance'];d=read_instance(path)
        mapping={int(v):i for i,v in enumerate(d['visits'])}
        profile_path=ROOT/'experiments/original_mdm/move_profiles_monotonic_elite_bounded'/f"{row['instance']}.seed42.profiled_baseline.json"
        profile=json.loads(profile_path.read_text());assert profile['completed'] and profile['iterations']==70
        incumbent=validate(d,[mapping[v] for v in profile['route']])
        assert row['warm']['version']==row['cold']['version']=='1.13.1'
        for mode in ('cold','warm'):
            r=row[mode];cost=validate(d,[mapping[v] for v in r['route']])
            assert abs(cost-r['cost'])<1e-6 and abs(incumbent-r['initial_cost'])<1e-6
            assert cost<=incumbent+1e-6
            assert r['improved']==(cost<incumbent-1e-7)
            assert r['dual_bound']<=cost+1e-4
            previous_time=0.;previous_cost=incumbent
            for event in r['events']:
                assert previous_time<=event['seconds']<=10.
                assert event['cost']<previous_cost-1e-7
                previous_time=event['seconds'];previous_cost=event['cost']
            if r['improved']:assert abs(previous_cost-cost)<1e-6
            if r['proven_optimal']:
                assert r['status']=='Optimal' and r['seconds']<=10.
                assert abs(r['dual_bound']-cost)<1e-4
        assert row['warm']['start_status']=='HighsStatus.kOk'
        rows.append(dict(instance=row['instance'],input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                         warm_improved=row['warm']['improved'],cold_improved=row['cold']['improved'],
                         warm_proven=row['warm']['proven_optimal'],cold_proven=row['cold']['proven_optimal']))
    count=sum(r['warm_improved'] for r in rows)
    assert count==results['warm_improvements'] and results['gate_passed']==(count>=2)
    report=dict(status='PASS',geometries=11,warm_improvements=count,
                cold_improvements=sum(r['cold_improved'] for r in rows),
                warm_proofs=sum(r['warm_proven'] for r in rows),cold_proofs=sum(r['cold_proven'] for r in rows),
                source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),rows=rows)
    (here/'highs_warm_start_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='rows'}))


if __name__=='__main__':audit()
