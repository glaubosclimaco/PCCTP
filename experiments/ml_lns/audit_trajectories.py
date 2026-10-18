"""Independently check final routes, best-so-far traces and training separation."""
import collections
import hashlib
import json
from pathlib import Path
from pilot import read_instance,validate

HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[1]

def main():
    folder=HERE/'hybrid_development'; protocols={}; failures=[]; grouped=collections.defaultdict(list)
    for path in sorted(folder.glob('protocol.*.json')):
        protocol=json.loads(path.read_text()); fingerprint=path.name.split('.')[1]
        assert hashlib.sha256(json.dumps(protocol,sort_keys=True).encode()).hexdigest()[:12]==fingerprint
        protocols[fingerprint]=protocol
    for path in sorted(folder.glob('g_*.json')):
        try:
            row=json.loads(path.read_text()); protocol=protocols[row['model_fingerprint']]
            assert row['geometry'] not in protocol['training_geometries']
            assert row['method'] in protocol['methods']
            d=read_instance(ROOT/'data'/row['instance']); mapping={int(v):i for i,v in enumerate(d['visits'])}
            assert abs(validate(d,[mapping[v] for v in row['route']])-row['cost'])<1e-6
            baseline=json.loads((HERE/'decoder_initialized_states'/f"{row['instance']}.seed{row['seed']}.json").read_text())
            assert baseline['role']=='development'
            initial=next(s for s in baseline['snapshots'] if s['iteration']==row['initial_phase'])
            assert abs(row['initial_cost']-initial['cost'])<1e-6
            if 'evaluation_cases' in protocol:
                case=next(c for c in protocol['evaluation_cases'] if c['instance']==row['instance'] and c['seed']==row['seed'])
                assert case['geometry']==row['geometry']
            previous=row['initial_cost']; clock=0.
            for point in row['trace']:
                assert point['cost']<=previous+1e-7
                assert clock<=point['elapsed']<=row['total_seconds']+1e-6
                previous=point['cost']; clock=point['elapsed']
            assert abs(previous-row['cost'])<1e-6
            assert abs(row['gain']-(row['initial_cost']-row['cost'])/max(1.,row['initial_cost']))<1e-12
            grouped[row['model_fingerprint']].append(row)
        except Exception as error: failures.append(dict(file=path.name,error=str(error)))
    summaries={}
    for fingerprint,rows in grouped.items():
        methods={}
        for method in sorted({r['method'] for r in rows}):
            subset=[r for r in rows if r['method']==method]
            geometries=sorted({r['geometry'] for r in subset})
            means=[sum(r['gain'] for r in subset if r['geometry']==g)/sum(r['geometry']==g for r in subset) for g in geometries]
            methods[method]=dict(instances=len(subset),geometries=len(geometries),mean_geometry_gain=sum(means)/len(means),
                                 matches_final_MDM=sum(r['cost']<=r['baseline_final_cost']+1e-6 for r in subset),
                                 beats_recorded_MDM=sum(r['cost']<r['recorded_MDM_cost_before_end']-1e-6 for r in subset))
        summaries[fingerprint]=dict(acceptance=protocols[fingerprint].get('acceptance','greedy'),methods=methods)
    result=dict(status='PASS' if not failures else 'FAIL',trajectories=sum(map(len,grouped.values())),
                failures=failures,scope='Development audit; protocol groups may be incomplete; no speed or significance claim',summaries=summaries)
    (HERE/'trajectory_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='summaries'},indent=2))
    if failures: raise SystemExit(1)

if __name__=='__main__': main()
