"""Recompute feasibility/objectives and split provenance for every collected route."""
import hashlib
import argparse
import json
from pathlib import Path
import numpy as np
from pilot import read_instance,validate

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--state-dir',type=Path,default=HERE/'original_states')
    parser.add_argument('--oracle-dir',type=Path)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    args.state_dir=args.state_dir.resolve()
    original=args.state_dir.resolve()==(HERE/'original_states').resolve()
    if args.oracle_dir is None: args.oracle_dir=HERE/('strong_oracle' if original else 'strong_oracle_'+args.state_dir.name)
    args.oracle_dir=args.oracle_dir.resolve()
    if args.output is None: args.output=HERE/('dataset_audit.json' if original else 'dataset_audit_'+args.state_dir.name+'.json')
    plan=json.loads((args.state_dir/'protocol.json').read_text())
    expected_method=('original_MDM_legacy_decoder_index_initialized' if plan.get('initialize_decoder_index',False)
                     else 'original_MDM_legacy_instrumented')
    cache={}; failures=[]; originals=0; snapshots=0; queries=0; cases=0
    manifest={}; incomplete=0; completed=0
    def load(name):
        if name not in cache:
            d=read_instance(ROOT/'data'/name)
            cache[name]=(d,{int(v):i for i,v in enumerate(d['visits'])})
        return cache[name]
    for path in sorted(args.state_dir.glob('g_*.json')):
        try:
            row=json.loads(path.read_text()); d,mapping=load(row['instance'])
            assert row['role']==plan['roles'][row['geometry']]
            assert row['method']==expected_method,'Mixed baseline methods'
            assert row['role'] in ('development','train'),'Reserved geometry labels generated'
            for state in [row,*row['snapshots']]:
                assert abs(validate(d,[mapping[v] for v in state['route']])-state['cost'])<1e-6
            if 'input_sha256' in row:
                assert row['input_sha256']==hashlib.sha256((ROOT/'data'/row['instance']).read_bytes()).hexdigest()
            originals+=1; snapshots+=len(row['snapshots'])
            if row.get('completed',True): completed+=1
            else: incomplete+=1
            manifest[str(path.relative_to(HERE))]=hashlib.sha256(path.read_bytes()).hexdigest()
        except Exception as error: failures.append(dict(file=path.name,error=str(error)))
    for path in sorted(args.oracle_dir.glob('*.json')):
        try:
            row=json.loads(path.read_text()); d,mapping=load(row['instance'])
            assert row['role']==plan['roles'][row['geometry']]
            assert row['baseline_method']==expected_method,'Mixed repair label provenance'
            assert row['role'] in ('development','train')
            assert abs(validate(d,[mapping[v] for v in row['initial_route']])-row['initial_cost'])<1e-6
            for r in row['records']:
                assert abs(validate(d,[mapping[v] for v in r['route']])-r['cost'])<1e-6
                assert r['cost']<=row['initial_cost']+1e-6
                expected=(row['initial_cost']-r['cost'])/max(1.,row['initial_cost'])
                assert abs(expected-r['normalized_gain'])<1e-12
                assert len(r['features'])==len(row['feature_names'])
                assert np.isfinite(r['features']).all()
                assert r['status'] in (0,1)
                assert r['proven_optimal']==(r['status']==0)
                queries+=1
            cases+=1
            manifest[str(path.relative_to(HERE))]=hashlib.sha256(path.read_bytes()).hexdigest()
        except Exception as error: failures.append(dict(file=path.name,error=str(error)))
    result=dict(status='PASS' if not failures else 'FAIL',baseline_method=expected_method,baseline_files=originals,
                completed_baselines=completed,checkpoint_only_baselines=incomplete,
                baseline_snapshots=snapshots,label_states=cases,label_queries=queries,
                geometries=len({d['geometry'] for d,_ in cache.values()}),
                reserved_geometry_labels_found=False if not failures else 'Inspect failures',
                failures=failures,source_hashes=manifest,
                scope='Snapshot of completed files only; generation failures/timeouts are separately recorded in manifests')
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))
    if failures: raise SystemExit(1)

if __name__=='__main__': main()
