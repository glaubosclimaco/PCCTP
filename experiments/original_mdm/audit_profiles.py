"""Independent snapshot audit of corrected, monotonic TRAIN operator profiles."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from experiments.ml_lns.pilot import read_instance,validate

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profile-dir',type=Path,default=HERE/'move_profiles_monotonic_elite_bounded')
    args=parser.parse_args();directory=args.profile_dir.resolve();rows=[];cache={}
    for path in sorted(directory.glob('g_*.json')):
        r=json.loads(path.read_text())
        if 'records' not in r:continue
        assert r['role']=='train' and r['elite_count_bounded'] is True
        assert r['clock_source']=='clock_gettime(CLOCK_MONOTONIC)'
        assert r['baseline_method']=='original_MDM_legacy_decoder_index_initialized_elite_count_bounded'
        source=ROOT/'experiments/ml_lns/decoder_initialized_states'/f"{r['instance']}.seed{r['seed']}.json"
        assert hashlib.sha256(source.read_bytes()).hexdigest()==r['source_sha256']
        baseline=json.loads(source.read_text());assert baseline['role']=='train' and baseline['geometry']==r['geometry']
        protocol=json.loads((directory/f"protocol.{r['protocol']}.json").read_text())
        assert protocol['elite_count_bounded'] and protocol['build_manifest']['clock_source']==r['clock_source']
        assert protocol['build_manifest']['instrumented_sha256']
        if r['instance'] not in cache:
            d=read_instance(ROOT/'data'/r['instance']);cache[r['instance']]=(d,{int(v):i for i,v in enumerate(d['visits'])},{})
        d,mapping,known=cache[r['instance']]
        def cost(route):
            key=tuple(route)
            if key not in known:known[key]=validate(d,[mapping[v] for v in route])
            return known[key]
        improving=0
        for i,record in enumerate(r['records']):
            assert record['call']==i and 0<=record['move']<15 and record['operator_seconds']>=0
            for prefix in ['before','after']:assert abs(cost(record[prefix+'_route'])-record[prefix+'_cost'])<1e-6
            improving+=record['after_cost']<record['before_cost']
        assert r['recorded_calls']==len(r['records'])<=2000
        assert r['sample_capped']==(len(r['records'])==2000)
        native=json.loads((directory/f"{r['instance']}.seed{r['seed']}.profiled_baseline.json").read_text())
        assert native['method']==r['baseline_method'] and native['elite_count_bounded']
        for state in [native,*native['snapshots']]:assert abs(cost(state['route'])-state['cost'])<1e-6
        assert native['cost']==r['final_cost'] and native['completed']==r['profiled_run_completed']
        rows.append(dict(file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                         instance=r['instance'],geometry=r['geometry'],calls=len(r['records']),improving=improving,
                         completed=r['profiled_run_completed'],sample_capped=r['sample_capped'],
                         final_cost=r['final_cost'],baseline_checkpoint_differences=r['baseline_checkpoint_differences']))
    plans=list(directory.glob('protocol.*.json'));assert len(plans)==1
    plan=json.loads(plans[0].read_text());completion=directory/f"completion.{plans[0].stem.split('.')[-1]}.json"
    complete=completion.exists();failure_count=len(json.loads(completion.read_text())['failures']) if complete else None
    if complete:assert len(rows)+failure_count==len(plan['cases'])
    result=dict(status='PASS',partial_snapshot=not complete,planned=len(plan['cases']),files=len(rows),
                calls=sum(r['calls'] for r in rows),improving=sum(r['improving'] for r in rows),
                completed_runs=sum(r['completed'] for r in rows),checkpoint_only=sum(not r['completed'] for r in rows),
                generation_failures=failure_count,rows=rows,
                code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (HERE/'monotonic_profile_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))

if __name__=='__main__':main()
