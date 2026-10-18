"""Observe unchanged MDM operators before attempting learned scheduling.

First 2000 calls per run only. This is an explicitly censored observational
sample, not counterfactual rewards for unexecuted operators. Profiling adds I/O;
profile run wall times cannot be used as baseline speed measurements.
"""
import hashlib
import argparse
import json
from pathlib import Path
import sys
from collections import defaultdict

ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from experiments.original_mdm.run import run
from experiments.ml_lns.pilot import read_instance,validate

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--train-limit',type=int,default=4)
    parser.add_argument('--dev-limit',type=int,default=2)
    parser.add_argument('--include-incomplete',action='store_true')
    parser.add_argument('--output-dir',type=Path,default=HERE/'move_profiles')
    parser.add_argument('--bound-elite-size',action='store_true')
    args=parser.parse_args()
    bank=ROOT/'experiments/ml_lns/decoder_initialized_states';out=args.output_dir.resolve();out.mkdir(exist_ok=True)
    lock=out/'.collector.lock'
    lock_fd=__import__('os').open(lock,__import__('os').O_CREAT|__import__('os').O_EXCL|__import__('os').O_WRONLY)
    __import__('os').write(lock_fd,str(__import__('os').getpid()).encode());__import__('os').close(lock_fd)
    __import__('atexit').register(lambda:lock.unlink(missing_ok=True))
    cases=[];seen=set();counts=defaultdict(int)
    for path in sorted(bank.glob('g_*.json')):
        b=json.loads(path.read_text());role=b['role']
        if role not in ('development','train') or b['geometry'] in seen:continue
        if counts[role]>=(args.dev_limit if role=='development' else args.train_limit):continue
        if not args.include_incomplete and not b.get('completed',True):continue
        counts[role]+=1;seen.add(b['geometry']);cases.append((path,b))
    build=HERE/('build_move_profile_elite_bounded' if args.bound_elite_size else 'build_move_profile')
    protocol=dict(scope=__doc__,elite_count_bounded=args.bound_elite_size,cases=[dict(instance=b['instance'],geometry=b['geometry'],role=b['role'],
                  baseline_sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p,b in cases],
                  build_manifest=json.loads((build/'manifest.json').read_text()),timeout=90,
                  code_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in
                               ['profile_moves.py','run.py','research_hook.cpp','research_hook.hpp','build.py']},
                  binary_sha256=hashlib.sha256((build/'original_mdm').read_bytes()).hexdigest())
    fingerprint=hashlib.sha256(json.dumps(protocol,sort_keys=True).encode()).hexdigest()[:12]
    (out/f'protocol.{fingerprint}.json').write_text(json.dumps(protocol,indent=2)+'\n')
    failures=[]
    for path,b in cases:
        target=out/f"{b['instance']}.seed{b['seed']}.{fingerprint}.json"
        if target.exists():continue
        # Reuse fully validated earlier profiles from the identical native build.
        existing=[]
        for candidate in out.glob(f"{b['instance']}.seed{b['seed']}.*.json"):
            if candidate.name.endswith('.profiled_baseline.json'):continue
            prior=json.loads(candidate.read_text())
            if 'records' not in prior:continue
            previous=json.loads((out/f"protocol.{prior['protocol']}.json").read_text())
            if (previous['build_manifest']==protocol['build_manifest'] and
                prior['source_sha256']==hashlib.sha256(path.read_bytes()).hexdigest()):existing.append(candidate)
        if existing:continue
        baseline_out=out/f"{b['instance']}.seed{b['seed']}.profiled_baseline.json"
        try:
            result=run(ROOT/'data'/b['instance'],b['seed'],baseline_out,timeout=90,initialize_decoder_index=True,profile_moves=True,bound_elite_size=args.bound_elite_size)
        except RuntimeError as error:
            failures.append(dict(instance=b['instance'],error=str(error)))
            (out/f'failures.{fingerprint}.json').write_text(json.dumps(failures,indent=2)+'\n')
            print(json.dumps(failures[-1]),flush=True)
            continue
        raw=baseline_out.with_suffix('.raw');profile=Path(str(raw)+'.moves.jsonl')
        records=[json.loads(line) for line in profile.read_text().splitlines()]
        d=read_instance(ROOT/'data'/b['instance']);mapping={int(v):i for i,v in enumerate(d['visits'])}
        for i,r in enumerate(records):
            assert r['call']==i and r['operator_seconds']>=0
            for prefix in ('before','after'):
                cost=validate(d,[mapping[v] for v in r[prefix+'_route']])
                assert abs(cost-r[prefix+'_cost'])<1e-6
        common={s['iteration']:s for s in b['snapshots']}
        differences=[dict(iteration=s['iteration'],before=common[s['iteration']]['cost'],after=s['cost'])
                     for s in result['snapshots'] if s['iteration'] in common and s['cost']!=common[s['iteration']]['cost']]
        if not args.bound_elite_size:assert not differences,'Profiling changed a recorded cost'
        grouped=defaultdict(list)
        for r in records:grouped[r['move']].append(r)
        summary={str(m):dict(calls=len(rr),improving=sum(r['after_cost']<r['before_cost'] for r in rr),
                           seconds=sum(r['operator_seconds'] for r in rr),
                           total_gain=sum(max(0.,r['before_cost']-r['after_cost']) for r in rr)) for m,rr in grouped.items()}
        row=dict(instance=b['instance'],geometry=b['geometry'],role=b['role'],seed=b['seed'],protocol=fingerprint,
                 clock_source=protocol['build_manifest'].get('clock_source','legacy_chrono_unverified'),
                 baseline_method=result['method'],elite_count_bounded=args.bound_elite_size,baseline_checkpoint_differences=differences,
                 source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),recorded_calls=len(records),
                 sample_capped=len(records)==2000,profiled_run_completed=result['completed'],
                 final_cost=result['cost'],baseline_cost=b['cost'],summary=summary,records=records)
        temp=target.with_suffix('.partial');temp.write_text(json.dumps(row,indent=2)+'\n');temp.replace(target)
        print(json.dumps({k:v for k,v in row.items() if k not in ('summary','records')}),flush=True)
    (out/f'completion.{fingerprint}.json').write_text(json.dumps(dict(planned=len(cases),failures=failures),indent=2)+'\n')

if __name__=='__main__':main()
