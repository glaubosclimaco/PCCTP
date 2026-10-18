"""Freeze geometry split and generate independently validated legacy strong states.

Uses the supplied original MDM sources with observational export hooks.
Keeps development geometries separate; split uses only geometry names, not labels.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import subprocess
import time
from pilot import read_instance, validate

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).with_name('original_states')
from experiments.original_mdm.run import run as original_run


def freeze_split():
    target=OUT/'protocol.json'
    if target.exists():
        return json.loads(target.read_text())
    OUT.mkdir(exist_ok=True)
    prior=json.loads(Path(__file__).with_name('geometry_pilot.json').read_text())
    development=sorted({row['geometry'] for row in prior['instances']})
    groups={}
    for path in sorted((ROOT/'data').glob('g_*')):
        d=read_instance(path); groups.setdefault(d['geometry'],[]).append(path.name)
    unseen=sorted(set(groups)-set(development),key=lambda g:hashlib.sha256(('pcctp-study-v1:'+g).encode()).hexdigest())
    roles={g:'development' for g in development}
    roles.update({g:'train' for g in unseen[:20]})
    roles.update({g:'validation' for g in unseen[20:27]})
    roles.update({g:'test' for g in unseen[27:]})
    plan=dict(scope='Exploratory development; new test labels must not guide model changes',
              baseline='Supplied original MDM + FPMax, isolated observational instrumentation',
              iterations=70,wall_timeout_seconds=180,support_absolute=2,stability_percent=5,seeds=[42,43,44],
              selection='Lexicographically first variant per geometry; split by SHA256 of name',
              roles=roles,instances={g:groups[g][0] for g in sorted(groups)},
              primary='Cross-geometry held-out benefit of joint feasibility descriptors beyond generic route descriptors',
              budget=dict(k=8,candidates=24,repair_seconds=.5),
              success='Positive held-out practical benefit with uncertainty and strong controls; no outcome-based stopping rule')
    target.write_text(json.dumps(plan,indent=2)+'\n')
    return plan


def generate(job):
    geometry, name, seed, plan=job
    dest=OUT/f'{name}.seed{seed}.json'
    if not dest.exists():
        original_run(ROOT/'data'/name,seed,dest,support=plan['support_absolute'],
                     stability=plan['stability_percent'],timeout=plan['wall_timeout_seconds'])
    row=json.loads(dest.read_text()); d=read_instance(ROOT/'data'/name)
    mapping={int(original):i for i,original in enumerate(d['visits'])}
    measured=validate(d,[mapping[v] for v in row['route']])
    if abs(measured-row['cost'])>1e-6:
        raise ValueError('Independent objective mismatch')
    for state in row['snapshots']:
        if abs(validate(d,[mapping[v] for v in state['route']])-state['cost'])>1e-6:
            raise ValueError('Invalid snapshot')
    row.update(instance=name,geometry=geometry,role=plan['roles'][geometry],independently_validated=True)
    dest.write_text(json.dumps(row,indent=2)+'\n')
    return dict(instance=name,geometry=geometry,seed=seed,cost=row['cost'],
                iterations=row['iterations'],mining_calls=row['mining_calls'],seconds=row['seconds'])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--roles',nargs='+',default=['development','train'])
    parser.add_argument('--seeds',nargs='+',type=int,default=[42])
    parser.add_argument('--workers',type=int,default=2)
    parser.add_argument('--limit',type=int)
    args=parser.parse_args(); plan=freeze_split()
    geoms=[g for g in plan['instances'] if plan['roles'][g] in args.roles]
    if args.limit: geoms=geoms[:args.limit]
    jobs=[(g,plan['instances'][g],seed,plan) for g in geoms for seed in args.seeds]
    began=time.perf_counter(); failures=[]
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures={pool.submit(generate,job):job for job in jobs}
        for future in as_completed(futures):
            try: print(json.dumps(future.result()),flush=True)
            except Exception as error:
                failure=dict(job=futures[future][:3],error=str(error)); failures.append(failure)
                print(json.dumps(failure),flush=True)
    (OUT/'last_generation.json').write_text(json.dumps(dict(jobs=len(jobs),failures=failures,
                                                         elapsed_seconds=time.perf_counter()-began),indent=2)+'\n')
    if failures: raise SystemExit(1)


if __name__=='__main__': main()
