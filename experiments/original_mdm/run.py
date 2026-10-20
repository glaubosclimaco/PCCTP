"""Run supplied MDM+FPMax with private temporary files and validate every export."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from experiments.ml_lns.pilot import read_instance, validate


def run(instance,seed,output,support=2,stability=5,timeout=180,initialize_decoder_index=False,profile_moves=False,operator_policy=None,bound_elite_size=False,export_elite=False):
    instance=Path(instance).resolve(); output=Path(output).resolve()
    if profile_moves and not initialize_decoder_index:raise ValueError('Profiling requires explicit decoder correction')
    if operator_policy is not None and (operator_policy not in ['original','static','size','generic','joint','adaptive'] or profile_moves or not initialize_decoder_index):raise ValueError('Invalid scheduled run')
    if export_elite and (operator_policy or profile_moves or not initialize_decoder_index or not bound_elite_size):
        raise ValueError('Elite export requires the corrected bounded build and no other instrumentation')
    build=Path(__file__).with_name('build_elite_export' if export_elite else 'build_scheduled' if operator_policy else 'build_move_profile' if profile_moves else 'build_decoder_initialized' if initialize_decoder_index else 'build').resolve()
    if bound_elite_size and not export_elite:build=build.with_name(build.name+'_elite_bounded')
    working=Path(__file__).with_name('runs')/(instance.name+f'.seed{seed}'+('.elite_export' if export_elite else '.scheduled.'+operator_policy if operator_policy else '.profile' if profile_moves else '.decoder_initialized' if initialize_decoder_index else ''))
    if bound_elite_size and not export_elite:working=working.with_name(working.name+'.elite_bounded')
    fp=working/'FPMax'; fp.mkdir(parents=True,exist_ok=True)
    for name in ['fim_maximal','fpmax_hnmp']:
        target=fp/name
        if not target.exists(): target.symlink_to(build/'FPMax'/name)
    output.parent.mkdir(parents=True,exist_ok=True)
    log=output.with_suffix('.log')
    raw=output.with_suffix('.raw')
    raw.unlink(missing_ok=True)
    if export_elite:Path(str(raw)+'.elites.jsonl').unlink(missing_ok=True)
    if profile_moves:Path(str(raw)+'.moves.jsonl').unlink(missing_ok=True)
    checkpoint=Path(str(raw)+'.checkpoint')
    checkpoint.unlink(missing_ok=True)
    timed_out=False
    process_began=time.perf_counter()
    environment=os.environ.copy()
    if operator_policy:environment['PCCTP_OPERATOR_POLICY']=operator_policy
    with log.open('w') as f:
        try:
            process=subprocess.run([str(build/'original_mdm'),str(instance),str(seed),str(support),str(stability),str(raw)],
                                   cwd=working,stdout=f,stderr=f,timeout=timeout,env=environment)
        except subprocess.TimeoutExpired:
            timed_out=True
    process_wall_seconds=time.perf_counter()-process_began
    if timed_out:
        if not checkpoint.exists(): raise RuntimeError(f'MDM timed out before its first checkpoint: {log}')
        raw.write_text(checkpoint.read_text())
    elif process.returncode or not raw.exists():
        raise RuntimeError(f'MDM failed, exit {process.returncode}: {log}')
    row=json.loads(raw.read_text()); d=read_instance(instance)
    mapping={int(original):i for i,original in enumerate(d['visits'])}
    for state in [row,*row['snapshots']]:
        measured=validate(d,[mapping[v] for v in state['route']])
        if abs(measured-state['cost'])>1e-6: raise ValueError('MDM reported cost mismatch')
    row.update(instance=instance.name,geometry=d['geometry'],independently_validated=True,
               move_profiling=profile_moves,
               operator_policy=operator_policy,
               elite_count_bounded=bound_elite_size,
               elite_export=export_elite,
               timed_out=timed_out,wall_timeout_seconds=timeout,
               process_wall_seconds=process_wall_seconds,
               input_sha256=__import__('hashlib').sha256(instance.read_bytes()).hexdigest())
    if operator_policy:row['method']='scheduled_MDM_'+operator_policy
    validated=output.with_suffix('.validated.tmp')
    validated.write_text(json.dumps(row,indent=2)+'\n')
    validated.replace(output)
    return row


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('instance',type=Path)
    parser.add_argument('--seed',type=int,default=42)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    result=run(args.instance,args.seed,args.output)
    print(json.dumps({k:v for k,v in result.items() if k not in ('route','snapshots')}))
