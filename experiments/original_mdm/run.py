"""Run supplied MDM+FPMax with private temporary files and validate every export."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from experiments.ml_lns.pilot import read_instance, validate


def run(instance,seed,output,support=2,stability=5,timeout=180,initialize_decoder_index=False):
    instance=Path(instance).resolve(); output=Path(output).resolve()
    build=Path(__file__).with_name('build_decoder_initialized' if initialize_decoder_index else 'build').resolve()
    working=Path(__file__).with_name('runs')/(instance.name+f'.seed{seed}'+('.decoder_initialized' if initialize_decoder_index else ''))
    fp=working/'FPMax'; fp.mkdir(parents=True,exist_ok=True)
    for name in ['fim_maximal','fpmax_hnmp']:
        target=fp/name
        if not target.exists(): target.symlink_to(build/'FPMax'/name)
    output.parent.mkdir(parents=True,exist_ok=True)
    log=output.with_suffix('.log')
    raw=output.with_suffix('.raw')
    raw.unlink(missing_ok=True)
    checkpoint=Path(str(raw)+'.checkpoint')
    checkpoint.unlink(missing_ok=True)
    timed_out=False
    process_began=time.perf_counter()
    with log.open('w') as f:
        try:
            process=subprocess.run([str(build/'original_mdm'),str(instance),str(seed),str(support),str(stability),str(raw)],
                                   cwd=working,stdout=f,stderr=f,timeout=timeout)
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
               timed_out=timed_out,wall_timeout_seconds=timeout,
               process_wall_seconds=process_wall_seconds,
               input_sha256=__import__('hashlib').sha256(instance.read_bytes()).hexdigest())
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
