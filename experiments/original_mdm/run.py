"""Run supplied MDM+FPMax with private temporary files and validate every export."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from experiments.ml_lns.pilot import read_instance, validate


def run(instance,seed,output,support=2,stability=5,timeout=180):
    instance=Path(instance).resolve(); output=Path(output).resolve()
    build=Path(__file__).with_name('build').resolve()
    working=Path(__file__).with_name('runs')/(instance.name+f'.seed{seed}')
    fp=working/'FPMax'; fp.mkdir(parents=True,exist_ok=True)
    for name in ['fim_maximal','fpmax_hnmp']:
        target=fp/name
        if not target.exists(): target.symlink_to(build/'FPMax'/name)
    output.parent.mkdir(parents=True,exist_ok=True)
    log=output.with_suffix('.log')
    checkpoint=Path(str(output)+'.checkpoint')
    checkpoint.unlink(missing_ok=True)
    timed_out=False
    with log.open('w') as f:
        try:
            process=subprocess.run([str(build/'original_mdm'),str(instance),str(seed),str(support),str(stability),str(output)],
                                   cwd=working,stdout=f,stderr=f,timeout=timeout)
        except subprocess.TimeoutExpired:
            timed_out=True
    if timed_out:
        if not checkpoint.exists(): raise RuntimeError(f'MDM timed out before its first checkpoint: {log}')
        output.write_text(checkpoint.read_text())
    elif process.returncode or not output.exists():
        raise RuntimeError(f'MDM failed, exit {process.returncode}: {log}')
    row=json.loads(output.read_text()); d=read_instance(instance)
    mapping={int(original):i for i,original in enumerate(d['visits'])}
    for state in [row,*row['snapshots']]:
        measured=validate(d,[mapping[v] for v in state['route']])
        if abs(measured-state['cost'])>1e-6: raise ValueError('MDM reported cost mismatch')
    row.update(instance=instance.name,geometry=d['geometry'],independently_validated=True,
               timed_out=timed_out,wall_timeout_seconds=timeout,
               input_sha256=__import__('hashlib').sha256(instance.read_bytes()).hexdigest())
    output.write_text(json.dumps(row,indent=2)+'\n')
    return row


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('instance',type=Path)
    parser.add_argument('--seed',type=int,default=42)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    result=run(args.instance,args.seed,args.output)
    print(json.dumps({k:v for k,v in result.items() if k not in ('route','snapshots')}))
