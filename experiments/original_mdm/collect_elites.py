"""Export the occupied elite at each mining call of the corrected MDM.

The search, the construction and the miner stay as in the bounded baseline.
Development geometries only. One seed. This is a bank for a later selector
diagnostic, not a performance claim.
"""
import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from experiments.original_mdm.run import run


def cases(limit):
    chosen=[]; seen=set()
    for path in sorted((ROOT/'experiments/ml_lns/decoder_initialized_states').glob('g_*.json')):
        row=json.loads(path.read_text())
        if row['role']!='development' or row['geometry'] in seen: continue
        seen.add(row['geometry']); chosen.append(row)
        if len(chosen)==limit: break
    return chosen


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--limit',type=int,default=9)
    parser.add_argument('--seed',type=int,default=42)
    parser.add_argument('--timeout',type=float,default=120)
    args=parser.parse_args()
    out=HERE/'elite_banks'; out.mkdir(exist_ok=True)
    rows=[]
    for case in cases(args.limit):
        target=out/f"{case['instance']}.seed{args.seed}.json"
        if target.exists():
            row=json.loads(target.read_text())
        else:
            row=run(ROOT/'data'/case['instance'],args.seed,target,timeout=args.timeout,
                     initialize_decoder_index=True,bound_elite_size=True,export_elite=True)
        sidecar=Path(str(target.with_suffix('.raw'))+'.elites.jsonl')
        banks=sidecar.read_text().splitlines() if sidecar.exists() else []
        summary=dict(instance=case['instance'],geometry=case['geometry'],seed=args.seed,
                     role='development',mining_calls=row.get('mining_calls'),completed=row.get('completed'),
                     final_cost=row.get('cost'),elite_lines=len(banks),sidecar=str(sidecar))
        rows.append(summary)
        print(json.dumps(summary),flush=True)
    (out/'collection.json').write_text(json.dumps(rows,indent=2)+'\n')


if __name__=='__main__':
    main()
