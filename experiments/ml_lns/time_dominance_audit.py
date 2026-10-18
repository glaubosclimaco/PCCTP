"""Conservative screen for gains already matched by continuing the original MDM.

Optimistically gives a repair only its own query cost: no ranking-pool feature
cost, model setup, or inference overhead. If MDM already achieved its cost at
an earlier recorded time, that local gain cannot demonstrate a speed advantage.
Absence of recorded dominance does NOT establish a repair advantage because
MDM checkpoints are sparse and original runtimes are separate executions.
"""
import argparse
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--state-dir',type=Path,default=HERE/'decoder_initialized_states')
    p.add_argument('--oracle-dir',type=Path,default=HERE/'strong_oracle_decoder_initialized_states')
    args=p.parse_args(); rows=[]; missing=0; assessed=0
    for path in sorted(args.oracle_dir.glob('*.json')):
        state=json.loads(path.read_text())
        baseline_path=args.state_dir/f"{state['instance']}.seed{state['seed']}.json"
        baseline=json.loads(baseline_path.read_text())
        assert baseline['method']==state['baseline_method']
        start=state.get('baseline_snapshot_seconds')
        if start is None: missing+=1; continue
        for r in state['records']:
            if 'query_wall_seconds' not in r: missing+=1; continue
            assessed+=1
            if r['normalized_gain']<=1e-12: continue
            earliest=start+r['query_wall_seconds']
            witnesses=[s for s in baseline['snapshots'] if s.get('seconds',float('inf'))<=earliest
                       and s['cost']<=r['cost']+1e-8]
            witness=min(witnesses,key=lambda s:s['seconds']) if witnesses else None
            rows.append(dict(instance=state['instance'],geometry=state['geometry'],phase=state['phase'],
                             gain=r['normalized_gain'],optimistic_repair_time=earliest,cost=r['cost'],
                             recorded_MDM_dominates=witness is not None,
                             witness=witness and {k:witness[k] for k in ('iteration','cost','seconds')}))
    result=dict(scope='Conservative development diagnostic, sparse checkpoints and historical timings; no confirmed speed claim',
                assessed_queries=assessed,missing_timing_cases_or_queries=missing,
                improving_queries=len(rows),already_matched_by_MDM=sum(r['recorded_MDM_dominates'] for r in rows),rows=rows)
    (HERE/'time_dominance_development.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))

if __name__=='__main__': main()
