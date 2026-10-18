"""Descriptive paired comparisons by geometry; no independence or novelty claim."""
import argparse
import hashlib
import json
from pathlib import Path
from collections import defaultdict
import numpy as np

HERE=Path(__file__).resolve().parent

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('protocol');p.add_argument('--partial',action='store_true');args=p.parse_args()
    out=HERE/'timed_development';plan=json.loads((out/f'protocol.{args.protocol}.json').read_text())
    invalid_path=HERE/'invalid_timing_protocols.json'
    if invalid_path.exists() and args.protocol in json.loads(invalid_path.read_text()):
        p.error('Protocol invalidated by clock audit; preserved for diagnostics only')
    completed=(out/f'completion.{args.protocol}.json').exists()
    if not completed and not args.partial:p.error('Collection unfinished; use --partial for an explicitly partial snapshot')
    records={};sources={}
    for path in sorted(out.glob(f'g_*.{args.protocol}.json')):
        r=json.loads(path.read_text());records[r['instance'],r['seed'],r['operator_policy']]=r
        sources[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
    expected=len(plan['cases'])*len(plan['seeds'])*len(plan['methods'])
    if completed:assert len(records)==expected
    comparisons={}
    for control in ['original','generic','adaptive']:
        for method in plan['methods']:
            if method==control:continue
            by_geometry=defaultdict(list);wins=ties=losses=missing=unavailable=0
            for case in plan['cases']:
                for seed in plan['seeds']:
                    a=records.get((case['instance'],seed,method));b=records.get((case['instance'],seed,control))
                    if a is None or b is None:missing+=1;continue
                    if a['cost'] is None or b['cost'] is None:unavailable+=1;continue
                    delta=(a['cost']-b['cost'])/max(1.,b['cost'])
                    by_geometry[case['geometry']].append(delta)
                    wins+=a['cost']<b['cost'];ties+=a['cost']==b['cost'];losses+=a['cost']>b['cost']
            means={g:float(np.mean(values)) for g,values in by_geometry.items()}
            comparisons[method+'_vs_'+control]=dict(mean_geometry_relative_cost_difference=float(np.mean(list(means.values()))) if means else None,
                geometry_differences=means,paired_runs_better=wins,tied=ties,worse=losses,missing=missing,unavailable=unavailable)
    result=dict(scope=__doc__,protocol=args.protocol,partial=not completed,expected=expected,observed=len(records),
                feasible=sum(r['cost'] is not None for r in records.values()),
                hard_timeouts=sum(r['hard_timeout'] for r in records.values()),
                abnormal_exits=sum(r['exit_code'] not in (0,None) for r in records.values()),
                source_hashes=sources,comparisons=comparisons,
                interpretation='Negative difference favors the named method. Seeds and variants are dependent within geometry; no hypothesis tests. Unavailable/missing pairs reported explicitly.')
    suffix='.partial_snapshot' if not completed else ''
    (HERE/f'timed_summary.{args.protocol}{suffix}.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('source_hashes','comparisons')},indent=2))
    for name in ['joint_vs_original','joint_vs_generic','joint_vs_adaptive']:print(name,json.dumps(comparisons[name]))

if __name__=='__main__':main()
