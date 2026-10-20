"""Passive native mining-exposure bounds; never changes an optimizer or its logs."""
import hashlib
import json
import sys
from pathlib import Path

HERE=Path(__file__).resolve().parent


def profile(folder):
    rows=[]
    for path in sorted(folder.glob('g_*.json')):
        if '.native.' in path.name:continue
        row=json.loads(path.read_text());native=row['native']
        budget=native['budget_seconds']
        exports=[];errors=[]
        raw=path.with_name(path.stem+'.native.raw')
        for source in (raw,Path(str(raw)+'.checkpoint')):
            if not source.exists():continue
            try:payload=json.loads(source.read_text())
            except (UnicodeDecodeError,json.JSONDecodeError) as error:
                errors.append(dict(file=source.name,error=str(error)));continue
            exports.append((source,payload))
        selected=native.get('payload_source')
        chosen=next((d for p,d in exports if p.name==selected),None)
        if selected is not None:assert chosen is not None
        lower=max([int(d['mining_calls']) for _,d in exports
                   if d['seconds']<=budget]+[0])
        # A completed cumulative count bounds timely calls, but can include
        # work in a cooperatively overrunning final iteration.
        finals=[d for _,d in exports if d['completed']]
        upper=max(int(d['mining_calls']) for d in finals) if finals else None
        if upper is not None:assert lower<=upper
        log=path.with_name(path.stem+'.native.log')
        returns=log.read_text(errors='replace').count('Minerador executado.') if log.exists() else None
        rows.append(dict(instance=row['instance'],seed=row['seed'],mode=row['mode'],
            native_budget_seconds=budget,
            observed_mining_calls=chosen['mining_calls'] if chosen else None,
            observed_iterations=chosen['iterations'] if chosen else None,
            final_export_present=bool(finals),timely_mining_calls_lower_bound=lower,
            timely_mining_calls_upper_bound=upper,logged_mining_returns=returns,
            completed_call_return_count_agrees=(returns==upper) if upper is not None and returns is not None else None,
            payload_read_errors=errors,
            export_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p,_ in exports}))
    modes={}
    for mode in sorted({r['mode'] for r in rows}):
        group=[r for r in rows if r['mode']==mode]
        modes[mode]=dict(runs=len(group),
            runs_with_observed_call=sum((r['observed_mining_calls'] or 0)>0 for r in group),
            runs_with_proven_timely_call=sum(r['timely_mining_calls_lower_bound']>0 for r in group),
            runs_with_final_export=sum(r['final_export_present'] for r in group),
            runs_with_unknown_upper_bound=sum(r['timely_mining_calls_upper_bound'] is None for r in group),
            call_return_count_mismatches=sum(r['completed_call_return_count_agrees'] is False for r in group))
    return dict(status='COMPLETE' if (folder/'completion.json').exists() else 'PARTIAL',
        runs=len(rows),modes=modes,rows=rows,
        scope='passive cumulative native exports only; un-timestamped calls are not automatically timely; exposure is not mining effectiveness')


if __name__=='__main__':
    batch=sys.argv[1] if len(sys.argv)>1 else 'ruin_total300_train'
    assert batch in ('ruin_total300_train','ruin_total60_absolute','ruin_validation60')
    report=profile(HERE/batch)
    (HERE/(batch+'_mining_exposure_bounds.json')).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='rows'}))
