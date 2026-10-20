"""Synthetic absolute-clock eligibility checks; no benchmark or optimizer executed."""
import hashlib,json,math
from pathlib import Path
from paired70_hybrid import select_events

if __name__=='__main__':
    rows=[]
    cases=[('exact_boundary',10,[dict(native_timestamp=110.,cost=8,route=[1,2,3])],8),
           ('late_better_rejected',10,[dict(native_timestamp=109.,cost=8,route=[1,2,3]),dict(native_timestamp=110.000001,cost=1,route=[1,3,2])],8),
           ('no_initial_feasible',None,[dict(native_timestamp=105.,cost=8,route=[1,2,3])],8),
           ('late_only',None,[dict(native_timestamp=111.,cost=8,route=[1,2,3])],None),
           ('fractional_horizon',10,[dict(native_timestamp=100.125,cost=8,route=[1,2,3]),dict(native_timestamp=100.126,cost=1,route=[1,3,2])],8)]
    for label,initial,events,expected in cases:
        budget=.125 if label=='fractional_horizon' else 10.
        cost,route,origin=select_events(initial,[0,1,2] if initial is not None else [],events,100.,budget)
        assert cost==expected,(label,cost,expected)
        rows.append(dict(case=label,expected_cost=expected,observed_cost=cost))
    rejected=0
    for budget in (0.,-1.,math.inf,math.nan):
        try:select_events(10,[0,1,2],[],100.,budget)
        except AssertionError:rejected+=1
    for timestamp in (99.,math.inf,math.nan):
        try:select_events(10,[0,1,2],[dict(native_timestamp=timestamp,cost=1,route=[1,2,3])],100.,10.)
        except AssertionError:rejected+=1
    assert rejected==7
    source=Path(__file__).with_name('paired70_hybrid.py')
    report=dict(status='PASS',valid_scenarios=len(rows),invalid_inputs_rejected=rejected,
        rows=rows,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        scope='finite selector checks only; main ten-seed study not launched')
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='rows'}))
