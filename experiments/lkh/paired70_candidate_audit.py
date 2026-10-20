"""Independent scalar audit for a reference-duration candidate record."""
import math
from raw_route_audit import check


def audit_candidate(r,d,expected_budget):
    assert math.isfinite(expected_budget) and expected_budget>0
    assert r['budget_seconds']==expected_budget
    assert r['mode'] in ('hybrid10','hybrid60')
    initial=min(float(r['mode'][6:]),expected_budget)
    assert r['initialization_budget_seconds']==initial
    began=r['native_started'];deadline=began+expected_budget
    assert math.isfinite(began) and r['wall_seconds']>=0
    checked=0;f=r['fallback'];check(d,f['route'],f['cost']);checked+=1
    assert f['seconds']>=0
    cost,route=(f['cost'],f['route']) if f['seconds']<=expected_budget else (None,[])
    origin='fallback' if route else None
    native=r['native'];assert native['budget_seconds']==initial
    assert native['feasible_before_deadline']==bool(native['anytime'])
    last_time=-1.;last_cost=float('inf')
    for event in native['anytime']:
        assert last_time<=event['seconds']<=initial+1e-9
        assert event['seconds']>=0 and event['cost']<last_cost
        check(d,event['route'],event['cost']);checked+=1
        last_time,last_cost=event['seconds'],event['cost']
        if began+event['seconds']<=deadline and (cost is None or event['cost']<cost):
            cost,route,origin=event['cost'],event['route'],'native'
    assert (native['cost'],native['route'])==(
        (native['anytime'][-1]['cost'],native['anytime'][-1]['route']) if native['anytime'] else (None,[]))
    assert (r['pre_cost'],r['pre_route'])==(cost,route)
    assert r['post_offset']>=0
    post=r['post']
    if post is not None:
        assert route and r['post_offset']<expected_budget and post['source_sha256'] is None
        assert post['initial_cost']==cost
        assert began+r['post_offset']-1e-5<=post['native_started']<=began+r['wall_seconds']
        assert post['lkh_calls']==post['lkh_completed']+post['lkh_timeouts']+post['lkh_failed']
        assert post['attempts']>=post['feasible_repairs']>=post['lkh_calls']
        previous=cost;timestamp=post['native_started']
        for event in post['trace']:
            check(d,event['route'],event['cost']);checked+=1
            assert event['cost']<previous;previous=event['cost']
            assert timestamp<=event['native_timestamp']<=post['native_started']+post['native_elapsed']+1e-5
            timestamp=event['native_timestamp']
            assert abs(timestamp-post['native_started']-event['native_seconds'])<1e-5
            if timestamp<=deadline and (cost is None or event['cost']<cost):
                cost,route,origin=event['cost'],event['route'],'post_event'
        assert post['cost']==previous
        check(d,post['route'],post['cost']);checked+=1
    assert (r['cost'],r['route'],r['origin'])==(cost,route,origin)
    assert r['feasible']==bool(route)
    if route:check(d,route,cost);checked+=1
    return dict(routes_checked=checked,cost=cost,feasible=bool(route))
