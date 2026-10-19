"""Constraint-only reductions preserving all feasible visit subsets (no metric assumption)."""
def reduce_constraints(visits, mandatory, prizes, quota, neighborhoods):
    visits=set(visits); forced=set(mandatory); rows={frozenset(s) for s in neighborhoods}
    assert forced<=visits and all(s<=visits for s in rows)
    if frozenset() in rows:raise ValueError('Infeasible empty coverage row')
    original=len(neighborhoods); history=[]
    while True:
        old=(frozenset(forced),frozenset(rows))
        # Insufficient prize without v proves v mandatory, irrespective of routing.
        total=sum(prizes[v] for v in visits)
        prize_forced={v for v in visits if total-prizes[v]<quota}
        singleton={next(iter(s)) for s in rows if len(s)==1}
        added=(prize_forced|singleton)-forced;forced|=added
        rows={s for s in rows if not s&forced}
        # A subset coverage row implies every superset row.
        rows={s for s in rows if not any(t<s for t in rows)}
        history.append(dict(new_forced=sorted(added),remaining_rows=len(rows)))
        if old==(frozenset(forced),frozenset(rows)):break
    return dict(mandatory=sorted(forced),neighborhoods=[sorted(s) for s in sorted(rows,key=lambda s:(len(s),sorted(s)))],
                removed_rows=original-len(rows),new_mandatory=len(forced-set(mandatory)),history=history)
