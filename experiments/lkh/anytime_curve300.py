"""Reconstruct descriptive anytime curves from audited absolute-clock events."""
import hashlib,json
from pathlib import Path
from raw_route_audit import load,check

HERE=Path(__file__).resolve().parent
GRID=(1,10,30,60,120,300)


def events(row):
    began=row['native_started'];fallback=row['fallback']
    result=[(began+fallback['seconds'],fallback['cost'],fallback['route'])]
    for event in row['native']['anytime']:
        assert event['seconds']<=row['native']['budget_seconds']+1e-9
        result.append((began+event['seconds'],event['cost'],event['route']))
    if row['post']:
        result.extend((e['native_timestamp'],e['cost'],e['route']) for e in row['post']['trace'])
    return result


def incumbent(stream,deadline):
    eligible=[(c,r) for timestamp,c,r in stream if timestamp<=deadline]
    return min(eligible,key=lambda item:item[0]) if eligible else (None,[])


if __name__=='__main__':
    folder=HERE/'ruin_total300_train'
    audit=json.loads((HERE/'ruin_total300_train_audit.json').read_text())
    assert audit['status']=='PASS'
    rows=[];checked=0;groups={};inputs={}
    for record in audit['rows']:
        path=folder/f"{record['instance']}.seed{record['seed']}.{record['mode']}.json"
        assert hashlib.sha256(path.read_bytes()).hexdigest()==record['result_sha256']
        row=json.loads(path.read_text());name=row['instance']
        if name not in inputs:inputs[name]=load(name)
        d=inputs[name];stream=events(row)
        for timestamp,cost,route in stream:
            check(d,route,cost);checked+=1
        costs={};last=float('inf')
        for horizon in GRID:
            cost,route=incumbent(stream,row['native_started']+horizon)
            if cost is not None:assert cost<=last;last=cost
            costs[str(horizon)]=cost
        assert costs['300']==row['cost']
        result=dict(instance=name,seed=row['seed'],mode=row['mode'],costs=costs)
        rows.append(result);groups.setdefault(name,{})[row['seed'],row['mode']]=result
    curve=[]
    for horizon in GRID:
        for mode in ('hybrid10','hybrid60'):
            effects=[];unavailable=[]
            for name,group in sorted(groups.items()):
                if len(group)!=9:continue
                values={m:[group[s,m]['costs'][str(horizon)] for s in (70,71,72)] for m in ('mdm',mode)}
                if any(v is None for costs in values.values() for v in costs):unavailable.append(name);continue
                baseline=sum(values['mdm'])/3;candidate=sum(values[mode])/3
                effects.append(dict(instance=name,gain_percent=100*(baseline-candidate)/baseline))
            curve.append(dict(horizon=horizon,mode=mode,complete_geometry_pairs=len(effects),
                mean_gain_percent=sum(e['gain_percent'] for e in effects)/len(effects) if effects else None,
                unavailable_due_to_feasibility=unavailable,effects=effects))
    report=dict(status='COMPLETE' if audit['batch_complete'] else 'PARTIAL',
        runs=len(rows),routes_checked=checked,grid=list(GRID),curve=curve,rows=rows,
        scope='descriptive curves only; primary endpoint remains300s; no additional gate or significance tests',
        audit_sha256=hashlib.sha256((HERE/'ruin_total300_train_audit.json').read_bytes()).hexdigest())
    (HERE/'ruin_total300_train_anytime.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('rows','curve')}))
