"""Compare DP and historical MIP in identical routing-release neighborhoods.

Different execution times and concurrent load prohibit a definitive speed claim;
cross-check every pair of certified optima and every incumbent against bounds.
"""
import json
from pathlib import Path
from dp_repair import DPRepair
from pilot import read_instance,validate

HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[1]

def main():
    out=HERE/'dp_development'; out.mkdir(exist_ok=True); allrows=[]
    for path in sorted((HERE/'routing_release_development').glob('*.json')):
        source=json.loads(path.read_text()); target=out/path.name
        if target.exists():
            allrows.extend(json.loads(target.read_text())['records']);continue
        d=read_instance(ROOT/'data'/source['instance'])
        baseline=json.loads((HERE/'decoder_initialized_states'/f"{source['instance']}.seed42.json").read_text())
        snapshot=next(s for s in baseline['snapshots'] if s['iteration']==1)
        mapping={int(v):i for i,v in enumerate(d['visits'])}
        route=[mapping[v] for v in snapshot['route']]; model=DPRepair(d); rows=[]
        for old in source['records']:
            u=set(old['optional_released'])|set(old['routing_released'])
            result=model.solve(route,u,.5)
            if old['proven_optimal'] and result['proven_optimal']:assert abs(old['cost']-result['cost'])<1e-6
            if old['proven_optimal']:assert result['cost']>=old['cost']-1e-6
            if result['proven_optimal']:assert old['cost']>=result['cost']-1e-6
            assert result['cost']<=validate(d,route)+1e-6
            rows.append(dict(optional_released=old['optional_released'],routing_released=old['routing_released'],
                             mip_cost=old['cost'],mip_optimal=old['proven_optimal'],mip_seconds=old['total_seconds'],dp=result))
        row=dict(instance=source['instance'],geometry=source['geometry'],optional_k=source['optional_k'],
                 routing_k=source['routing_k'],records=rows,dp_optimal=sum(r['dp']['proven_optimal'] for r in rows),
                 dp_best_cost=min(r['dp']['cost'] for r in rows),mip_best_cost=min(r['mip_cost'] for r in rows),
                 dp_strictly_better=sum(r['dp']['cost']<r['mip_cost']-1e-6 for r in rows))
        target.write_text(json.dumps(row,indent=2)+'\n');allrows.extend(rows)
        print(json.dumps({k:v for k,v in row.items() if k!='records'}),flush=True)
    summary=dict(scope='Development same-neighborhood diagnostic; no final speed or novelty claim',queries=len(allrows),
                 dp_optimal=sum(r['dp']['proven_optimal'] for r in allrows),mip_optimal=sum(r['mip_optimal'] for r in allrows),
                 dp_strictly_better=sum(r['dp']['cost']<r['mip_cost']-1e-6 for r in allrows),
                 dp_size_limits=sum(r['dp']['reason']=='component_limit' for r in allrows))
    (HERE/'dp_development_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary),flush=True)

if __name__=='__main__':main()
