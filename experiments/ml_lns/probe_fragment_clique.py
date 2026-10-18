"""Same-neighborhood exact repair comparison using a complete boundary graph."""
import json
import argparse
from pathlib import Path
from dp_repair import DPRepair
from fragment_repair import FragmentRepair
from pilot import read_instance,validate

HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[1]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,default=HERE/'fragment_clique_development')
    args=parser.parse_args()
    out=args.output_dir; out.mkdir(exist_ok=True); allrows=[]
    for path in sorted((HERE/'routing_release_development').glob('*.json')):
        source=json.loads(path.read_text()); target=out/path.name
        if target.exists():allrows.extend(json.loads(target.read_text())['records']);continue
        d=read_instance(ROOT/'data'/source['instance'])
        baseline=json.loads((HERE/'decoder_initialized_states'/f"{source['instance']}.seed42.json").read_text())
        snapshot=next(s for s in baseline['snapshots'] if s['iteration']==1)
        mapping={int(v):i for i,v in enumerate(d['visits'])}
        route=[mapping[v] for v in snapshot['route']]
        dp=DPRepair(d); mip=FragmentRepair(d,'highspy',True,'cuts'); rows=[]
        for i,record in enumerate(source['records']):
            u=set(record['optional_released'])|set(record['routing_released'])
            # Alternate backend order to reduce a fixed ordering effect.
            if i%2:
                a=dp.solve(route,u,.5,'fragment_clique'); b=mip.solve(route,u,.5,'fragment_clique')
            else:
                b=mip.solve(route,u,.5,'fragment_clique'); a=dp.solve(route,u,.5,'fragment_clique')
            if a['proven_optimal'] and b['proven_optimal']:assert abs(a['cost']-b['cost'])<1e-6
            if a['proven_optimal']:assert b['cost']>=a['cost']-1e-6
            if b['proven_optimal']:assert a['cost']>=b['cost']-1e-6
            assert a['contracted_nodes']<=2*len(u)
            rows.append(dict(optional_released=record['optional_released'],routing_released=record['routing_released'],dp=a,mip=b))
        row=dict(instance=source['instance'],geometry=source['geometry'],optional_k=source['optional_k'],routing_k=source['routing_k'],
                 initial_cost=validate(d,route),records=rows,dp_optimal=sum(r['dp']['proven_optimal'] for r in rows),
                 mip_optimal=sum(r['mip']['proven_optimal'] for r in rows),
                 dp_best_cost=min(r['dp']['cost'] for r in rows),mip_best_cost=min(r['mip']['cost'] for r in rows))
        target.write_text(json.dumps(row,indent=2)+'\n');allrows.extend(rows)
        print(json.dumps({k:v for k,v in row.items() if k!='records'}),flush=True)
    summary=dict(scope='Development same-neighborhood comparison; no final speed or novelty claim',queries=len(allrows),
                 dp_optimal=sum(r['dp']['proven_optimal'] for r in allrows),mip_optimal=sum(r['mip']['proven_optimal'] for r in allrows),
                 dp_strictly_better=sum(r['dp']['cost']<r['mip']['cost']-1e-6 for r in allrows),
                 size_limited=sum(r['dp']['reason']=='component_limit' for r in allrows))
    summary_path=HERE/('fragment_clique_summary.json' if out.name=='fragment_clique_development' else out.name+'.summary.json')
    summary_path.write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary),flush=True)

if __name__=='__main__':main()
