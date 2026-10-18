"""Diagnose optional-selection versus mandatory-routing release under fixed total size.

Releasing a mandatory vertex frees incident route edges; its visit stays fixed.
This is a mechanism probe, not evidence of ML or equal-time superiority.
"""
import json
from pathlib import Path
import time
from fragment_repair import FragmentRepair
from pilot import candidates,read_instance,validate

HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[1]

def main():
    out=HERE/'routing_release_development'; out.mkdir(exist_ok=True)
    plan=json.loads((HERE/'decoder_initialized_states/protocol.json').read_text())
    for geometry in ('g_att48','g_berlin52','g_eil51'):
        name=plan['instances'][geometry]
        baseline=json.loads((HERE/'decoder_initialized_states'/f'{name}.seed42.json').read_text())
        assert baseline['role']=='development'
        snapshot=next(s for s in baseline['snapshots'] if s['iteration']==1)
        d=read_instance(ROOT/'data'/name); mapping={int(v):i for i,v in enumerate(d['visits'])}
        route=[mapping[v] for v in snapshot['route']]; initial=validate(d,route)
        model=FragmentRepair(d,'highspy',True,'cuts')
        for optional_k,routing_k in ((8,0),(12,0),(8,4)):
            target=out/f'{name}.k{optional_k}.r{routing_k}.json'
            if target.exists():continue
            records=[]; began=time.perf_counter()
            for kind,u in candidates(d,route,optional_k,12,42):
                anchors=sorted(d['mandatory'],key=lambda v:(min(d['costs'][v,w] for w in u),v))[:routing_k]
                released=u|set(anchors)
                result=model.solve(route,released,.5,'route_context',seed=42)
                records.append(dict(kind=kind,optional_released=sorted(u),routing_released=anchors,**result))
            row=dict(instance=name,geometry=geometry,phase=1,optional_k=optional_k,routing_k=routing_k,
                     initial_cost=initial,elapsed=time.perf_counter()-began,queries=len(records),
                     best_cost=min(r['cost'] for r in records),improving=sum(r['cost']<initial for r in records),
                     optimal=sum(r['proven_optimal'] for r in records),records=records,
                     note='Mandatory anchors remain visited; total release12 compared with optional-only12; developmental mechanism diagnostic')
            target.write_text(json.dumps(row,indent=2)+'\n')
            print(json.dumps({k:v for k,v in row.items() if k!='records'}),flush=True)

if __name__=='__main__':main()
