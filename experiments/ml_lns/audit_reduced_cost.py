"""Recompute saved certificates with scalar Python integer arithmetic, no LP."""
import json,hashlib
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from experiments.ml_lns.pilot import Repair,read_instance,validate


if __name__=='__main__':
    here=Path(__file__).resolve().parent;path=here/'reduced_cost_results.json'
    results=json.loads(path.read_text());assert len(results['rows'])==11
    checked=0;fixed=0;rows=[]
    for result in results['rows']:
        cert=result['fixing_certificate'];scale=cert['scale']
        d=read_instance(ROOT/'data'/result['instance']);repair=Repair(d)
        matrix=repair.constraints.A.tocsr();lower=list(repair.constraints.lb);upper=list(repair.constraints.ub)
        coefficients=[]
        for r in range(matrix.shape[0]):
            coefficients.append({int(matrix.indices[k]):int(matrix.data[k]) for k in range(matrix.indptr[r],matrix.indptr[r+1])})
            assert all(float(matrix.data[k])==int(matrix.data[k]) for k in range(matrix.indptr[r],matrix.indptr[r+1]))
        root=min(d['mandatory'])
        for cut in result['cuts']:
            side=set(cut['side']);v=cut['representative'];assert root not in side and v in side
            row={i:1 for i,(a,b) in enumerate(repair.edges) if (a in side)!=(b in side)}
            row[repair.m+v]=-2;coefficients.append(row);lower.append(0);upper.append(float('inf'))
        multipliers=dict(cert['multipliers']);r=[scale*int(v) for v in repair.obj]
        assert all(float(v)==int(v) for v in repair.obj)
        bound=0
        for i,row in enumerate(coefficients):
            p=multipliers.get(i,0)
            if not p:continue
            b=lower[i] if p>0 else upper[i]
            assert abs(b)!=float('inf') and b==int(b)
            bound+=p*int(b)
            for j,a in row.items():r[j]-=p*a
        for j,value in enumerate(r):bound+=min(value*int(repair.lb[j]),value*int(repair.ub[j]))
        assert bound==cert['bound_numerator']
        expected=[]
        for j in range(repair.f0):
            if repair.lb[j]==repair.ub[j]:continue
            margin=bound+abs(r[j])-scale*cert['incumbent']
            if margin>0:expected.append(dict(column=j,value=0 if r[j]>0 else 1,reduced_numerator=r[j],margin_numerator=margin))
        assert expected==cert['fixings']
        checked+=1;fixed+=len(expected)
        rows.append(dict(instance=result['instance'],fixed=len(expected),exact_bound_numerator=bound))
    assert results['passing']==sum(row['fixing_passes'] for row in results['rows'])
    assert results['gate_passed']==(results['passing']>=6)
    report=dict(status='PASS',certificates=checked,total_certified_fixings=fixed,
                passing=results['passing'],gate_passed=results['gate_passed'],
                result_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),rows=rows)
    (here/'reduced_cost_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='rows'}))
