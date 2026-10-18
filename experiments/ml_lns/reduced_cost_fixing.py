"""Exact integer weak-duality certificates built from rounded LP multipliers."""
import json,time
from pathlib import Path
import sys
import numpy as np
from scipy.sparse import csr_matrix,csc_matrix
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
SCALE=1000000


def integer_array(values):
    values=np.asarray(values,dtype=float)
    assert np.all(np.isfinite(values)) and np.all(values==np.rint(values))
    assert np.max(np.abs(values),initial=0)<2**52
    return np.rint(values).astype(np.int64)


def certificate(matrix,c,lo,up,rowlo,rowup,duals,binary_count,incumbent,initial):
    matrix=matrix.tocsr().copy();matrix.data=integer_array(matrix.data)
    c=integer_array(c);lo=integer_array(lo);up=integer_array(up)
    assert np.all(np.isfinite(duals))
    p=np.rint(np.asarray(duals)*SCALE)
    assert np.max(abs(p),initial=0)<2**52
    p=p.astype(np.int64)
    b=[]
    for i in range(len(p)):
        if not np.isfinite(rowlo[i]) and p[i]>0:p[i]=0
        if not np.isfinite(rowup[i]) and p[i]<0:p[i]=0
        value=rowlo[i] if p[i]>=0 else rowup[i]
        b.append(0 if p[i]==0 else int(integer_array([value])[0]))
    # Bound every intermediate sparse dot-product before using int64.
    assert sum(abs(int(v)) for v in matrix.data)<2**60
    colnorm=np.asarray(abs(matrix).sum(axis=0)).ravel()
    assert max((int(v) for v in colnorm),default=0)*max((abs(int(v)) for v in p),default=0)+SCALE*max((abs(int(v)) for v in c),default=0)<2**60
    reduced=SCALE*c-matrix.T@p
    bound=sum(int(v)*w for v,w in zip(p,b))+sum(min(int(r)*int(l),int(r)*int(u)) for r,l,u in zip(reduced,lo,up))
    incumbent_int=int(integer_array([incumbent])[0]);assert bound<=SCALE*incumbent_int
    fixes=[];free=0
    for j in range(binary_count):
        if lo[j]==up[j]:continue
        assert lo[j]==0 and up[j]==1;free+=1
        r=int(reduced[j]);margin=bound+abs(r)-SCALE*incumbent_int
        if margin>0:
            value=0 if r>0 else 1
            assert initial[j]==value
            fixes.append(dict(column=j,value=value,reduced_numerator=r,margin_numerator=margin))
    return dict(scale=SCALE,bound_numerator=bound,incumbent=incumbent_int,
                multipliers=[[i,int(v)] for i,v in enumerate(p) if v],
                initial_free_binary=free,fixings=fixes,fixed_fraction=len(fixes)/max(1,free))


def from_highs(h,repair,incumbent,initial):
    import highspy
    lp=h.getLp();a=lp.a_matrix_
    constructor=csr_matrix if a.format_==highspy.MatrixFormat.kRowwise else csc_matrix
    matrix=constructor((a.value_,a.index_,a.start_),shape=(lp.num_row_,lp.num_col_))
    return certificate(matrix,lp.col_cost_,lp.col_lower_,lp.col_upper_,lp.row_lower_,lp.row_upper_,
                       h.getSolution().row_dual,repair.f0,incumbent,initial)


if __name__=='__main__':
    from experiments.ml_lns.gsec_root import probe
    directory=ROOT/'experiments/original_mdm/move_profiles_monotonic_elite_bounded'
    rows=[]
    for path in sorted(directory.glob('*.profiled_baseline.json')):
        profile=json.loads(path.read_text())
        if not profile['completed'] or profile['iterations']!=70:continue
        result=probe(path,with_fixing=True);rows.append(result)
        cert=result.get('fixing_certificate')
        print(json.dumps(dict(instance=result['instance'],seconds=result['seconds'],fraction=None if cert is None else cert['fixed_fraction'],passes=result.get('fixing_passes',False))),flush=True)
    passing=sum(r.get('fixing_passes',False) for r in rows)
    Path(__file__).with_name('reduced_cost_results.json').write_text(json.dumps(dict(rows=rows,passing=passing,gate_passed=len(rows)==11 and passing>=6),indent=2)+'\n')
