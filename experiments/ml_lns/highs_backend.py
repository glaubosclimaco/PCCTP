"""Native installed HiGHS adapter with optional partial integer MIP start.

Kept separate from scipy's bundled HiGHS: version and start are explicit
experimental factors, never silently substituted into historical results.
"""
from types import SimpleNamespace
import highspy
import numpy as np

def solve(obj,integrality,bounds,constraints,seconds,start=None):
    solver=highspy.Highs()
    for name,value in [('output_flag',False),('threads',1),('random_seed',42),
                       ('time_limit',seconds),('mip_rel_gap',0.),('mip_abs_gap',0.)]:
        if solver.setOptionValue(name,value)!=highspy.HighsStatus.kOk:
            raise RuntimeError('Failed HiGHS option '+name)
    lp=highspy.HighsLp(); lp.num_col_=len(obj); lp.num_row_=len(constraints.lb)
    lp.col_cost_=np.asarray(obj,dtype=float)
    lp.col_lower_=bounds.lb; lp.col_upper_=bounds.ub
    lp.row_lower_=constraints.lb; lp.row_upper_=constraints.ub
    lp.integrality_=[highspy.HighsVarType.kInteger if i else highspy.HighsVarType.kContinuous for i in integrality]
    matrix=constraints.A.tocsr()
    lp.a_matrix_.format_=highspy.MatrixFormat.kRowwise
    lp.a_matrix_.start_=matrix.indptr.astype(np.int32)
    lp.a_matrix_.index_=matrix.indices.astype(np.int32)
    lp.a_matrix_.value_=matrix.data
    if solver.passModel(lp)!=highspy.HighsStatus.kOk: raise RuntimeError('Failed HiGHS model')
    if start is not None:
        indices=np.array(sorted(start),dtype=np.int32)
        values=np.array([start[int(i)] for i in indices],dtype=float)
        if solver.setSolution(len(indices),indices,values)!=highspy.HighsStatus.kOk:
            raise RuntimeError('Failed HiGHS MIP start')
    solver.run(); status=solver.getModelStatus(); solution=solver.getSolution(); info=solver.getInfo()
    mapped=0 if status==highspy.HighsModelStatus.kOptimal else (
        1 if status in (highspy.HighsModelStatus.kTimeLimit,highspy.HighsModelStatus.kIterationLimit) else 2)
    return SimpleNamespace(status=mapped,message=solver.modelStatusToString(status),
                           x=np.array(solution.col_value) if solution.value_valid else None,
                           fun=info.objective_function_value,mip_dual_bound=info.mip_dual_bound,
                           version=solver.version())
