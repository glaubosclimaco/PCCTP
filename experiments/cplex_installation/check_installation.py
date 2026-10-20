"""Small isolated CPLEX/DOcplex installation check; no PCCTP experiment changes."""
import importlib.metadata,json,math,platform,sys
from pathlib import Path
import cplex
from docplex.mp.model import Model


def quiet(model):
    for setter in (model.set_log_stream,model.set_results_stream,
                   model.set_warning_stream,model.set_error_stream):setter(None)
    model.parameters.threads.set(1)
    model.parameters.timelimit.set(5)


if __name__=='__main__':
    with cplex.Cplex() as model:
        quiet(model)
        model.objective.set_sense(model.objective.sense.maximize)
        model.variables.add(obj=[3.,2.],types='BB',names=['x','y'])
        model.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=[0,1],val=[1.,1.])],
                                     senses='L',rhs=[1.])
        model.solve();value=model.solution.get_objective_value()
        assert math.isclose(value,3.)
        low_level=dict(status=model.solution.get_status_string(),objective=value,
                       values=model.solution.get_values(),engine_version=model.get_version())
    with Model(name='installation_check') as model:
        x=model.binary_var(name='x');y=model.binary_var(name='y')
        model.maximize(3*x+2*y);model.add_constraint(x+y<=1)
        model.parameters.threads=1;model.parameters.timelimit=5
        solution=model.solve(log_output=False)
        assert solution is not None and math.isclose(solution.objective_value,3.)
        modeling=dict(status=str(model.solve_details.status),objective=solution.objective_value)
    limit=dict(variables=1001,constraints=1,restricted_at_1000=None)
    with cplex.Cplex() as model:
        quiet(model)
        model.variables.add(obj=[1.]*1001,types='B'*1001)
        model.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=list(range(1001)),val=[1.]*1001)],
                                     senses='L',rhs=[500.])
        try:
            model.solve();limit.update(restricted_at_1000=False,status=model.solution.get_status_string())
        except cplex.exceptions.CplexSolverError as error:
            message=str(error)
            if '1016' not in message:raise
            limit.update(restricted_at_1000=True,error_code=1016,message=message)
    report=dict(status='PASS',python=sys.executable,python_version=platform.python_version(),
        architecture=platform.machine(),cplex=importlib.metadata.version('cplex'),
        docplex=importlib.metadata.version('docplex'),low_level=low_level,modeling=modeling,
        size_limit_probe=limit,scope='isolated Python runtime; not full Optimization Studio IDE or C++ SDK')
    Path(__file__).with_name('installation_check.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report))
