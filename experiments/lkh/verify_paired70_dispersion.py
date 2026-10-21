"""Synthetic verification only; does not launch optimization or read PCCTP inputs."""
import hashlib,json,math
from copy import deepcopy
from pathlib import Path
from paired70_dispersion import report
HERE=Path(__file__).resolve().parent
protocol=dict(seeds=list(range(81,91)),cases=[dict(instance=f'toy{i}',geometry=f'toy{i}') for i in range(9)])
rows=[dict(instance='toy0',seed=s,status='PAIR_COMPLETED',reference_cost=100+s-81,duration=10+s-81,candidate_costs=dict(hybrid10=99+s-81,hybrid60=100+s-81)) for s in range(81,91)]
a=dict(status='PASS',rows=rows);r=report(a,protocol);g=r['geometries'][0]
assert g['reference_cost']['mean']==g['reference_cost']['median']==104.5
assert abs(g['reference_cost']['sample_sd']-math.sqrt(55/6))<1e-12
assert g['candidates']['hybrid10']['wins']==10
assert g['candidates']['hybrid10']['paired_cost_reduction']['sample_sd']==0
assert g['candidates']['hybrid60']['ties']==10
assert r['completed_geometries']==1 and len(r['geometries'])==9
b=deepcopy(a);b['rows'].pop();assert 'reference_cost' not in report(b,protocol)['geometries'][0]
b=deepcopy(a);b['rows'][0]['status']='REFERENCE_INCOMPLETE';assert not report(b,protocol)['geometries'][0]['complete']
for altered in [dict(status='PASS',rows=rows+[rows[0]]),dict(status='FAIL',rows=rows)]:
    try:report(altered,protocol)
    except AssertionError:pass
    else:raise AssertionError('Invalid audit accepted')
result=dict(status='PASS_SYNTHETIC_ONLY',scope='Hand-calculated ten-seed toy data; not PCCTP outcomes. No primary gate change.',checks=['mean/median104.5','sample variance55/6','constant paired gain SD0','10ties','9geometry denominator','incomplete suppressed','failed reference suppressed','duplicate rejected','failed audit rejected'],source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (HERE/'paired70_dispersion.py',Path(__file__))})
(HERE/'paired70_dispersion_verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
