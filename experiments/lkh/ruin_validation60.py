"""Frozen validation-only run of the corrected hybrid; no test data access."""
import sys,json,hashlib,random
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import ruin_total60_absolute as core
if __name__=='__main__':
 training=json.loads((HERE/'ruin_total60_absolute_summary.json').read_text());assert training['status']=='COMPLETE' and training['gate_passed']
 component=json.loads((HERE/'fixed_support_ablation_audit.json').read_text());assert component['batch_complete'] and component['status']=='PASS' and component['component_gate_passed']
 split_path=ROOT/'experiments/ml_lns/decoder_initialized_states/protocol.json';split=json.loads(split_path.read_text())
 cases=[dict(geometry=g,instance=split['instances'][g]) for g,role in sorted(split['roles'].items()) if role=='validation'];assert len(cases)==7
 out=HERE/'ruin_validation60';out.mkdir(exist_ok=True);core.OUT=out
 names=['experiments/lkh/ruin_validation60.py','experiments/lkh/ruin_total60_absolute.py','experiments/lkh/ruin_repair_absolute.py','experiments/lkh/feasible_fallback.py','experiments/lkh/screen.py','experiments/original_mdm/timed_experiment.py','experiments/ml_lns/pilot.py','experiments/ml_branching/oracle_probe.py']
 binaries=[ROOT/'experiments/original_mdm/build_timed_elite_bounded/original_mdm',ROOT/'experiments/original_mdm/build_timed_elite_bounded/FPMax/fim_maximal',ROOT/'experiments/original_mdm/build_timed_elite_bounded/FPMax/fpmax_hnmp',Path('/private/tmp/pcctp-lkh3-research/LKH')]
 protocol=dict(scope='validation screening; same geometries used by earlier different method; not independent confirmation',cases=cases,seeds=list(range(45,50)),budget=60,split_sha256=hashlib.sha256(split_path.read_bytes()).hexdigest(),source_hashes={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names},binary_hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in binaries})
 p=out/'protocol.json'
 if p.exists():assert json.loads(p.read_text())==protocol
 else:p.write_text(json.dumps(protocol,indent=2)+'\n')
 for c in cases:
  assert split['roles'][c['geometry']]=='validation'
  for seed in range(45,50):
   order=['mdm','hybrid'];random.Random(c['instance']+str(seed)).shuffle(order)
   for mode in order:
    row=core.run(c['instance'],seed,mode);print(json.dumps({k:row[k] for k in ('instance','seed','mode','cost','wall_seconds')}),flush=True)
 (out/'completion.json').write_text(json.dumps(dict(runs=70))+'\n')
