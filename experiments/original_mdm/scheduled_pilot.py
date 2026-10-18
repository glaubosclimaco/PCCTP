"""Development fixed-70-iteration pilot, NOT a matched-time benchmark.

Native timers and one-off preparation are charged to full process wall time.
Timeout retains only legacy sparse iteration checkpoints, not best at deadline.
"""
import hashlib
import json
import os
from pathlib import Path
import random
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from experiments.original_mdm.run import run
from experiments.original_mdm.operator_prediction import state_features
from experiments.ml_lns.pilot import read_instance

def main():
    out=HERE/'scheduled_development';out.mkdir(exist_ok=True)
    states=ROOT/'experiments/ml_lns/decoder_initialized_states';cases=[];seen=set()
    for path in sorted(states.glob('g_*.json')):
        b=json.loads(path.read_text())
        if b['role']!='development' or b['geometry'] in seen:continue
        cases.append((path,b));seen.add(b['geometry'])
        if len(cases)==3:break
    # Verify the raw-ID legacy adapter against Python before the policy pilot.
    verified=[]
    for path,b in cases:
        feature_file=out/(b['instance']+'.features.jsonl');feature_file.unlink(missing_ok=True)
        os.environ['PCCTP_SCHEDULE_VERIFY']=str(feature_file.resolve())
        try:
            run(ROOT/'data'/b['instance'],b['seed'],out/(b['instance']+'.adapter_check.json'),
                timeout=90,initialize_decoder_index=True,operator_policy='joint')
        finally:os.environ.pop('PCCTP_SCHEDULE_VERIFY',None)
        d=read_instance(ROOT/'data'/b['instance']);mapping={int(v):i for i,v in enumerate(d['visits'])}
        rows=[json.loads(line) for line in feature_file.read_text().splitlines()];assert rows
        error=0.
        for row in rows:
            expected=state_features(d,[mapping[v] for v in row['route']])
            delta=float(np.max(abs(np.array(expected)-row['features'])));assert delta<1e-9;error=max(error,delta)
        verified.append(dict(instance=b['instance'],states=len(rows),max_error=error))
    methods=['original','static','size','generic','joint','adaptive']
    protocol=dict(scope=__doc__,methods=methods,iterations=70,timeout=90,
                  adapter_verification=verified,
                  cases=[dict(instance=b['instance'],geometry=b['geometry'],role=b['role'],seed=b['seed'],
                              source_sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p,b in cases],
                  code_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in
                               ['scheduled_pilot.py','scheduler.hpp','operator_state.cpp','operator_inference.cpp','run.py','build.py']},
                  build_manifest=json.loads((HERE/'build_scheduled/manifest.json').read_text()))
    fingerprint=hashlib.sha256(json.dumps(protocol,sort_keys=True).encode()).hexdigest()[:12]
    (out/f'protocol.{fingerprint}.json').write_text(json.dumps(protocol,indent=2)+'\n')
    for path,b in cases:
        order=methods[:];random.Random(b['instance']).shuffle(order)
        for policy in order:
            target=out/f"{b['instance']}.seed{b['seed']}.{policy}.{fingerprint}.json"
            if target.exists():continue
            result=run(ROOT/'data'/b['instance'],b['seed'],target,timeout=90,
                       initialize_decoder_index=True,operator_policy=policy)
            if policy=='original' and result['completed']:assert result['cost']==b['cost'],'Original policy changed baseline cost'
            print(json.dumps(dict(instance=b['instance'],policy=policy,cost=result['cost'],
                                  seconds=result['process_wall_seconds'],completed=result['completed'],
                                  baseline_cost=b['cost'],protocol=fingerprint)),flush=True)

if __name__=='__main__':main()
