"""Strict-time baseline for ALL nine incomplete TRAIN profiles; no profiling."""
import json,re,hashlib
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from experiments.original_mdm.timed_experiment import run_one
from experiments.ml_lns.pilot import read_instance,validate
HERE=Path(__file__).resolve().parent
OUT=HERE/'hard_training_baseline'


def published_rows():
    lines=(ROOT/'artigos/mdm/main_review.tex').read_text().splitlines()
    start=next(i for i,s in enumerate(lines) if r'\label{tab:mdm-grasp}' in s)
    end=next(i for i in range(start,len(lines)) if r'\end{longtable}' in lines[i])
    rows={}
    for index in range(start,end):
        line=lines[index].replace(r'\_','_')
        cells=[s.replace(r'\textbf{','').replace('}','').strip() for s in line.split('&')]
        if len(cells)!=9 or not re.fullmatch(r'\w+_R\d+_T\d+_W\d+_\d+',cells[0]):continue
        rows[cells[0]]=dict(best=float(cells[4]),mean=float(cells[5]),
                             published_average_seconds=float(cells[6]),source_line=index+1)
    return rows


if __name__=='__main__':
    OUT.mkdir(exist_ok=True);published=published_rows();profiles=[]
    for path in sorted((HERE/'move_profiles_monotonic_elite_bounded').glob('*.profiled_baseline.json')):
        profile=json.loads(path.read_text())
        if not profile['completed']:profiles.append(profile)
    assert len(profiles)==9
    protocol=dict(budget_seconds=10.,seed=42,policy='original',
                  cohort_rule='All nine TRAIN geometries with incomplete observational move profiles',
                  cases=[p['instance'] for p in profiles],
                  binary_sha256=hashlib.sha256((HERE/'build_timed_elite_bounded/original_mdm').read_bytes()).hexdigest(),
                  publication_source='artigos/mdm/main_review.tex, table tab:mdm-grasp',publication_times_comparable=False)
    protocol_path=OUT/'protocol.json'
    if protocol_path.exists():assert json.loads(protocol_path.read_text())==protocol
    else:protocol_path.write_text(json.dumps(protocol,indent=2)+'\n')
    rows=[]
    for profile in profiles:
        name=profile['instance'];alias=name[2:].replace('_VT','_R');reference=published.get(alias)
        d=read_instance(ROOT/'data'/name);mapping={int(v):i for i,v in enumerate(d['visits'])}
        assert validate(d,[mapping[v] for v in profile['route']])==profile['cost']
        target=OUT/f'{name}.seed42.original.json'
        timed=json.loads(target.read_text()) if target.exists() else run_one(name,42,'original',10.,target,True)
        row=dict(instance=name,geometry=profile['geometry'],profile_completed=False,
                 profile_iterations=profile['iterations'],profile_checkpoint_cost=profile['cost'],
                 profile_checkpoint_seconds=profile['seconds'],profile_wall_timeout=profile['wall_timeout_seconds'],
                 profile_instrumented=profile['move_profiling'],published=reference,
                 timed_cost=timed['cost'],feasible_before_deadline=timed['feasible_before_deadline'],
                 timed_overrun=timed['overrun_seconds'],hard_timeout=timed['hard_timeout'])
        row['gap_to_published_best_percent']=None if reference is None or timed['cost'] is None else 100*(timed['cost']/reference['best']-1)
        rows.append(row)
        (OUT/'summary.json').write_text(json.dumps(dict(rows=rows,completed_cases=len(rows),target_cases=9),indent=2)+'\n')
        print(json.dumps(row),flush=True)
