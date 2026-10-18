"""Paired total-time comparison of fixed root GSEC strengthening."""
import time,json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from experiments.ml_lns.gsec_root import probe
from experiments.ml_lns.highs_warm_start import solve,read_instance


if __name__=='__main__':
    profiles=ROOT/'experiments/original_mdm/move_profiles_monotonic_elite_bounded'
    rows=[];output=Path(__file__).with_name('gsec_integer_results.json')
    paths=[p for p in sorted(profiles.glob('*.profiled_baseline.json')) if json.loads(p.read_text())['completed']]
    for k,path in enumerate(paths):
        profile=json.loads(path.read_text());assert profile['iterations']==70
        row=dict(instance=profile['instance'],geometry=profile['geometry'])
        order=('plain','cuts') if k%2==0 else ('cuts','plain')
        row['order']=order
        for mode in order:
            began=time.perf_counter();root=probe(path) if mode=='cuts' else None
            profile=json.loads(path.read_text());d=read_instance(ROOT/'data'/profile['instance'])
            mapping={int(v):i for i,v in enumerate(d['visits'])}
            route=[mapping[v] for v in profile['route']]
            result=solve(d,route,True,began=began,root_cuts=None if root is None else root['cuts'])
            result['separation_seconds']=0. if root is None else root['seconds']
            result['root_cut_count']=0 if root is None else root['cut_count']
            row[mode]=result
        rows.append(row)
        gain=sum(r['cuts']['improved'] for r in rows);basegain=sum(r['plain']['improved'] for r in rows)
        extra=sum(int(r['cuts']['proven_optimal'])-int(r['plain']['proven_optimal']) for r in rows)
        no_loss=all(r['cuts']['cost']<=r['plain']['cost']+1e-6 for r in rows)
        result=dict(rows=rows,cut_gains=gain,plain_gains=basegain,net_extra_proofs=extra,
                    gate_passed=len(rows)==11 and ((gain>=2 and gain>basegain) or (extra>=2 and no_loss)))
        output.write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps(dict(instance=row['instance'],plain_cost=row['plain']['cost'],
                    cut_cost=row['cuts']['cost'],plain_proven=row['plain']['proven_optimal'],
                    cut_proven=row['cuts']['proven_optimal'],plain_seconds=row['plain']['seconds'],
                    cut_seconds=row['cuts']['seconds'],separation=row['cuts']['separation_seconds'])),flush=True)
