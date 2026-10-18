"""Observed one-state scheduling opportunity, never an algorithm speedup bound.

At an input route that later yields an improving executed action, count time
spent on earlier failed actions. A different ordering changes future states;
this diagnostic neither measures nor bounds its whole-trajectory performance.
"""
import hashlib
import json
from pathlib import Path
from collections import defaultdict

HERE=Path(__file__).resolve().parent

def main():
    rows=[];hashes={}
    for path in sorted((HERE/'move_profiles').glob('g_*.json')):
        s=json.loads(path.read_text())
        if 'records' not in s:continue
        epochs=[];current=[];key=None
        for record in s['records']:
            state=tuple(record['before_route'])
            if state!=key and current:epochs.append(current);current=[]
            current.append(record);key=state
            if record['after_cost']<record['before_cost']:
                epochs.append(current);current=[];key=None
        if current:epochs.append(current)
        total=sum(r['operator_seconds'] for r in s['records'])
        useful=[e for e in epochs if e[-1]['after_cost']<e[-1]['before_cost']]
        failed=sum(sum(r['operator_seconds'] for r in e[:-1]) for e in useful)
        terminal=[e for e in epochs if e[-1]['after_cost']>=e[-1]['before_cost'] and
                  set(r['move'] for r in e)==set(range(15))]
        rows.append(dict(instance=s['instance'],geometry=s['geometry'],role=s['role'],calls=len(s['records']),
                         operator_seconds=total,improving_epochs=len(useful),
                         failed_prefix_seconds=failed,failed_prefix_fraction=failed/max(total,1e-12),
                         complete_failed_sweeps=len(terminal),complete_failed_sweep_seconds=sum(r['operator_seconds'] for e in terminal for r in e),
                         sample_capped=s['sample_capped'],profiled_run_completed=s['profiled_run_completed']))
        hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
    grouped=defaultdict(list)
    for r in rows:grouped[r['geometry']].append(r)
    mean=sum(sum(r['failed_prefix_fraction'] for r in rr)/len(rr) for rr in grouped.values())/len(grouped)
    source=Path(__file__).read_bytes();digest=hashlib.sha256(source).hexdigest()
    versions=HERE/'operator_prediction_sources';versions.mkdir(exist_ok=True)
    (versions/(digest+'.py')).write_bytes(source)
    result=dict(scope=__doc__,snapshot=True,source_hashes=hashes,code_sha256=digest,rows=rows,
                geometries=len(grouped),mean_geometry_failed_prefix_fraction=mean,
                definition='Only earlier failures on the same exact input route before an observed successful call; excludes construction and profiling overhead')
    fingerprint=hashlib.sha256(json.dumps(hashes,sort_keys=True).encode()).hexdigest()[:12]
    (HERE/f'operator_opportunity.{fingerprint}.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('source_hashes','rows')},indent=2))

if __name__=='__main__':main()
