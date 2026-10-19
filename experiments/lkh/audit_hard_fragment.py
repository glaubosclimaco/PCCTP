"""Audit saved exact repairs and their fixed visits/edges."""
import sys,json,statistics,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'experiments/ml_lns'));sys.path.insert(0,str(Path(__file__).resolve().parent))
from pilot import read_instance,free_edge_set
from raw_route_audit import load,check
HERE=Path(__file__).resolve().parent
if __name__=='__main__':
    folder=HERE/'hard_fragment';assert (folder/'completion.json').exists();groups={};hashes={};proofs=0
    for p in sorted(folder.glob('*.neighbor*.json')):
        r=json.loads(p.read_text());name=r['instance'];d=load(name);rr=r['result'];check(d,r['initial_route'],r['initial_cost']);check(d,rr['route'],rr['cost']);assert rr['cost']<=r['initial_cost']
        parsed=read_instance(ROOT/'data'/name);mapping={int(v):i for i,v in enumerate(parsed['visits'])};initial=[mapping[v] for v in r['initial_route']];final=[mapping[v] for v in rr['route']];released={mapping[v] for v in r['released']}
        assert len(released)==8 and not released.intersection(parsed['mandatory'])
        assert set(initial)-released==set(final)-released
        edges=lambda route:{tuple(sorted((route[i-1],v))) for i,v in enumerate(route)}
        free=free_edge_set(parsed,initial,released,'route_context',seed=42)
        assert not (edges(initial)^edges(final))-free
        assert rr['proven_optimal']==(rr['status']==0);proofs+=rr['proven_optimal']
        hashes[p.name]=hashlib.sha256(p.read_bytes()).hexdigest();groups.setdefault(name,[]).append(r)
    assert len(groups)==9 and all(len(g)==9 for g in groups.values())
    rows=[]
    for name,rs in sorted(groups.items()):
        initial=rs[0]['initial_cost'];best=min(r['result']['cost'] for r in rs)
        rows.append(dict(instance=name,initial_cost=initial,best_cost=best,gain_percent=100*(initial-best)/initial,by_kind={kind:min(r['result']['cost'] for r in rs if r['kind']==kind) for kind in ('random','coverage','prize')}))
    mean=statistics.mean(r['gain_percent'] for r in rows);wins=sum(r['gain_percent']>0 for r in rows)
    report=dict(status='PASS',scope='offline TRAIN exact repairs; optimality certificates only local to fixed neighborhoods',neighborhoods=81,local_optimal_status=proofs,rows=rows,wins=wins,mean_gain_percent=mean,gate_passed=wins>=3 and mean>=.2,result_hashes=hashes)
    (HERE/'hard_fragment_audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k not in ('rows','result_hashes')}))
