"""Independent scalar audit and geometry-level candidate headroom."""
import json,hashlib,statistics
from pathlib import Path
from raw_route_audit import load,check
HERE=Path(__file__).resolve().parent
if __name__=='__main__':
    out=HERE/'candidate_oracle';assert (out/'completion.json').exists();rows=[];checked=0;hashes={}
    for p in sorted(out.glob('g_*.json')):
        r=json.loads(p.read_text());d=load(r['instance']);check(d,r['initial_route'],r['initial_cost']);checked+=1
        initial=HERE/'compiled_control'/f"{r['instance']}.fixed_lkh.json"
        assert hashlib.sha256(initial.read_bytes()).hexdigest()==r['initial_sha256']
        ranks=[e['raw_rank'] for e in r['records']];assert len(ranks)==len(set(ranks)) and set(range(min(8,r['pool_size'])))<=set(ranks)
        assert len(ranks)==min(48,r['pool_size']) and all(0<=i<r['pool_size'] for i in ranks)
        for e in r['records']:
            check(d,e['raw_route'],e['raw_cost']);check(d,e['result']['route'],e['result']['cost']);checked+=2
            assert set(e['raw_route'])==set(e['result']['route']) and e['result']['cost']<=e['raw_cost']
        top=min([r['initial_cost']]+[e['result']['cost'] for e in r['records'] if e['raw_rank']<8]);allbest=min([r['initial_cost']]+[e['result']['cost'] for e in r['records']])
        rows.append(dict(instance=r['instance'],pool_size=r['pool_size'],evaluated=len(ranks),best_top8=top,best_sample=allbest,gain_percent=100*(top-allbest)/top,outside_improves=allbest<top))
        hashes[p.name]=hashlib.sha256(p.read_bytes()).hexdigest()
    assert len(rows)==9
    mean=statistics.mean(r['gain_percent'] for r in rows);wins=sum(r['outside_improves'] for r in rows)
    report=dict(status='PASS',scope='offline TRAIN sampled oracle; not proven neighborhood optimum or end-to-end',rows=rows,routes_checked=checked,wins=wins,mean_gain_percent=mean,gate_passed=wins>=3 and mean>=.2,result_hashes=hashes)
    (HERE/'candidate_oracle_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('rows','result_hashes')}))
