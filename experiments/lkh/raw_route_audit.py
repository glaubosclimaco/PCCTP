"""Independent scalar route audit without project readers or NumPy."""
import json,hashlib,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent

def load(name):
    p=ROOT/'data'/name;t=p.read_text().split();n=int(t[0]);assert len(t)==4+2*n+n*n
    return dict(n=n,quota=float(t[1]),radius=float(t[2]),header_target=float(t[3]),types=list(map(int,t[4:4+n])),prizes=list(map(float,t[4+n:4+2*n])),matrix=list(map(float,t[4+2*n:])),input_hash=hashlib.sha256(p.read_bytes()).hexdigest())

def check(d,route,cost):
    n=d['n'];s=set(route)
    assert len(route)==len(s)>=3 and all(type(v)==int and 0<=v<n for v in s)
    assert all(d['types'][v]!=2 for v in s)
    assert all(v in s for v in range(n) if d['types'][v]==1)
    assert sum(d['prizes'][v] for v in s)>=d['quota']
    assert all(any(d['matrix'][v*n+w]<=d['radius'] for v in s) for w in range(n) if d['types'][w]==2)
    actual=sum(d['matrix'][route[i-1]*n+v] for i,v in enumerate(route))
    assert actual==cost,(actual,cost)
    return actual

if __name__=='__main__':
    assert len(sys.argv)==2 and sys.argv[1] in ('support_total','support_diagnostic','support_common_fallback','compiled_control','direct_support')
    folder=HERE/sys.argv[1];instances={};hashes={};checked=0;empty=0
    for p in sorted(folder.glob('*.json')):
        if p.name in ('protocol.json','completion.json'):continue
        r=json.loads(p.read_text());name=r['instance'];d=instances.setdefault(name,load(name))
        hashes[p.name]=hashlib.sha256(p.read_bytes()).hexdigest()
        if r['cost'] is None:
            assert r['route']==[];empty+=1
        else:check(d,r['route'],r['cost']);checked+=1
        for e in r.get('anytime',[])+r.get('trace',[]):check(d,e['route'],e['cost']);checked+=1
        if r.get('fallback'):
            f=r['fallback'];check(d,f['route'],f['cost']);checked+=1
            assert f['accepted']==(f['accepted_seconds']<=15.)
            if r.get('pre_post_cost') is not None:
                check(d,r['pre_post_route'],r['pre_post_cost']);checked+=1
        if r.get('post'):
            post=r['post'];check(d,post['route'],post['cost']);checked+=1
            for e in post['trace']:check(d,e['route'],e['cost']);checked+=1
    report=dict(status='PASS',scope='independent raw-input scalar route feasibility and cost; not deadline audit',batch=folder.name,batch_complete=(folder/'completion.json').exists(),files=len(hashes),routes_checked=checked,empty_incumbent_records=empty,input_hashes={n:d['input_hash'] for n,d in instances.items()},result_hashes=hashes)
    (HERE/(folder.name+'_raw_audit.json')).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('input_hashes','result_hashes')}))
