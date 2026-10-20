"""Scalar original-input, matched destruction, slots and absolute-deadline audit."""
import hashlib
import json
import random
from pathlib import Path
from raw_route_audit import load, check

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


if __name__ == '__main__':
    folder = HERE/'joint_insertion_train'
    if not (folder/'protocol.json').exists():
        print('NO_RESULTS')
        raise SystemExit(0)
    protocol = json.loads((folder/'protocol.json').read_text())
    for source, hashed in protocol['source_hashes'].items():
        assert digest(ROOT/source) == hashed
    rows = []
    checked = 0
    data = {}
    keys = set()
    for path in sorted(folder.glob('g_*.json')):
        r = json.loads(path.read_text())
        name, seed, sample, arm = (r[k] for k in ('instance','seed','sample','arm'))
        assert seed in (42,43,44) and sample in range(4)
        assert arm in ('greedy','mip_nearest','mip_random')
        key = (name,seed,sample,arm)
        assert key not in keys
        keys.add(key)
        if name not in data:
            data[name] = load(name)
        d = data[name]
        source = HERE/'ruin_total60_absolute'/r['input_file']
        assert digest(source) == r['input_sha256'] == protocol['input_hashes'][source.name]
        start = json.loads(source.read_text())
        assert start['instance'] == name and start['seed'] == seed and start['mode']=='hybrid'
        original = start['pre_route']
        check(d, original, start['pre_cost'])
        checked += 1
        assert r['fallback_cost'] == start['pre_cost']
        label = f'{name}:{seed}:{sample}:20261020'
        rng = random.Random(int(hashlib.sha256(label.encode()).hexdigest(),16))
        # Project visit indices are sorted original non-demand indices.
        visits = [v for v,t in enumerate(d['types']) if t != 2]
        optional = sorted(v for v in original if d['types'][v] != 1)
        removed = set(rng.sample(optional,min(6,len(optional),len(original)-3)))
        assert set(r['removed']) == removed
        retained = [v for v in original if v not in removed]
        assert r['retained'] == retained and len(retained) >= 3
        assert abs(r['deadline']-r['native_started']-1.) < 1e-8
        best_cost, best_route = start['pre_cost'], original
        previous = r['native_started']
        for event in r['events']:
            assert previous <= event['native_timestamp'] <= r['native_started']+r['native_elapsed']
            previous = event['native_timestamp']
            check(d, event['route'], event['cost'])
            checked += 1
            if event['native_timestamp'] <= r['deadline'] and event['cost'] < best_cost:
                best_cost, best_route = event['cost'], event['route']
        assert r['route'] == best_route and r['cost'] == best_cost
        check(d, r['route'], r['cost'])
        checked += 1
        assert r['feasible_repair_returned'] == bool(r['events'])
        assert r['native_started'] <= r['repair_returned'] <= r['native_started']+r['native_elapsed']
        if arm != 'greedy':
            details = r['repair_details']
            slots = {s['vertex']:s['gaps'] for s in details['slots']}
            absent = sorted(set(visits)-set(retained))
            assert set(slots) == set(absent)
            assert details['slot_variables'] == 3*len(absent)
            slot_rng = random.Random(int(hashlib.sha256((label+':'+arm).encode()).hexdigest(),16))
            for v in absent:
                gaps = list(range(len(retained)))
                if arm == 'mip_nearest':
                    n=d['n'];m=d['matrix']
                    gaps.sort(key=lambda i:(m[retained[i-1]*n+v]+m[v*n+retained[i]]
                                           -m[retained[i-1]*n+retained[i]],i))
                else:
                    slot_rng.shuffle(gaps)
                assert slots[v] == gaps[:3]
            if r['events']:
                assignments = [(visits[v],i) for v,i in details['assignments']]
                assert len({v for v,i in assignments}) == len(assignments)
                assert len({i for v,i in assignments}) == len(assignments)
                assert all(i in slots[v] for v,i in assignments)
                at = {i:v for v,i in assignments}
                expanded=[]
                for i,v in enumerate(retained):
                    if i in at: expanded.append(at[i])
                    expanded.append(v)
                assert r['events'][0]['route'] == expanded
        rows.append(dict(instance=name,seed=seed,sample=sample,arm=arm,
                         cost=r['cost'],feasible_repair=r['feasible_repair_returned'],
                         result_sha256=digest(path)))
    complete=(folder/'completion.json').exists()
    if complete:
        assert len(rows)==324 and len(data)==9
        assert all((n,s,k,a) in keys for n in data for s in (42,43,44)
                   for k in range(4) for a in ('greedy','mip_nearest','mip_random'))
    report=dict(status='PASS' if rows else 'NO_RESULTS',batch_complete=complete,
                runs=len(rows),routes_checked=checked,rows=rows,
                scope='conditional TRAIN scalar, source, destruction, slots and deadline audit')
    (HERE/'joint_insertion_train_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='rows'}))
