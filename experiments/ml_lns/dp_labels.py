"""Fresh labels for fragment_clique/DP, restricted to training geometries.

Old route_context MIP targets are not reused. Size-limited queries stay in the
dataset. Collector uses an exclusive lock and atomic file publication.
"""
import hashlib
import json
import os
from pathlib import Path
import random
import time
from dp_repair import DPRepair
from exchange_signature import signature,NAMES
from pilot import candidates,features,read_instance,validate

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]

def footprint(route,released):
    flags=[v in released for v in route]
    if not any(flags):return len(released)+1
    runs=sum(flags[i] and not flags[i-1] for i in range(len(flags)))
    return len(released)+runs

def main():
    out=HERE/'dp_training_labels';out.mkdir(exist_ok=True)
    lock=out/'.collector.lock'
    fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY);os.close(fd)
    try:
        seen=set()
        for path in sorted((HERE/'decoder_initialized_states').glob('g_*.json')):
            baseline=json.loads(path.read_text());geometry=baseline['geometry']
            if baseline['role']!='train' or geometry in seen:continue
            snapshot=next((s for s in baseline['snapshots'] if s['iteration']==1),None)
            if snapshot is None:continue
            seen.add(geometry);target=out/(path.stem+'.phase1.json')
            provenance={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in
                        ['dp_labels.py','dp_repair.py','path_dp.cpp','pilot.py','exchange_signature.py']}
            source_hash=hashlib.sha256(path.read_bytes()).hexdigest()
            if target.exists():
                previous=json.loads(target.read_text())
                assert previous['source_sha256']==source_hash and previous['code_hashes']==provenance
                continue
            d=read_instance(ROOT/'data'/baseline['instance'])
            mapping={int(v):i for i,v in enumerate(d['visits'])};route=[mapping[v] for v in snapshot['route']]
            cost=validate(d,route);model=DPRepair(d);records=[]
            rng=random.Random(42)
            for kind,u in candidates(d,route,8,12,42):
                began=time.perf_counter()
                base=features(d,route,u);joint=signature(d,route,u)
                shared_seconds=time.perf_counter()-began
                near=sorted(d['mandatory'],key=lambda v:(min(d['costs'][v,w] for w in u),v))[:4]
                random_anchors=rng.sample(sorted(d['mandatory']),min(4,len(d['mandatory'])))
                for anchor_kind,anchors in [('near',near),('random',random_anchors)]:
                    began=time.perf_counter();j=u|set(anchors);q=footprint(route,j)
                    # Anchor deletion saving is a routing proxy only; anchors remain visited.
                    anchor_saving=sum(d['costs'][route[route.index(v)-1],v]+
                                      d['costs'][v,route[(route.index(v)+1)%len(route)]]-
                                      d['costs'][route[route.index(v)-1],route[(route.index(v)+1)%len(route)]]
                                      for v in anchors)/max(1.,cost)
                    proximity=sum(min(d['costs'][v,w] for w in u) for v in anchors)/max(1.,cost)
                    x=base+[q/18.,len(anchors)/4.,anchor_saving,proximity]+joint
                    feature_seconds=shared_seconds+time.perf_counter()-began
                    result=model.solve(route,j,.2,'fragment_clique',seed=42)
                    assert result['contracted_nodes']==q
                    records.append(dict(kind=kind,anchor_kind=anchor_kind,optional_released=sorted(u),mandatory_released=anchors,
                                        features=x,feature_seconds=feature_seconds,
                                        normalized_gain=(cost-result['cost'])/max(1.,cost),**result))
            row=dict(instance=baseline['instance'],geometry=geometry,role='train',phase=1,seed=baseline['seed'],
                     source_sha256=source_hash,code_hashes=provenance,initial_cost=cost,initial_route=snapshot['route'],
                     feature_names=['selected_fraction','prize_fraction','residual_prize','affected_coverage',
                                    'repairable_coverage','coverage_overlap','deletion_saving','distance_spread',
                                    'component_fraction','anchor_fraction','anchor_deletion_proxy','anchor_proximity']+NAMES,
                     backend='DP',operator='fragment_clique',query_seconds=.2,component_cap=18,records=records)
            temporary=target.with_suffix('.partial');temporary.write_text(json.dumps(row,indent=2)+'\n');temporary.replace(target)
            print(json.dumps(dict(instance=row['instance'],queries=len(records),
                                  improving=sum(r['normalized_gain']>0 for r in records),
                                  size_limits=sum(r['reason']=='component_limit' for r in records))),flush=True)
    finally:lock.unlink(missing_ok=True)

if __name__=='__main__':main()
