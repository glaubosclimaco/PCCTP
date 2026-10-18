"""Independent route/deadline audit and geometry-level exploratory analysis."""
import json,hashlib,random,statistics
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from experiments.ml_lns.pilot import read_instance,validate
HERE=Path(__file__).resolve().parent;OUT=HERE/'total_time'


def audit(validation=False):
    global OUT
    if validation: OUT=HERE/'validation'
    assert (OUT/'completion.json').exists(),'Batch still running'
    cases=json.loads((ROOT/'experiments/original_mdm/hard_training_baseline/summary.json').read_text())['rows']
    if validation:
        cases=json.loads((OUT/'protocol.json').read_text())['cases']
    n=len(cases)
    assert n==(7 if validation else 9)
    geometries=[];pair_rows=[];hashes={};failures=0
    for case in cases:
        name=case['instance'];d=read_instance(ROOT/'data'/name)
        mapping={int(v):i for i,v in enumerate(d['visits'])}
        gains=[];costs={'mdm':[],'hybrid':[]}
        for seed in (42,43,44):
            pair={}
            for mode in ('mdm','hybrid'):
                path=OUT/f'{name}.seed{seed}.{mode}.json'
                native_path=OUT/f'{name}.seed{seed}.{mode}.native.json'
                r=json.loads(path.read_text());native=json.loads(native_path.read_text())
                hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
                assert r['instance']==name and r['seed']==seed and r['mode']==mode
                assert native['input_sha256']==hashlib.sha256((ROOT/'data'/name).read_bytes()).hexdigest()
                budget=10. if mode=='mdm' else 9.
                assert native['budget_seconds']==budget and native['elite_count_bounded']
                assert native['clock_source']=='clock_gettime(CLOCK_MONOTONIC)'
                assert native['feasible_before_deadline']==bool(native['anytime'])==r['feasible']
                previous=float('inf');last_time=-1.
                for event in native['anytime']:
                    assert last_time<=event['seconds']<=budget+1e-9 and event['cost']<previous
                    assert validate(d,[mapping[v] for v in event['route']])==event['cost']
                    previous=event['cost'];last_time=event['seconds']
                if not r['feasible']:
                    failures+=1;pair[mode]=r;continue
                assert validate(d,[mapping[v] for v in r['route']])==r['cost']
                assert r['cost']<=native['cost']
                if mode=='mdm':
                    assert r['cost']==native['cost'] and r['route']==native['route']
                if r['lkh_used']:
                    lkh=r['lkh']
                    assert lkh['initial_cost']==native['cost']
                    assert set(lkh['route'])==set(native['route'])
                    assert validate(d,[mapping[v] for v in lkh['route']])==lkh['cost']
                if r['cost']<native['cost']:
                    assert mode=='hybrid' and r['lkh_used'] and r['wall_seconds']<=10.
                    assert r['cost']==r['lkh']['cost']
                costs[mode].append(r['cost']);pair[mode]=r
            if all(pair[m]['feasible'] for m in pair):
                gain=100*(pair['mdm']['cost']-pair['hybrid']['cost'])/pair['mdm']['cost']
                gains.append(gain)
                pair_rows.append(dict(instance=name,seed=seed,gain_percent=gain,
                                      mdm_cost=pair['mdm']['cost'],hybrid_cost=pair['hybrid']['cost']))
        assert len(gains)==3,'Feasibility failure requires explicit censored analysis'
        avg=statistics.mean(gains)
        mean_difference=statistics.mean(costs['mdm'])-statistics.mean(costs['hybrid'])
        geometries.append(dict(instance=name,geometry=case['geometry'],seeds=3,
                               mean_paired_gain_percent=avg,median_paired_gain_percent=statistics.median(gains),
                               mdm_mean=statistics.mean(costs['mdm']),hybrid_mean=statistics.mean(costs['hybrid']),
                               classification='win' if mean_difference>1e-9 else 'loss' if mean_difference < -1e-9 else 'tie'))
    means=[r['mean_paired_gain_percent'] for r in geometries]
    wins=sum(r['classification']=='win' for r in geometries)
    losses=sum(r['classification']=='loss' for r in geometries)
    # Resample geometries, retaining their aggregated paired-seed result.
    # Exploratory only: method/cohort were selected during training.
    rng=random.Random(20261018)
    boot=sorted(statistics.mean(rng.choices(means,k=n)) for _ in range(10000))
    average=statistics.mean(means)
    gate=average>0 and wins>=4 and losses<=1 and failures==0 if validation else average>=.5 and wins>=6 and losses<=1 and failures==0
    report=dict(status='PASS',scope='frozen reserved validation' if validation else 'exploratory hard-TRAIN pilot; not confirmatory',
        geometries=n,paired_seeds=3,pairs=3*n,mean_geometry_gain_percent=average,
        median_geometry_gain_percent=statistics.median(means),
        geometry_bootstrap_95_percentile=[boot[249],boot[9749]],
        geometry_wins=wins,geometry_losses=losses,geometry_ties=n-wins-losses,
        feasibility_failures=failures,gate_passed=gate,rows=geometries,pair_rows=pair_rows,
        result_hashes=hashes)
    (HERE/('validation_audit.json' if validation else 'total_time_audit.json')).write_text(json.dumps(report,indent=2)+'\n')
    lines=['# Piloto MDM10s versus MDM9s+LKH','',
           f'{n} geometrias, três sementes pareadas; '+('validação reservada com método congelado.' if validation else 'treino exploratório, não confirmatório.')+' Ambos recebem dez segundos totais. Overruns são registrados; melhorias tardias do LKH são descartadas.',
           '', '| Geometria | Média MDM | Média híbrido | Ganho pareado médio % |','|---|---:|---:|---:|']
    for row in geometries:
        lines.append(f"| {row['geometry']} | {row['mdm_mean']:.2f} | {row['hybrid_mean']:.2f} | {row['mean_paired_gain_percent']:.3f} |")
    lines+=['',f'Média entre geometrias: {average:.3f}%; mediana: {statistics.median(means):.3f}%.',
            f'Vitórias/empates/derrotas por geometria: {wins}/{n-wins-losses}/{losses}. Gate pré-fixado: {gate}.',
            f'Intervalo bootstrap exploratório por geometria: [{boot[249]:.3f}, {boot[9749]:.3f}]%.',
            'O intervalo não corrige a seleção de método/coorte durante a pesquisa. Não alegar confirmação estatística ou publicabilidade só com este piloto.']
    (ROOT/'research_notes'/('lkh_validation_results.md' if validation else 'lkh_total_time_results.md')).write_text('\n'.join(lines)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('rows','pair_rows','result_hashes')}))


if __name__=='__main__':
    assert sys.argv[1:] in ([],['--validation'])
    audit(validation='--validation' in sys.argv)
