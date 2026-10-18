"""Accounting decomposition, not a causal attribution of runtime effects."""
import json,hashlib,statistics,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent

def analyze(folder):
    assert (folder/'completion.json').exists()
    groups={};hashes={}
    for path in sorted(folder.glob('*.hybrid.json')):
        h=json.loads(path.read_text());name=h['instance'];seed=h['seed']
        paths=[path,folder/f'{name}.seed{seed}.mdm.json',folder/f'{name}.seed{seed}.hybrid.native.json']
        _,m,n=[json.loads(p.read_text()) for p in paths]
        for p in paths:hashes[p.name]=hashlib.sha256(p.read_bytes()).hexdigest()
        assert h['feasible'] and m['feasible'] and n['feasible_before_deadline']
        ordering=n['cost']-h['cost'];difference=n['cost']-m['cost'];net=m['cost']-h['cost']
        assert ordering>=0 and net==ordering-difference
        row=dict(instance=name,seed=seed,selected_visits=len(h['route']),baseline_cost=m['cost'],hybrid_native_cost=n['cost'],hybrid_final_cost=h['cost'],accepted_ordering_gain=ordering,native_cost_difference=difference,net_gain=net,ordering_gain_percent=100*ordering/m['cost'],native_difference_percent=100*difference/m['cost'],net_gain_percent=100*net/m['cost'],lkh_called=h['lkh_used'])
        groups.setdefault(name,[]).append(row)
    rows=[]
    for name,rs in groups.items():
        assert len(rs)==3
        rows.append(dict(instance=name,pairs=rs,mean_visits=statistics.mean(r['selected_visits'] for r in rs),mean_ordering_gain_percent=statistics.mean(r['ordering_gain_percent'] for r in rs),mean_native_difference_percent=statistics.mean(r['native_difference_percent'] for r in rs),mean_net_gain_percent=statistics.mean(r['net_gain_percent'] for r in rs)))
    result=dict(scope='descriptive accounting; native difference is not a causal estimate of time sacrifice',batch=folder.name,geometries=len(rows),rows=rows,input_hashes=hashes)
    for key in ('ordering_gain_percent','native_difference_percent','net_gain_percent'):
        result['geometry_mean_'+key]=statistics.mean(r['mean_'+key] for r in rows)
    (HERE/f'{folder.name}_decomposition.json').write_text(json.dumps(result,indent=2)+'\n')
    text=['# Decomposição contábil do ganho: '+folder.name,'','Ganho líquido = melhoria de ordenação aceita − diferença entre o custo nativo do híbrido e o custo do braço MDM. A identidade é exata por par. A diferença nativa combina tempo, trajetória e variação de execução; não estima causalmente o custo de retirar um segundo. Cada geometria tem peso igual.','','| Instância | Visitas médias | Ordenação % | Diferença nativa % | Ganho líquido % |','|---|---:|---:|---:|---:|']
    for r in rows:text.append(f"| {r['instance']} | {r['mean_visits']:.1f} | {r['mean_ordering_gain_percent']:.3f} | {r['mean_native_difference_percent']:.3f} | {r['mean_net_gain_percent']:.3f} |")
    text+=['', 'Uma melhoria de ordenação mantém exatamente os vértices visitados, preservando cobertura e prêmio. O algoritmo não fornece certificado de otimalidade. A decomposição explica os resultados observados, sem estabelecer novidade ou desempenho fora da amostra.']
    (HERE.parents[1]/'research_notes'/f'lkh_{folder.name}_decomposition.md').write_text('\n'.join(text)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('rows','input_hashes')}))

if __name__=='__main__':
    assert len(sys.argv)==2 and sys.argv[1] in ('total_time','validation','horizon30')
    analyze(HERE/sys.argv[1])
