"""Outcome-free instance descriptors and elementary visit-count lower bounds."""
import hashlib,json,math
from pathlib import Path
from raw_route_audit import load
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]

if __name__=='__main__':
    protocol_path=HERE/'ruin_total300_awake_train/protocol.json'
    protocol=json.loads(protocol_path.read_text())
    roles=json.loads((ROOT/'experiments/ml_lns/decoder_initialized_states/protocol.json').read_text())['roles']
    rows=[]
    for case in protocol['cases']:
        assert roles[case['geometry']]=='train'
        d=load(case['instance']);n=d['n']
        mandatory={v for v,t in enumerate(d['types']) if t==1}
        optional={v for v,t in enumerate(d['types']) if t==0}
        demands={v for v,t in enumerate(d['types']) if t==2}
        assert len(mandatory|optional|demands)==n and len(mandatory|optional)>=3
        assert all(p>=0 for p in d['prizes']) and d['quota']>=0
        covers={w:{v for v in mandatory|optional if d['matrix'][v*n+w]<=d['radius']} for w in demands}
        assert all(covers.values()),'Uncoverable benchmark demand'
        uncovered={w for w,c in covers.items() if not c&mandatory}
        prize_mandatory=sum(d['prizes'][v] for v in mandatory)
        deficit=max(0.,d['quota']-prize_mandatory)
        optional_prizes=sorted((d['prizes'][v] for v in optional),reverse=True)
        total=0.;prize_lb=0
        while total<deficit and prize_lb<len(optional_prizes):
            total+=optional_prizes[prize_lb];prize_lb+=1
        assert total>=deficit,'All visitable prizes cannot satisfy quota'
        max_cover=max((sum(v in covers[w] for w in uncovered) for v in optional),default=0)
        cover_lb=math.ceil(len(uncovered)/max_cover) if uncovered else 0
        optional_lb=max(prize_lb,cover_lb,max(0,3-len(mandatory)))
        rows.append(dict(instance=case['instance'],geometry=case['geometry'],n=n,
            mandatory=len(mandatory),optional=len(optional),demands=len(demands),
            quota=d['quota'],radius=d['radius'],visitable_prize=prize_mandatory+sum(optional_prizes),
            mandatory_prize=prize_mandatory,prize_deficit_after_mandatory=deficit,
            demands_not_covered_by_mandatory=len(uncovered),
            optional_count_prize_lower_bound=prize_lb,optional_count_coverage_lower_bound=cover_lb,
            selected_count_lower_bound=len(mandatory)+optional_lb,
            input_sha256=d['input_hash']))
    report=dict(scope='nine selected hardTRAIN instances only; no outcomes, solver, metric assumption or new-method claim',
        protocol_sha256=hashlib.sha256(protocol_path.read_bytes()).hexdigest(),rows=rows,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (HERE/'hard_train_instance_profile.json').write_text(json.dumps(report,indent=2)+'\n')
    lines=['# Estrutura dos nove casos difíceis de treino','',
        'Descritores calculados exclusivamente das entradas; nenhum resultado de busca entra nos cálculos. '
        'Não representam os107 casos do artigo nem incluem validação ou teste.','',
        '| Instância | Obrigatórios | Opcionais | Demandas | Déficit de prêmio após obrigatórios | Demandas descobertas | LB opcionais por prêmio | LB opcionais por cobertura | LB visitados |',
        '|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    fields=['mandatory','optional','demands','prize_deficit_after_mandatory','demands_not_covered_by_mandatory',
        'optional_count_prize_lower_bound','optional_count_coverage_lower_bound','selected_count_lower_bound']
    for r in rows:lines.append('| '+r['instance']+' | '+' | '.join(str(r[k]) for k in fields)+' |')
    lines.extend(['','O limite por prêmio soma os maiores prêmios opcionais até cobrir o déficit. '
        'Nenhum conjunto com menos opcionais alcança esse total. O limite por cobertura é '
        'o teto do número de demandas ainda descobertas dividido pela maior quantidade dessas '
        'demandas que um opcional pode cobrir. A união de k opcionais cobre no máximo k vezes '
        'essa quantidade. Combina-se o máximo desses dois limites e da necessidade de pelo '
        'menos três vértices no ciclo; não se somam os limites, pois um vértice pode cumprir '
        'ambas as funções. São limites elementares de cardinalidade, sem pretensão de novidade '
        'ou de ótimo para o custo da rota.',''])
    (ROOT/'research_notes/hard_train_instance_profile.md').write_text('\n'.join(lines))
    print(json.dumps(dict(instances=len(rows),scope=report['scope'])))
