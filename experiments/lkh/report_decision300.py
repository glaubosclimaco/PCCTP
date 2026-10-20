"""Decision report from audited summaries; partial batches cannot yield a verdict."""
import hashlib,json
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]

if __name__=='__main__':
    names=['ruin_total300_train_audit.json','ruin_total300_train_summary.json',
           'ruin_total300_train_anytime.json','ruin_total300_train_mining_exposure_bounds.json']
    documents={n:json.loads((HERE/n).read_text()) for n in names}
    audit,summary,curve,mining=[documents[n] for n in names]
    assert audit['status']=='PASS'
    assert len({d['runs'] for d in documents.values()})==1, 'Refresh observers to the same audit snapshot'
    assert curve['audit_sha256']==hashlib.sha256((HERE/names[0]).read_bytes()).hexdigest()
    complete=audit['batch_complete'] and all(d['status']=='COMPLETE' for d in (summary,curve,mining))
    candidates=[m for m,c in summary['comparisons'].items() if c['gate_passed']] if complete else []
    decision='INCOMPLETE' if not complete else ('TRAIN_SCREEN_PASS_WITH_LIMITATIONS' if candidates else 'TRAIN_SCREEN_FAILED')
    lines=['# Avaliação experimental PCCTP: orçamento total de 300 segundos','',
        f"Estado: **{decision}**. {audit['runs']}/81 execuções auditadas; "
        f"{summary['completed_geometries']}/9 geometrias completas.",'',
        'Nenhuma conclusão de publicação, significância ou generalização é automática.',
        'Dados de treino selecionados por dificuldade; teste reservado permanece fechado.','',
        '| Método contra MDM300 | Ganho médio (%) | Vitórias | Perdas | Critério |',
        '|---|---:|---:|---:|---|']
    for mode,c in summary['comparisons'].items():
        effect='—' if c['mean_gain_percent'] is None else f"{c['mean_gain_percent']:.6f}"
        status='pendente' if not complete else ('passou' if c['gate_passed'] else 'não passou')
        lines.append(f"| {mode} | {effect} | {c['wins']} | {c['losses']} | {status} |")
    lines += ['', '## Custos e consistência nas geometrias completas', '',
        '| Instância | Média MDM | Média H10 | Média H60 | H10: sementes melhores/iguais/piores | H60: sementes melhores/iguais/piores |',
        '|---|---:|---:|---:|---|---|']
    lookup={(r['instance'],r['seed'],r['mode']):r for r in audit['rows']}
    for row in summary['rows']:
        name=row['instance'];counts=[]
        for mode in ('hybrid10','hybrid60'):
            differences=[lookup[name,s,'mdm']['cost']-lookup[name,s,mode]['cost'] for s in (70,71,72)]
            counts.append('/'.join(str(sum(test(d) for d in differences)) for test in
                                   (lambda d:d>0,lambda d:d==0,lambda d:d<0)))
        means=row['means']
        lines.append(f"| {name} | {means['mdm']:.3f} | {means['hybrid10']:.3f} | {means['hybrid60']:.3f} | {counts[0]} | {counts[1]} |")
    lines += ['', 'As contagens por semente são descritivas dentro da mesma instância;',
        'não multiplicam o número de geometrias independentes.',
        f"Geometrias com falha de factibilidade: {len(summary['feasibility_failures'])}."]
    lines += ['', 'O desfecho primário é a média das reduções relativas das médias por',
        'geometria. Curvas em outros tempos são descritivas e não substituem300s.',
        'As contagens acima referem-se apenas a geometrias completas; resultados',
        'parciais não fundamentam aprovação ou rejeição do método.','',
        '| Método | Execuções | Mineração observada | Mineração pontual comprovada | Limite superior desconhecido |',
        '|---|---:|---:|---:|---:|']
    for mode,c in sorted(mining['modes'].items()):
        lines.append(f"| {mode} | {c['runs']} | {c['runs_with_observed_call']} | {c['runs_with_proven_timely_call']} | {c['runs_with_unknown_upper_bound']} |")
    lines += ['', 'Mineração observada em export final pode incluir sobretempo cooperativo.',
        'A prova de ativação pontual usa export/checkpoint anterior ao orçamento nativo.',
        'Ativação não demonstra eficácia da mineração. Tempo de CPU do processo pai',
        'não deve ser apresentado como CPU total: LKH e mineradores são subprocessos.',
        'A comparação principal usa tempo de parede debitado e incumbentes pontuais.','',
        'Este orçamento não reproduz as70iterações dos artigos. A decisão sobre',
        'publicação ou abandono exige confronto com a referência completa;',
        'published70_calibration_protocol.md define a calibração inicial, que',
        'também não substitui a avaliação final com dez sementes.', '',
        'Se o filtro de treino passar, a próxima etapa é desenhar confirmação',
        'prospectiva mais ampla, com parâmetros congelados e controles pertinentes.',
        'Se ambos falharem, preservar o diagnóstico e reavaliar essa linha antes',
        'de consumir o teste reservado ou continuar a redação como se houvesse ganho.']
    target=ROOT/'research_notes/DECISION300_REPORT.md';target.write_text('\n'.join(lines)+'\n')
    provenance=dict(status=decision,runs=audit['runs'],candidates=candidates,
        source_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names},
        report_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
        scope='training investment screen, never a publication acceptance criterion')
    (HERE/'decision300_report_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print(json.dumps(dict(status=decision,runs=audit['runs'],report=str(target))))
