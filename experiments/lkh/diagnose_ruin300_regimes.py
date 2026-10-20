"""Descriptive time/support witnesses on complete TRAIN geometries only."""
import hashlib,json
from pathlib import Path
from raw_route_audit import load,check
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

if __name__=='__main__':
    audit_path=HERE/'ruin_total300_train_audit.json'
    curve_path=HERE/'ruin_total300_train_anytime.json'
    summary_path=HERE/'ruin_total300_train_summary.json'
    audit=json.loads(audit_path.read_text());curve=json.loads(curve_path.read_text())
    summary=json.loads(summary_path.read_text())
    assert audit['status']=='PASS' and curve['audit_sha256']==digest(audit_path)
    assert audit['runs']==curve['runs']==summary['runs']
    complete={r['instance'] for r in summary['rows']}
    assert len(complete)==summary['completed_geometries']
    rows=[]
    for item in audit['rows']:
        if item['instance'] not in complete or item['mode']=='mdm':continue
        path=HERE/'ruin_total300_train'/f"{item['instance']}.seed{item['seed']}.{item['mode']}.json"
        assert digest(path)==item['result_sha256']
        row=json.loads(path.read_text());data=load(row['instance']);n=data['n']
        check(data,row['pre_route'],row['pre_cost']);check(data,row['route'],row['cost'])
        initial=set(row['pre_route']);final=set(row['route'])
        removed=initial-final;added=final-initial;retained=initial-removed
        mandatory={v for v,t in enumerate(data['types']) if t==1}
        assert mandatory<=retained
        uncovered=[w for w,t in enumerate(data['types']) if t==2 and not any(
            data['matrix'][v*n+w]<=data['radius'] for v in retained)]
        retained_prize=sum(data['prizes'][v] for v in retained)
        deadline=row['native_started']+300.
        timely=[e for e in row['post']['trace'] if e['native_timestamp']<=deadline] if row['post'] else []
        best_events=[e for e in timely if e['cost']==row['cost']]
        final_time=min(e['native_timestamp']-row['native_started'] for e in best_events) if best_events else None
        rows.append(dict(instance=row['instance'],seed=row['seed'],mode=row['mode'],
            initialization_cost=row['pre_cost'],final_cost=row['cost'],
            post_initialization_gain_percent=100*(row['pre_cost']-row['cost'])/row['pre_cost'],
            removed=sorted(removed),added=sorted(added),support_changed=bool(removed or added),
            net_removed_without_added_uncovered_vertices=uncovered,
            net_removed_without_added_prize_deficit=max(0.,data['quota']-retained_prize),
            final_incumbent_post_event_seconds=final_time,
            source_sha256=digest(path),input_sha256=data['input_hash']))
    assert len(rows)==len(complete)*6
    modes={}
    for mode in ('hybrid10','hybrid60'):
        selected=[r for r in rows if r['mode']==mode]
        modes[mode]=dict(runs=len(selected),
            support_changed_runs=sum(r['support_changed'] for r in selected),
            net_coverage_repair_witnesses=sum(bool(r['net_removed_without_added_uncovered_vertices']) for r in selected),
            net_prize_repair_witnesses=sum(r['net_removed_without_added_prize_deficit']>0 for r in selected),
            improved_after_initialization=sum(r['final_cost']<r['initialization_cost'] for r in selected),
            final_incumbent_first_seen_after120=sum(r['final_incumbent_post_event_seconds'] is not None
                and r['final_incumbent_post_event_seconds']>120 for r in selected))
    report=dict(status=summary['status'],runs=summary['runs'],completed_geometries=len(complete),
        modes=modes,rows=rows,curve=curve['curve'],
        audit_sha256=digest(audit_path),curve_sha256=digest(curve_path),summary_sha256=digest(summary_path),
        scope='descriptive complete TRAIN geometries only; net-support witnesses do not reconstruct destruction paths or establish causality; no new selection gate')
    (HERE/'ruin_total300_train_regimes.json').write_text(json.dumps(report,indent=2)+'\n')
    lines=['# Diagnóstico descritivo do lote de 300 segundos','',
        f"Estado: {summary['status']}; {summary['runs']}/81 execuções; {len(complete)}/9 geometrias completas.",'',
        'Somente geometrias com as nove execuções entram nas comparações abaixo. '
        'Estas observações de treino não substituem o MDM completo de 70 iterações nem alteram o critério final de 300 segundos.','',
        '| Horizonte (s) | H10: ganho médio (%) | H60: ganho médio (%) | Geometrias |',
        '|---:|---:|---:|---:|']
    for horizon in curve['grid']:
        entries={r['mode']:r for r in curve['curve'] if r['horizon']==horizon}
        values=[entries[m]['mean_gain_percent'] for m in ('hybrid10','hybrid60')]
        formatted=['indisponível' if v is None else f'{v:.4f}' for v in values]
        lines.append(f"| {horizon} | {formatted[0]} | {formatted[1]} | {entries['hybrid10']['complete_geometry_pairs']} |")
    lines.extend(['','Valores positivos indicam menor custo que o MDM no mesmo horizonte. '
        'Médias são calculadas por geometria antes da agregação.','',
        '| Variante | Execuções | Suporte alterado | Testemunho de cobertura | Testemunho de prêmio | Melhorou após inicialização | Melhor incumbente após 120 s |',
        '|---|---:|---:|---:|---:|---:|---:|'])
    for mode,r in modes.items():
        lines.append('| '+mode+' | '+' | '.join(str(r[k]) for k in (
            'runs','support_changed_runs','net_coverage_repair_witnesses','net_prize_repair_witnesses',
            'improved_after_initialization','final_incumbent_first_seen_after120'))+' |')
    lines.extend(['','Os testemunhos verificam se remover do suporte inicial todos os vértices '
        'ausentes na solução final, sem inserir os novos, violaria cobertura ou prêmio. '
        'Não representam a trajetória efetiva de destruição/reparação e não provam a causa do ganho. '
        'As colunas de cobertura e prêmio podem se sobrepor. O horário do melhor incumbente '
        'refere-se ao primeiro evento registrado com o custo final durante a fase posterior à inicialização.',''])
    (ROOT/'research_notes/ruin300_regimes.md').write_text('\n'.join(lines))
    print(json.dumps({k:v for k,v in report.items() if k not in ('rows','curve')}))
