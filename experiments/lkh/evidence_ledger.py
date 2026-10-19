"""Trace exploratory outcomes without conflating different baselines or metrics."""
import json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
SPECS=[
('total_time_audit.json','MDM9s + LKH versus MDM10s','mean_geometry_gain_percent','% cost reduction, mean by geometry','gate_passed'),
('validation_audit.json','Frozen 9+1 validation','mean_geometry_gain_percent','% cost reduction, mean by geometry','gate_passed'),
('horizon30_audit.json','MDM29s + LKH versus MDM30s','mean_geometry_gain_percent','% cost reduction, TRAIN persistence; no test-unlocking gate',None),
('support_audit.json','Support + LKH from fixed incumbent','wins_vs_fixed','geometries beating fixed-support LKH','gate_passed'),
('support_total_audit.json','Support + LKH, 15s total',None,'feasibility failures retained; conditional cost comparison','gate_passed'),
('support_common_fallback_audit.json','Support + LKH, common fallback',None,'15s total; implementation confound in Python 2-opt control','gate_passed'),
('direct_support_summary.json','Cheap construction + support descent','mean_gain_percent','% relative to historical cached MDM; not paired','gate_passed'),
('diversification_audit.json','Diversified support descent','mean_gain_percent','% relative to historical cached MDM; not paired','gate_passed'),
('candidate_oracle_audit.json','Candidate-selection offline oracle','mean_gain_percent','% over best top8, one state per TRAIN geometry','gate_passed'),
('hard_fragment_audit.json','Contracted exact-neighborhood offline oracle','mean_gain_percent','% over fixed LKH incumbent; extra offline work','gate_passed'),
('fragment_selection_results.json','Ridge selection, held-out TRAIN geometry','mean_delta_vs_random_pp','percentage-point advantage over random expectation','gate_passed'),
('exact_random_total_audit.json','Random exact repair, 15s total',None,'major return-time anomalies and one technical retry; exploratory','gate_passed')]
if __name__=='__main__':
    rows=[]
    for filename,title,metric,meaning,gate in SPECS:
        p=HERE/filename
        if not p.exists():continue
        r=json.loads(p.read_text());rows.append(dict(title=title,artifact=str(p.relative_to(ROOT)),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),metric_name=metric,metric_value=r.get(metric) if metric else None,metric_interpretation=meaning,gate=r.get(gate) if gate else None))
    report=dict(scope='Exploratory decision ledger, not a common-baseline ranking or confirmatory meta-analysis',rows=rows,publication_ready=False,reserved_test_opened=False,invalid_instrumentation=['experiments/lkh/timing_diagnostic/invalid_clock_protocol.json'])
    (HERE/'evidence_ledger.json').write_text(json.dumps(report,indent=2)+'\n')
    lines=['# Registro de evidência e decisões exploratórias','','As métricas abaixo usam referências, custos de coleta e escopos diferentes. Não comparar seus números como um ranking único. Passar em um critério de triagem não demonstra publicabilidade, novidade ou significância confirmatória. As hipóteses foram investigadas sequencialmente; essa seleção impede apresentar seus intervalos exploratórios como confirmação independente.','','| Investigação | Métrica registrada | Critério de triagem | Artefato |','|---|---|---|---|']
    for r in rows:
        value=f"{r['metric_value']:.6g}; " if isinstance(r['metric_value'],(int,float)) else ''
        decision='Passou' if r['gate'] is True else 'Falhou' if r['gate'] is False else 'Sem liberação de teste'
        lines.append(f"| {r['title']} | {value}{r['metric_interpretation']} | {decision} | `{r['artifact']}` |")
    lines+=['','O diagnóstico inicial de CPU/relógio foi invalidado por origem temporal incompatível entre Python e C++; sua repetição corrigida não explica retrospectivamente os atrasos anteriores. Os comparadores principais já passavam CLOCK_MONOTONIC ao executável nativo. Não converter resultados limitados por atrasos em demonstrações controladas de velocidade.','','O diagnóstico de 70 iterações está em andamento. Casos censurados serão mantidos e o LKH recebe tempo adicional explicitamente separado. O teste reservado continua fechado; o artigo permanece incompleto. Resultados de uma execução parcial não substituem o conjunto completo.']
    (ROOT/'research_notes/EVIDENCE_LEDGER.md').write_text('\n'.join(lines)+'\n')
    print(f'Recorded {len(rows)} exploratory studies')
