# Registro de evidência e decisões exploratórias

As métricas abaixo usam referências, custos de coleta e escopos diferentes. Não comparar seus números como um ranking único. Passar em um critério de triagem não demonstra publicabilidade, novidade ou significância confirmatória. As hipóteses foram investigadas sequencialmente; essa seleção impede apresentar seus intervalos exploratórios como confirmação independente.

| Investigação | Métrica registrada | Critério de triagem | Artefato |
|---|---|---|---|
| MDM9s + LKH versus MDM10s | 1.37394; % cost reduction, mean by geometry | Passou | `experiments/lkh/total_time_audit.json` |
| Frozen 9+1 validation | 0.89319; % cost reduction, mean by geometry | Falhou | `experiments/lkh/validation_audit.json` |
| MDM29s + LKH versus MDM30s | 0.479637; % cost reduction, TRAIN persistence; no test-unlocking gate | Sem liberação de teste | `experiments/lkh/horizon30_audit.json` |
| Support + LKH from fixed incumbent | 5; geometries beating fixed-support LKH | Passou | `experiments/lkh/support_audit.json` |
| Support + LKH, 15s total | feasibility failures retained; conditional cost comparison | Falhou | `experiments/lkh/support_total_audit.json` |
| Support + LKH, common fallback | 15s total; implementation confound in Python 2-opt control | Falhou | `experiments/lkh/support_common_fallback_audit.json` |
| Cheap construction + support descent | -2.57287; % relative to historical cached MDM; not paired | Falhou | `experiments/lkh/direct_support_summary.json` |
| Diversified support descent | -6.13308; % relative to historical cached MDM; not paired | Falhou | `experiments/lkh/diversification_audit.json` |
| Candidate-selection offline oracle | 0.158958; % over best top8, one state per TRAIN geometry | Falhou | `experiments/lkh/candidate_oracle_audit.json` |
| Contracted exact-neighborhood offline oracle | 0.500945; % over fixed LKH incumbent; extra offline work | Passou | `experiments/lkh/hard_fragment_audit.json` |
| Ridge selection, held-out TRAIN geometry | 0.0986079; percentage-point advantage over random expectation | Falhou | `experiments/lkh/fragment_selection_results.json` |
| Random exact repair, 15s total | major return-time anomalies and one technical retry; exploratory | Falhou | `experiments/lkh/exact_random_total_audit.json` |

O diagnóstico inicial de CPU/relógio foi invalidado por origem temporal incompatível entre Python e C++; sua repetição corrigida não explica retrospectivamente os atrasos anteriores. Os comparadores principais já passavam CLOCK_MONOTONIC ao executável nativo. Não converter resultados limitados por atrasos em demonstrações controladas de velocidade.

O diagnóstico de 70 iterações está em andamento. Casos censurados serão mantidos e o LKH recebe tempo adicional explicitamente separado. O teste reservado continua fechado; o artigo permanece incompleto. Resultados de uma execução parcial não substituem o conjunto completo.
