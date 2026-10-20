# Diagnóstico passivo de inicialização no lote300s

Snapshot: 69/81 execuções auditadas; 7 geometrias completas. Lote completo: False.

Nenhuma alteração de política ou lançamento de otimizador. Todas as geometrias completas são mantidas.

| Instância | Ganho H10 (%) | Ganho H60 (%) | Sem solução MDM no prazo H10/H60 |
|---|---:|---:|---:|
| g_a280_VT56_T168_W56_50 | 1.544021 | 1.379763 | 0/3; 0/3 |
| g_gil262_VT51_T158_W53_75 | 1.669697 | 1.553976 | 0/3; 0/3 |
| g_gr229_VT45_T138_W46_50 | 1.586263 | 1.570661 | 0/3; 0/3 |
| g_kroA200_VT40_T120_W40_25 | -0.032627 | 0.011349 | 0/3; 0/3 |
| g_pr226_VT44_T136_W46_50 | 1.209808 | 1.209808 | 1/3; 0/3 |
| g_rd400_VT80_T240_W80_25 | -19.479595 | 0.959129 | 2/3; 0/3 |
| g_si175_VT35_T105_W35_50 | 1.949281 | 2.091219 | 0/3; 0/3 |

## Detalhe de rd400

| Semente | Método | MDM tempestivo | Inicial | Final | Suporte inicial | Chamadas/exit0/timeout/other |
|---|---|---|---:|---:|---:|---|
| 70 | hybrid10 | False | 15499 | 13975 | 306 | 1065/0/1032/33 |
| 70 | hybrid60 | True | 10919 | 10630 | 257 | 1081/852/199/30 |
| 71 | hybrid10 | True | 11453 | 10667 | 262 | 1482/1037/445/0 |
| 71 | hybrid60 | True | 11022 | 10631 | 262 | 1002/596/383/23 |
| 72 | hybrid10 | False | 14143 | 13975 | 262 | 1187/0/1186/1 |
| 72 | hybrid60 | True | 11486 | 10750 | 259 | 985/0/984/1 |

## Interpretação limitada

A inicialização de10s falha em fornecer um incumbente MDM tempestivo em duas sementes de rd400. A reparação então começa do fallback. Nessas duas fases, o contador exit0 é zero. Mesmo H60/seed72 tem zero exit0 e ainda reduz o custo: reparo guloso pode produzir melhorias sem reordenação LKH bem-sucedida. Isso não identifica uma causa única nem prova benefício do cache.

A perda de H10 em rd400 é preservada; não excluir a geometria, trocar sementes ou retunar a fase. O protocolo de300s permanece o mesmo. Comparação MDM70 completa continua necessária.

Contraponto: pr226/seed70/H10 também parte do fallback, mas cai de 88536 para 67286, com 3189 registros exit0. A ausência de MDM tempestivo, isoladamente, não explica todas as diferenças de qualidade.

- Exit code zero does not prove a timely validated improving tour.
- Counter buckets can include the final call after the outer event deadline.
- Other status combines nonzero exit and exhausted setup allowance; no per-call logs separate them.
- Poor initialization and ordering statuses are associated observations, not isolated causal effects.
