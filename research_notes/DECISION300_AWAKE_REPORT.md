> Auditoria ambiental independente: 27/27 blocos validados; elegibilidade ambiental do lote completo: True. Execuções no bloco ainda aberto aguardam seu diagnóstico ambiental.

# Avaliação experimental PCCTP: orçamento total de 300 segundos

Estado: **TRAIN_SCREEN_PASS_WITH_LIMITATIONS**. 81/81 execuções auditadas; 9/9 geometrias completas.

Nenhuma conclusão de publicação, significância ou generalização é automática.
Dados de treino selecionados por dificuldade; teste reservado permanece fechado.

| Método contra MDM300 | Ganho médio (%) | Vitórias | Perdas | Critério |
|---|---:|---:|---:|---|
| hybrid10 | -1.279343 | 6 | 3 | não passou |
| hybrid60 | 0.979441 | 8 | 1 | passou |

## Custos e consistência nas geometrias completas

| Instância | Média MDM | Média H10 | Média H60 | H10: sementes melhores/iguais/piores | H60: sementes melhores/iguais/piores |
|---|---:|---:|---:|---|---|
| g_a280_VT56_T168_W56_50 | 2029.333 | 1998.000 | 2001.333 | 3/0/0 | 3/0/0 |
| g_gil262_VT51_T158_W53_75 | 2016.333 | 1982.667 | 1985.000 | 3/0/0 | 3/0/0 |
| g_gr229_VT45_T138_W46_50 | 98281.333 | 96722.333 | 96737.667 | 3/0/0 | 3/0/0 |
| g_kroA200_VT40_T120_W40_25 | 23497.667 | 23505.333 | 23495.000 | 1/1/1 | 2/1/0 |
| g_pr226_VT44_T136_W46_50 | 68110.000 | 67286.000 | 67286.000 | 3/0/0 | 3/0/0 |
| g_rd400_VT80_T240_W80_25 | 10773.667 | 12872.333 | 10670.333 | 1/0/2 | 2/0/1 |
| g_si175_VT35_T105_W35_50 | 7045.333 | 6908.000 | 6898.000 | 3/0/0 | 3/0/0 |
| g_ts225_VT45_T135_W45_50 | 108998.000 | 108866.333 | 108866.333 | 1/2/0 | 1/2/0 |
| g_tsp225_VT45_T45_W135_75 | 2447.000 | 2449.000 | 2449.000 | 0/0/3 | 0/0/3 |

As contagens por semente são descritivas dentro da mesma instância;
não multiplicam o número de geometrias independentes.
Geometrias com falha de factibilidade: 0.

O desfecho primário é a média das reduções relativas das médias por
geometria. Curvas em outros tempos são descritivas e não substituem300s.
As contagens acima referem-se apenas a geometrias completas; resultados
parciais não fundamentam aprovação ou rejeição do método.

| Método | Execuções | Mineração observada | Mineração pontual comprovada | Limite superior desconhecido |
|---|---:|---:|---:|---:|
| hybrid10 | 27 | 0 | 0 | 0 |
| hybrid60 | 27 | 2 | 1 | 0 |
| mdm | 27 | 11 | 9 | 0 |

Mineração observada em export final pode incluir sobretempo cooperativo.
A prova de ativação pontual usa export/checkpoint anterior ao orçamento nativo.
Ativação não demonstra eficácia da mineração. Tempo de CPU do processo pai
não deve ser apresentado como CPU total: LKH e mineradores são subprocessos.
A comparação principal usa tempo de parede debitado e incumbentes pontuais.

Este orçamento não reproduz as70iterações dos artigos. A decisão sobre
publicação ou abandono exige confronto com a referência completa;
published70_calibration_protocol.md define a calibração inicial, que
também não substitui a avaliação final com dez sementes.

Se o filtro de treino passar, a próxima etapa é desenhar confirmação
prospectiva mais ampla, com parâmetros congelados e controles pertinentes.
Se ambos falharem, preservar o diagnóstico e reavaliar essa linha antes
de consumir o teste reservado ou continuar a redação como se houvesse ganho.
