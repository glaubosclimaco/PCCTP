# MDM30s versus MDM29s+LKH

9 geometrias, três sementes pareadas; treino exploratório, não confirmatório. Ambos recebem 30 segundos totais. Overruns são registrados; melhorias tardias do LKH são descartadas.

| Geometria | Média MDM | Média híbrido | Ganho pareado médio % |
|---|---:|---:|---:|
| g_a280 | 2124.67 | 2131.67 | -0.359 |
| g_gil262 | 2068.00 | 2056.33 | 0.570 |
| g_gr229 | 99658.67 | 99115.33 | 0.541 |
| g_kroA200 | 23516.33 | 23497.67 | 0.079 |
| g_pr226 | 68340.33 | 68223.00 | 0.172 |
| g_rd400 | 11374.00 | 11087.67 | 2.474 |
| g_si175 | 7079.33 | 7038.00 | 0.582 |
| g_ts225 | 109540.67 | 109257.00 | 0.258 |
| g_tsp225 | 2449.00 | 2449.00 | 0.000 |

Média entre geometrias: 0.480%; mediana: 0.258%.
Vitórias/empates/derrotas por geometria: 7/1/1. Gate pré-fixado: False.
Intervalo bootstrap exploratório por geometria: [0.071, 1.043]%.
O intervalo não corrige a seleção de método/coorte durante a pesquisa. Não alegar confirmação estatística ou publicabilidade só com este piloto.

## Tempo de retorno e prazo de aceitação

O prazo de 30 segundos limita a aceitação de soluções. O processo nativo verifica o prazo cooperativamente e pode retornar depois dele; seus incumbentes tardios são descartados. Assim, os tempos reais de retorno devem acompanhar os custos.
mdm: retorno médio 30.5904s; máximo 32.3813s; excesso médio 0.5904s; hard timeouts 0; chamadas LKH 0.
hybrid: retorno médio 29.6911s; máximo 34.3527s; excesso médio 0.1941s; hard timeouts 0; chamadas LKH 23.
