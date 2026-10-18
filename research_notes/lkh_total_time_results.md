# Piloto MDM10s versus MDM9s+LKH

9 geometrias, três sementes pareadas; treino exploratório, não confirmatório. Ambos recebem dez segundos totais. Overruns são registrados; melhorias tardias do LKH são descartadas.

| Geometria | Média MDM | Média híbrido | Ganho pareado médio % |
|---|---:|---:|---:|
| g_a280 | 2131.33 | 2045.33 | 3.963 |
| g_gil262 | 2091.00 | 2042.67 | 2.306 |
| g_gr229 | 101028.67 | 99929.33 | 1.082 |
| g_kroA200 | 23613.67 | 23593.67 | 0.085 |
| g_pr226 | 69306.00 | 69009.67 | 0.422 |
| g_rd400 | 11683.67 | 11402.33 | 2.471 |
| g_si175 | 7258.00 | 7129.00 | 1.780 |
| g_ts225 | 109592.67 | 109309.00 | 0.258 |
| g_tsp225 | 2454.00 | 2454.00 | 0.000 |

Média entre geometrias: 1.374%; mediana: 1.082%.
Vitórias/empates/derrotas por geometria: 8/1/0. Gate pré-fixado: True.
Intervalo bootstrap exploratório por geometria: [0.601, 2.244]%.
O intervalo não corrige a seleção de método/coorte durante a pesquisa. Não alegar confirmação estatística ou publicabilidade só com este piloto.

## Tempo de retorno e prazo de aceitação

Os dez segundos limitam a aceitação de soluções. O processo nativo verifica o prazo cooperativamente e pode retornar depois dele; seus incumbentes tardios são descartados. Assim, os tempos reais de retorno devem acompanhar os custos.
mdm: retorno médio 10.3013s; máximo 11.3296s; excesso médio 0.3013s; hard timeouts 0; chamadas LKH 0.
hybrid: retorno médio 9.4706s; máximo 10.1059s; excesso médio 0.0044s; hard timeouts 0; chamadas LKH 26.
