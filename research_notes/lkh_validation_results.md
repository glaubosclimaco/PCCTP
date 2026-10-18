# MDM10s versus MDM9s+LKH

7 geometrias, três sementes pareadas; validação reservada com método congelado. Ambos recebem 10 segundos totais. Overruns são registrados; melhorias tardias do LKH são descartadas.

| Geometria | Média MDM | Média híbrido | Ganho pareado médio % |
|---|---:|---:|---:|
| g_d198 | 13724.67 | 13680.67 | 0.321 |
| g_kroB100 | 18406.00 | 18406.00 | 0.000 |
| g_kroE100 | 10586.00 | 10586.00 | 0.000 |
| g_pr124 | 34359.00 | 34359.00 | 0.000 |
| g_pr152 | 49198.67 | 49198.67 | 0.000 |
| g_pr264 | 47589.33 | 45442.00 | 4.433 |
| g_pr299 | 31911.00 | 31433.00 | 1.499 |

Média entre geometrias: 0.893%; mediana: 0.000%.
Vitórias/empates/derrotas por geometria: 3/4/0. Gate pré-fixado: False.
Intervalo bootstrap exploratório por geometria: [0.046, 2.160]%.
O intervalo não corrige a seleção de método/coorte durante a pesquisa. Não alegar confirmação estatística ou publicabilidade só com este piloto.

## Tempo de retorno e prazo de aceitação

O prazo de 10 segundos limita a aceitação de soluções. O processo nativo verifica o prazo cooperativamente e pode retornar depois dele; seus incumbentes tardios são descartados. Assim, os tempos reais de retorno devem acompanhar os custos.
mdm: retorno médio 10.1520s; máximo 10.6243s; excesso médio 0.1520s; hard timeouts 0; chamadas LKH 0.
hybrid: retorno médio 9.2742s; máximo 9.8056s; excesso médio 0.0000s; hard timeouts 0; chamadas LKH 21.
