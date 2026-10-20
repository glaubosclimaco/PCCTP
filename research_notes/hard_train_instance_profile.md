# Estrutura dos nove casos difíceis de treino

Descritores calculados exclusivamente das entradas; nenhum resultado de busca entra nos cálculos. Não representam os107 casos do artigo nem incluem validação ou teste.

| Instância | Obrigatórios | Opcionais | Demandas | Déficit de prêmio após obrigatórios | Demandas descobertas | LB opcionais por prêmio | LB opcionais por cobertura | LB visitados |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| g_a280_VT56_T168_W56_50 | 168 | 56 | 56 | 1004.0 | 45 | 10 | 9 | 178 |
| g_gil262_VT51_T158_W53_75 | 158 | 51 | 53 | 1923.0 | 38 | 20 | 10 | 178 |
| g_gr229_VT45_T138_W46_50 | 138 | 45 | 46 | 95968.0 | 0 | 15 | 0 | 153 |
| g_kroA200_VT40_T120_W40_25 | 120 | 40 | 40 | 3237.0 | 0 | 7 | 0 | 127 |
| g_pr226_VT44_T136_W46_50 | 136 | 44 | 46 | 185185.0 | 0 | 20 | 0 | 156 |
| g_rd400_VT80_T240_W80_25 | 240 | 80 | 80 | 6026.0 | 0 | 15 | 0 | 255 |
| g_si175_VT35_T105_W35_50 | 105 | 35 | 35 | 1356.0 | 19 | 17 | 10 | 122 |
| g_ts225_VT45_T135_W45_50 | 135 | 45 | 45 | 125640.0 | 34 | 17 | 12 | 152 |
| g_tsp225_VT45_T45_W135_75 | 45 | 45 | 135 | 1328.0 | 127 | 28 | 12 | 73 |

O limite por prêmio soma os maiores prêmios opcionais até cobrir o déficit. Nenhum conjunto com menos opcionais alcança esse total. O limite por cobertura é o teto do número de demandas ainda descobertas dividido pela maior quantidade dessas demandas que um opcional pode cobrir. A união de k opcionais cobre no máximo k vezes essa quantidade. Combina-se o máximo desses dois limites e da necessidade de pelo menos três vértices no ciclo; não se somam os limites, pois um vértice pode cumprir ambas as funções. São limites elementares de cardinalidade, sem pretensão de novidade ou de ótimo para o custo da rota.

Verificação adicional: os nove casos foram relidos pelo leitor de produção
pilot.read_instance; contagem de obrigatórios, déficit de prêmio e demandas
não cobertas pelos obrigatórios coincidiram com o leitor escalar independente.
Em gr229, kroA200, pr226 e rd400, os obrigatórios já cobrem todas as demandas.
Portanto as restrições de cobertura são redundantes nessas quatro instâncias
inteiras; não atribuir nelas eventual ganho à restauração de cobertura.
Isso não demonstra a causa dos ganhos nem altera operadores ou o protocolo.
