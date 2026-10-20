# Diagnóstico descritivo do lote de 300 segundos

Estado: PARTIAL; 72/81 execuções; 8/9 geometrias completas.

Somente geometrias com as nove execuções entram nas comparações abaixo. Estas observações de treino não substituem o MDM completo de 70 iterações nem alteram o critério final de 300 segundos.

| Horizonte (s) | H10: ganho médio (%) | H60: ganho médio (%) | Geometrias |
|---:|---:|---:|---:|
| 1 | -0.1019 | 0.2884 | 8 |
| 10 | -0.0220 | 1.1764 | 8 |
| 30 | 1.0913 | 0.0835 | 8 |
| 60 | -0.0091 | 0.0700 | 8 |
| 120 | -0.7257 | 1.7751 | 8 |
| 300 | -1.4290 | 1.1121 | 8 |

Valores positivos indicam menor custo que o MDM no mesmo horizonte. Médias são calculadas por geometria antes da agregação.

| Variante | Execuções | Suporte alterado | Testemunho de cobertura | Testemunho de prêmio | Melhorou após inicialização | Melhor incumbente após 120 s |
|---|---:|---:|---:|---:|---:|---:|
| hybrid10 | 24 | 23 | 6 | 20 | 24 | 2 |
| hybrid60 | 24 | 21 | 5 | 18 | 24 | 2 |

Os testemunhos verificam se remover do suporte inicial todos os vértices ausentes na solução final, sem inserir os novos, violaria cobertura ou prêmio. Não representam a trajetória efetiva de destruição/reparação e não provam a causa do ganho. As colunas de cobertura e prêmio podem se sobrepor. O horário do melhor incumbente refere-se ao primeiro evento registrado com o custo final durante a fase posterior à inicialização.
