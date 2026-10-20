# Diagnóstico descritivo do lote de 300 segundos

Estado: PARTIAL; 30/81 execuções; 3/9 geometrias completas.

Somente geometrias com as nove execuções entram nas comparações abaixo. Estas observações de treino não substituem o MDM completo de 70 iterações nem alteram o critério final de 300 segundos.

| Horizonte (s) | H10: ganho médio (%) | H60: ganho médio (%) | Geometrias |
|---:|---:|---:|---:|
| 1 | 0.0000 | 0.0000 | 3 |
| 10 | -1.6989 | 0.0767 | 3 |
| 30 | 3.9764 | -0.6542 | 3 |
| 60 | 3.2961 | 0.1947 | 3 |
| 120 | 1.9500 | 1.9738 | 3 |
| 300 | 1.4642 | 1.5076 | 3 |

Valores positivos indicam menor custo que o MDM no mesmo horizonte. Médias são calculadas por geometria antes da agregação.

| Variante | Execuções | Suporte alterado | Testemunho de cobertura | Testemunho de prêmio | Melhorou após inicialização | Melhor incumbente após 120 s |
|---|---:|---:|---:|---:|---:|---:|
| hybrid10 | 9 | 9 | 4 | 9 | 9 | 0 |
| hybrid60 | 9 | 9 | 2 | 9 | 9 | 2 |

Os testemunhos verificam se remover do suporte inicial todos os vértices ausentes na solução final, sem inserir os novos, violaria cobertura ou prêmio. Não representam a trajetória efetiva de destruição/reparação e não provam a causa do ganho. As colunas de cobertura e prêmio podem se sobrepor. O horário do melhor incumbente refere-se ao primeiro evento registrado com o custo final durante a fase posterior à inicialização.
