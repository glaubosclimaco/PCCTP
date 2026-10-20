# Diagnóstico descritivo do lote de 300 segundos

Estado: COMPLETE; 81/81 execuções; 9/9 geometrias completas.

Somente geometrias com as nove execuções entram nas comparações abaixo. Estas observações de treino não substituem o MDM completo de 70 iterações nem alteram o critério final de 300 segundos.

| Horizonte (s) | H10: ganho médio (%) | H60: ganho médio (%) | Geometrias |
|---:|---:|---:|---:|
| 1 | -0.4841 | -0.1626 | 9 |
| 10 | -0.0195 | 1.0457 | 9 |
| 30 | 0.9413 | 0.0757 | 9 |
| 60 | -0.0156 | 0.0668 | 9 |
| 120 | -0.6466 | 1.5764 | 9 |
| 300 | -1.2793 | 0.9794 | 9 |

Valores positivos indicam menor custo que o MDM no mesmo horizonte. Médias são calculadas por geometria antes da agregação.

| Variante | Execuções | Suporte alterado | Testemunho de cobertura | Testemunho de prêmio | Melhorou após inicialização | Melhor incumbente após 120 s |
|---|---:|---:|---:|---:|---:|---:|
| hybrid10 | 27 | 25 | 8 | 22 | 26 | 2 |
| hybrid60 | 27 | 22 | 5 | 19 | 25 | 2 |

Os testemunhos verificam se remover do suporte inicial todos os vértices ausentes na solução final, sem inserir os novos, violaria cobertura ou prêmio. Não representam a trajetória efetiva de destruição/reparação e não provam a causa do ganho. As colunas de cobertura e prêmio podem se sobrepor. O horário do melhor incumbente refere-se ao primeiro evento registrado com o custo final durante a fase posterior à inicialização.
