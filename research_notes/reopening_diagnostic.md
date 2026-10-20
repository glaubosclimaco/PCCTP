# Diagnóstico de reabertura no padrão usado pelo MDM

Registrado antes de ler os custos desta comparação. A utilidade que escolhe o padrão inteiro já foi rejeitada. A pergunta agora é outra: soltar duas arestas desse padrão, por uma regra fixa, melhora o ótimo restrito mais do que soltar duas arestas ao acaso.

O padrão é o último do arquivo do minerador na última chamada de cada geometria, semente 42. É a primeira semente de construção do MDM depois da mineração. São liberadas exatamente 2 arestas. O restante do padrão continua obrigatório no modelo.

A prioridade de liberação é o custo da aresta dividido pela mediana dos custos positivos entre visitas, mais a média dos indicadores de redundância dos dois extremos. Um vértice opcional é redundante quando sai da união do padrão sem reduzir o número de demandas cobertas. Obrigatório não é redundante. Empate: maior custo, depois o par de índices menor. A escolha aleatória usa o SHA-256 da instância, do índice da mineração, do número 2 e de 20261018.

O critério principal só inclui geometria em que três condições valem ao mesmo tempo: o modelo sem arestas forçadas e os dois modelos com reabertura estão provados, os dois conjuntos liberados são diferentes, e o padrão inteiro fixo tem custo estritamente pior que o modelo sem força. Sem essa folga, reabrir não pode melhorar a qualidade. A regra continua somente se a liberação informada for estritamente melhor que a aleatória em mais geometrias do que o contrário. Empate não conta como sucesso. Se menos de três geometrias de desenvolvimento tiverem folga, este lote é inconclusivo e a mesma regra, sem mudança de k nem da fórmula, será repetida nas geometrias de treino. Validação e teste permanecem fechados. Não há ajuste depois de ver os números.

## Resultado

Nove geometrias, todas com mineração. Três entram na conta: att48, berlin52 e hk48. A liberação informada foi estritamente melhor em hk48, 7593 contra 7605, e empatou nas outras duas. Nas duas empatadas, soltar duas arestas não moveu o custo do padrão fixo: att48 ficou em 2676 e berlin52 em 5263. O modelo sem força provou 2650 e 5186, iguais ao custo final do MDM.

Onde o modelo sem força provou o ótimo, o MDM já o atingiu: att48 2650, berlin52 5186, brazil58 20466, eil76 332, gr48 3571, hk48 7593 e st70 419. A reabertura informada nunca ficou abaixo do MDM. Em eil51 e pr76 o modelo sem força não provou em 20 s. O padrão fixo custou 275 e 66012, e o MDM já estava em 273 e 66007.

A bandeira numérica de continuidade acende, porque houve uma vitória e nenhuma derrota. A decisão é não construir a matheurística: o único ganho iguala o MDM, e aumentar k depois de ver que duas arestas não fecharam att48 nem berlin52 seria ajuste. k permanece 2.
