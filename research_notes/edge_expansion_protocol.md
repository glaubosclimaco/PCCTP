# Expansão de arestas — protocolo exploratório, 19/10/2026

Antes da coleta: nove geometrias hardTRAIN, mesmos conjuntos de rotas e mesma
melhor partida do edge_pool_study original, sem incorporar sua saída melhorada.
Duas expansões: união com cinco vizinhos mais baratos por vértice (empates por
índice); controle adiciona exatamente a mesma quantidade de arestas ausentes,
sorteadas uniformemente sem reposição com seed42. Nenhuma variável de visita
é fixada. Arestas originais sempre preservadas. HiGHS seed42, uma thread,60s
por braço incluindo leitura e construção; registrar CPU e CLOCK_MONOTONIC.
Alternar ordem dos braços entre geometrias. Mesmos dados offline: custo de
geração das rotas não incluído. Comparação ao pool puro já coletado é diagnóstica.

Gate: expansão geométrica melhora pelo menos três geometrias e média>=0.2%
sobre melhor entrada, e supera expansão aleatória em pelo menos cinco geometrias,
sem perda de factibilidade. Caso contrário preservar negativo; não abrir teste.
Ótimos e bounds são restritos ao respectivo grafo. Resultado após prazo descartado
com registro bruto. Não alegar método novo baseado apenas em vizinhos próximos.
