# Diagnóstico de seletores de padrão no mesmo banco elite

Registrado antes de ler os custos deste lote. Geometrias de desenvolvimento, semente 42, última mineração de cada execução. O MDM corrigido, com decodificador inicializado e elite delimitada, só é observado: a construção e a busca não mudam.

A comparação de qualidade usa apenas padrões com a mesma quantidade de arestas, igual à do maior padrão devolvido pelo minerador naquela chamada. A utilidade é a fração de demandas cobertas pela união dos vértices, mais a fração de prêmio limitada a 1, menos a fração de vértices opcionais que podem sair sem perder cobertura, menos o custo do padrão dividido pelo número de arestas e pela mediana dos custos positivos entre visitas. Cada escore desempata pelo índice mais cedo no arquivo, sem usar o escore de outro seletor.

Os controles são tamanho, soma das frequências direcionadas das arestas, suporte do padrão inteiro, o último padrão do arquivo (primeira escolha do MDM) e uma escolha aleatória com semente derivada da instância, do índice da mineração e de 20261018. O custo do último padrão só entra na comparação quando ele tem a mesma quantidade de arestas.

Um modelo inviável comprovado é pior do que um custo finito comprovado. O par só conta se os dois desfechos estão comprovados. Se a utilidade escolher o mesmo padrão que a frequência marginal, ela não muda a decisão. Se escolher outro e o custo restrito não for estritamente melhor na maioria das geometrias comparáveis, o sinal não justifica montar o método completo. Este lote não é teste reservado e não serve para ajustar a fórmula.

## Resultado

As nove geometrias completaram a busca e mineraram. Quatro não tinham dois padrões do tamanho máximo, então não entram na comparação de qualidade: berlin52, eil51, eil76 e gr48. gr48 devolveu um único padrão.

Nas cinco geometrias comparáveis, todos os modelos restritos foram resolvidos até a otimalidade. Utilidade contra frequência marginal: 1 melhor, 3 empates, 1 pior. Contra suporte, tamanho e aleatório: 0 melhores, 4 empates, 1 pior. Contra a primeira escolha do MDM: 1 melhor, 3 empates, 1 pior.

O único custo estritamente melhor é att48, 2652 contra 2676 da frequência marginal e da primeira escolha do MDM. Nesse caso a utilidade escolheu o mesmo padrão que tamanho e suporte. O único padrão diferente desses controles é pr76: a utilidade sobe de 0,5914 para 0,6262 e o custo restrito piora de 66012 para 66020. Nos demais empates o custo ótimo coincide mesmo quando o padrão muda. Todos os padrões do tamanho máximo têm suporte 2, o mínimo pedido ao minerador, então o suporte conjunto não separa esses candidatos.

A regra registrada não é satisfeita. Esta utilidade não justifica implementar a matheurística completa. A fórmula não será ajustada neste lote.
