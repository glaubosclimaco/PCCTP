# Projeção exata de seleção para o operador route_context

Status: lema específico do operador implementado, com prova para revisão e
teste finito; **novidade na literatura não estabelecida**. Não se aplica
automaticamente a outros operadores ou a grafos incompletos.

Seja C um ciclo simples incumbente em grafo completo não direcionado, com
vértices obrigatórios não vazios. Libere visitas de U, apenas opcionais,
|U|=k>=1, todas as arestas incidentes em U e as arestas adicionais escolhidas
por `route_context`, com orçamento 2k. Fora de U, mantenha as visitas de C.

**Afirmação.** A projeção das rotas admissíveis desse reparo sobre as visitas
de U coincide com as atribuições que satisfazem cobertura, prêmio e pelo menos
três visitas originais, mantendo todas as visitas externas. A hipótese pressupõe
as restrições do PCCTP usado aqui, sem outras restrições de roteamento.

Necessidade é imediata. Para suficiência, chame B=V(C) menos U, não vazio pois
U não contém obrigatórios. A seleção desejada é B união A, A contido em U.

Se |B|>=3, remova provisoriamente todos os vértices de U do ciclo. Cada corrida
maximal removida requer uma aresta de desvio entre extremos retidos distintos.
Há no máximo |U interseção V(C)|<=k dessas arestas. Todas recebem prioridade
zero no operador; portanto entram no orçamento 2k, mesmo se o número de arestas
elegíveis for menor, pois nesse caso entram todas. As demais arestas do ciclo
reduzido já pertenciam ao incumbente. Se havia uma corrida removida, existe
ao menos um desvio livre que pode receber inserções. Insira os vértices de A
nesse desvio em qualquer ordem: cada aresta nova é incidente em U e está livre.
Inserções subsequentes podem usar as novas arestas, também livres.

Se nenhum vértice incumbente foi removido, U contém apenas candidatos ausentes.
Se A é vazio, conserve C. Caso contrário o operador atribui prioridade um às
arestas incumbentes, liberando pelo menos uma. Insira todos os vértices de A
nesse intervalo; todas as novas arestas são incidentes em U.

Se |B|=2, a condição de três visitas exige |A|>=1. A aresta entre os dois
vértices de B é incumbente ou é um desvio prioritário livre. Ela e um vértice
de A formam um triângulo; insira os demais pelas arestas incidentes em U.
Se |B|=1, |A|>=2, e um triângulo com dois vértices de A usa somente arestas
incidentes em U; os demais entram da mesma forma. Em todos os casos, arestas
fixadas em um são preservadas e nenhuma aresta fixada em zero é introduzida.

Cobertura e prêmio dependem da seleção e são preservados pela construção.
A prova não fornece limite válido para o custo de reconexão. O ganho aditivo
usado na heurística é um substituto, não um bound do custo real.

`selection_lift.py` implementa a construção e verifica todas as fixações.
`verify_selection_lift.py` testou977 atribuições viáveis em80contextos aleatórios,
incluindo |B|=1,2,3. O teste não substitui a revisão da prova.

Consequências a investigar: fornecer soluções iniciais viáveis ao reparo,
separar exatidão da viabilidade e heurística do custo, gerar trocas conjuntas
antes do MIP, e treinar seleção de vizinhanças/atribuições com controles fortes.
Tudo ainda deve ser confrontado com a literatura e com o algoritmo completo.

## Refinamento: orçamento k suficiente e necessário no pior caso

A mesma prova funciona com orçamento **k**, em vez de2k. Se há vértices
incumbentes liberados e pelo menos três retidos, existem no máximo k desvios
prioritários. Sem incumbentes liberados, basta liberar uma aresta incumbente
para inserir candidatos ausentes. Os casos com um ou dois retidos continuam
pela construção de triângulos descrita acima. O operador existente aceita
`extra_edges=k`, sem alteração da ordem de prioridades.

O limite k é necessário no pior caso para garantir toda a projeção de seleção.
Construa um ciclo em grafo completo com k opcionais U e k+2 obrigatórios B.
Intercale cada opcional entre os primeiros k pares consecutivos de B. Os dois
últimos arcos de fechamento são arestas entre obrigatórios. Remova todos os
opcionais; essa seleção é viável se B cobre tudo e satisfaz a quota (por exemplo,
quota zero). O ciclo final contém k+2 arestas entre obrigatórios; apenas duas
eram arestas incumbentes fora de U. Todas as outras estão fixadas em zero se não
forem liberadas. Logo, qualquer aumento que permita essa seleção precisa liberar
ao menos k arestas não incidentes em U. O argumento inclui k=1, pois |B|=3.

Assim k é o orçamento universal mínimo no pior caso, sob as hipóteses acima.
Isso não diz que k é mínimo em toda instância nem que preserva o mesmo ótimo
de roteamento obtido com2k arestas extras. A contração passa a admitir o limite
de3k componentes ativos em vez de4k.

O teste foi ampliado para os dois orçamentos:1954 construções verificadas
sobre977 atribuições viáveis nos mesmos80contextos.
