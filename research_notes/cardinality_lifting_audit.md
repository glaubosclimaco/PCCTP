# Cardinalidade e lifting: confronto de hipóteses

18 de outubro de 2026. A busca não estabelece novidade por ausência de resultados.

## Fonte primária e trechos localizados

Rüdiger Stephan, *Polyhedral Aspects of Cardinality Constrained Combinatorial
Optimization Problems*, tese de doutorado, TU Berlin, 2009.
PDF institucional: https://webdoc.sub.gwdg.de/ebook/dissts/Berlin_TU/Stephan2009.pdf

Seção 3.6, corolário 3.53(f), página impressa 113: o limite inferior de
cardinalidade é faceta do politopo não dirigido de ciclos com cardinalidades
permitidas. Seção 3.2, teorema 3.11, página 60: clonagem de nó interno em
caminhos dirigidos de cardinalidade fixa, sob regularidade e condição sobre
bowties, preserva faceta com coeficientes calculados. A seção 3.5 também
contém partições com RHS 1, distintas da fórmula de posto aqui estudada.

## O que já podemos atribuir

Com blocos singleton e q opcionais, nossa seleção exige pelo menos q-1
opcionais e inclui a raiz: os comprimentos permitidos são q e q+1.
O corte é x(E)>=q. Seu formato é precisamente o limite clássico de
cardinalidade; não o apresentar como nova desigualdade.

A raiz obrigatória define adicionalmente a face x(delta(r))=2 do politopo
de ciclos sem raiz fixa. A propriedade de faceta relativa a essa face NÃO
decorre automaticamente da propriedade de faceta no politopo sem raiz.
A prova direta escrita no projeto cobre essa restrição; sua novidade,
porém, continua não estabelecida.

## O que não foi deduzido da fonte

O teorema 3.11 mantém cardinalidade fixa do caminho e usa orientação e
condições explícitas de regularidade/bowtie. Na construção PCCTP, inserir
vértices consecutivos dentro de um bloco aumenta a cardinalidade da rota,
mas mantém o número de blocos ocupados e os cruzamentos. Não é legítimo
invocar aquele teorema apenas porque ambos os argumentos duplicam nós.
Também não é legítimo afirmar que essa diferença prova novidade.

O lifting de visitas fixadas na proposição 3.1 de Sundar–Rathinam
arXiv:1508.01813, já examinado, libera y num conjunto de vértices fixo.
A equivalência com a duplicação de blocos não foi demonstrada.

## Decisão editorial

Apresentar, se usado, o caso singleton como antecedente conhecido. O resultado
candidato é a caracterização de faceta na extensão por blocos e suas hipóteses,
não a cardinalidade em si, nem a realização mista do mesmo politopo.
Não chamar a família inédita antes de conferir o lifting de GTSP1995 ou
outra fonte primária que cubra duplicação com múltiplas visitas por bloco.
Nenhum desses fatos altera os resultados negativos de ML ou de aceleração.
