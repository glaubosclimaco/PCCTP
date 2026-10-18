# Redução exata do reparo por contração de caminhos fixos

Status: argumento para revisão; não há alegação de novidade desta redução.
Implementação: `experiments/ml_lns/fragment_repair.py`.

Considere um ciclo incumbente simples C e um conjunto de arestas livres E_U.
As arestas de C fora de E_U são fixadas em um; as demais arestas fora de E_U
são fixadas em zero. Somente visitas de U podem mudar; visitas obrigatórias
permanecem selecionadas. Custos são os do problema original.

## Equivalência inteira

Se todas as arestas de C forem fixas, qualquer ciclo simples conectado que as
contenha é o próprio C. Caso contrário, F = E(C) menos E_U forma uma floresta
de caminhos, incluindo vértices isolados. Contraia cada caminho não trivial
em um supervértice obrigatório, conservando os dois extremos identificados.
Conserve também cada vértice selecionado isolado e cada vértice ausente de U.
Vértices ausentes fora de U não podem ser incluídos e são eliminados.

Uma aresta livre só pode ser usada se cada extremo original tiver grau em F
menor que dois. Para cada extremo de caminho, imponha exatamente uma aresta
livre incidente. Para cada vértice isolado v, imponha grau livre igual a 2 y_v.
O supervértice de cada caminho é obrigatório. Um fluxo de uma mercadoria no
multigrafo contraído exige conectividade entre todos os supervértices
selecionados, tendo como raiz o componente de um vértice obrigatório.

Laços são necessários: fechar um único caminho com pelo menos três vértices
originais pode ser um ciclo válido. Arestas paralelas também são necessárias
quando dois caminhos são conectados pelos dois pares de extremos. O fluxo
ignora laços e conserva as arestas paralelas. Exigir pelo menos três vértices
**originais**, usando o tamanho dos componentes, impede ciclos degenerados.

Qualquer ciclo original admissível contém F, satisfaz os graus residuais e
induz um multigrafo conectado de grau dois nos componentes selecionados; ele
admite o fluxo exigido. Reciprocamente, uma solução inteira reduzida expande
para um grafo conectado, simples e de grau dois nos vértices originais. Portanto
é um único ciclo, preserva todas as fixações e satisfaz o mínimo de três visitas.
Cobertura agrega a união das incidências dos vértices de cada componente; prêmio
agrega sua soma. O custo é o custo constante de F mais o das arestas livres.
Logo, as duas formulações têm o mesmo conjunto de ciclos e valor ótimo inteiro.
Isso não afirma equivalência das relaxações lineares ou do desempenho do solver.

## Tamanho em função da vizinhança

Se U possui k vértices opcionais, a deles selecionados no incumbente, e E_U
contém todas as arestas incidentes em U mais no máximo p arestas adicionais,
no máximo 2a+p arestas do ciclo podem ser livres. Remover r>0 arestas de um
ciclo produz exatamente r componentes sobre seus vértices. Os k-a candidatos
ausentes acrescentam vértices isolados. Assim, o número de componentes ativos
é no máximo r+k-a <= k+a+p <= 2k+p. Para p=2k, são no máximo 4k.

Somente extremos de caminhos e candidatos isolados podem receber novas arestas.
Seu número é O(k+p); logo, existem O(k(k+p)+p) arestas livres utilizáveis para
esse esquema. O número de variáveis de arestas e fluxo deixa de depender
quadraticamente do tamanho total da instância. Ainda é necessário processar
os dados de cobertura/prêmio e montar as restrições correspondentes.

## Verificações disponíveis

- 15 vizinhanças confrontadas com enumeração independente de ciclos.
- 20 casos aleatórios confrontados com a formulação SCF original completa.
- Dois casos de liberação apenas de vértices ausentes.
- Backends e warm starts são fatores separados. Tempos preliminares sob carga
  concorrente não sustentam uma afirmação definitiva de aceleração.
