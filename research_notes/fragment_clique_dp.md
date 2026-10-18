# Reparo exato por componentes e programação dinâmica

Status: implementação e prova para revisão; novidade não estabelecida. Held–Karp
e contração de caminhos são conhecidos. Não confundir este operador com
`route_context`, que libera outras arestas incumbentes e pode gerar muitos
componentes, tornando a DP inviável.

## Operador

Seja J um conjunto de k vértices de liberação de roteamento, que pode incluir
opcionais U e obrigatórios R. Obrigatórios continuam visitados; somente opcionais
de U podem mudar de seleção. Remova temporariamente todas as arestas incumbentes
incidentes em J, conservando as demais arestas do ciclo como F. Se J contém
apenas vértices ausentes, abra uma aresta incumbente para permitir inserções.

Considere os vértices de J e todos os extremos dos caminhos de F. Libere todas
as arestas entre esses vértices de fronteira, **exceto as que pertencem a F**.
Essa exceção importa para caminhos de duas visitas, cujos dois extremos são
adjacentes: sua aresta interna deve continuar fixa. Os demais vértices de F
têm grau dois fixo. A construção é linear no tamanho da rota mais O(k²), além
do acesso aos custos; não precisa enumerar todas as arestas do grafo original.

## Limite de tamanho

Se a=|J interseção V(C)|>0 e existem r corridas de visitas liberadas no ciclo,
remover suas arestas deixa a vértices liberados isolados e r caminhos retidos
(se houver retidos). Os k-a candidatos ausentes também são isolados. Portanto
há q=k+r<=2k componentes ativos. Se todos os visitados estiverem em J, q=k.
Se a=0, a abertura adicional deixa um único caminho retido e q=k+1<=2k.
A fronteira possui no máximo k+2r<=3k vértices. Isso limita o tamanho combinatório
do reparo, mas processamento de cobertura/prêmio ainda depende dos dados de entrada.

Todo componente não trivial é obrigatório porque contém arestas fixas. Cada
obrigatório liberado é um componente obrigatório isolado. Visitas opcionais
liberadas são componentes isolados selecionáveis. A projeção de seleção é exata:
uma seleção viável de componentes pode ser conectada em qualquer ordem pela
fronteira completa, respeitando as arestas internas e o mínimo de três visitas.
Laços de fechamento de um único caminho e conexão de dois caminhos precisam
ser tratados como no argumento de equivalência anterior.

## DP

Cada caminho admite duas orientações; um isolado admite uma. Fixe uma orientação
do componente raiz, sem perda em grafo não direcionado com custos simétricos.
D[S,j,o] guarda o menor custo de conectar componentes em S, começando na raiz
e terminando no componente j com orientação o. A transição acrescenta a aresta
livre entre a saída de j e a entrada do próximo componente. Fechamentos são
avaliados apenas para máscaras que incluem todos os componentes obrigatórios,
cuja seleção satisfaz cobertura, prêmio e pelo menos três visitas originais.
Somam-se os custos constantes de F ao final.

A recorrência enumera todos os ciclos da vizinhança por sua sequência de
componentes e orientações. Complexidade O(q²2^q), memória O(q2^q), mais o custo
de enumerar a viabilidade das seleções opcionais. Não chamar isso de kernel
polinomial: o número de máscaras pode ser exponencial no parâmetro.

Custos não negativos permitem podar prefixos cujo custo já alcança o melhor
fechamento conhecido. A implementação começa com o incumbente como limite e
fallback. Só retorna certificado de ótimo ao terminar a enumeração. Limites
de tempo ou de componentes preservam uma rota viável, sem certificado.

## Verificações e diagnóstico

- `verify_dp.py`:30comparações com SCF completo,18vizinhanças por enumeração
  independente, teste de timeout com a raiz na posição zero. Todos passaram.
- Variante inicial de DP sobre `route_context`:108/108casos excederam o limite
  de18componentes. Resultado negativo preservado em `dp_development/`.
- `fragment_clique`: em108consultas de três geometrias, DP certificou106 e
  MIP104; duas consultas excederam o limite de componentes da DP.
- A vetorização da seleção reduziu a mediana observada de38ms para17ms; a
  mediana do MIP pareado foi36ms. São medições de desenvolvimento sob carga
  concorrente, não uma afirmação definitiva de aceleração do algoritmo completo.
- Outputs anteriores à correção das arestas internas estão segregados em
  `fragment_clique_pre_boundary_fix/` e não sustentam a avaliação atual.

Próximo passo: avaliar trajetórias completas com seleções e âncoras de rota,
medir o custo de prever/escolher vizinhanças e comparar modelos de aprendizado
com controles que também conheçam o número de componentes e o limite de DP.
