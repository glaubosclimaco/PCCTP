# Desigualdades estruturais: descobertas e limites

## 1. Resultado negativo: corte único incondicional

Considere S sem a raiz, e o teste “nenhuma seleção em V\S satisfaz obrigatórios, cobertura e prêmio”. Como os prêmios são não negativos e a cobertura é monotônica, se existir uma seleção viável em V\S, a seleção de TODOS os vértices de V\S também satisfaz essas três condições. Portanto, a inviabilidade conjunta ocorre exatamente quando falta um obrigatório, existe w com R_w contido em S, ou p(V\S)<Q. Esses são precisamente os três certificados já disponíveis. Logo, um oráculo conjunto de seleção não cria nova família x(delta(S))>=2. Atenção: tamanho mínimo do ciclo, esparsidade, restrições superiores e custos mudam a premissa; não estão nesse resultado negativo.

Também é inválido inferir “preciso de k vértices dentro de S, então preciso cruzar delta(S) pelo menos 2k vezes”: uma excursão única pode visitar todos os k vértices e retornar com apenas duas cruzadas.

## 2. Teorema de posto de partição com alternativas de cobertura/prêmio

Particione V em P_0,...,P_k, com raiz r em P_0. Defina rho(P) como o menor número de blocos não raiz que uma seleção satisfazendo T, cobertura e quota deve ativar. Para rho>=1, é válida

    x(E_entre(P)) >= rho(P)+1.

Se rho=0, a RHS segura é 0, não 1.

Prova: contraia cada bloco visitado do ciclo, remova laços e preserve multiplicidades de arestas. Se b>=2 blocos são visitados, o multigrafo resultante é conexo e euleriano; cada um dos b vértices tem grau par pelo menos2. Logo há pelo menos b arestas entre blocos. Como a raiz está visitada e b>=rho+1, segue a tese. Mesmo se o ciclo reentrar muitas vezes num bloco, multiplicidades tornam a prova válida.

Oracle para rho: variável binária z_i por bloco não raiz. Obrigatórios ativam seus blocos; para cada demanda não coberta por P0, impor soma de z_i para blocos que intersectam Rw >=1; prêmio p(P0)+soma p(Pi)z_i >=Q; minimizar soma z_i. O oráculo considera todos os vértices de um bloco ativado disponíveis, o que é exato para seleção monotônica (desconsiderando tamanho mínimo do circuito). Um limite inferior certificado do oráculo também gera corte válido; incumbente do oráculo não fornece o rho inferior necessário. Arredondar um bound LP para cima é válido por integralidade de rho.

Status científico: teorema válido; novidade bibliográfica NÃO estabelecida. É próximo a rank/partition inequalities de roteamento por grupos. Seu diferencial candidato é o rank conjunto induzido por cobertura sobreposta e quota. Não chamar nova família antes de confrontar GTSP, group TSP, covering salesman e generalized subtour/polyhedral cuts.

### Testemunha estrita, mesmo com todo o casco de seleção y

Nove vértices: raiz0 e quatro blocos B1={1,2}, B2={3,4}, B3={5,6}, B4={7,8}. T={0}; há seis demandas, uma para cada par de blocos, com Rw=Bi união Bj. Quota1; prêmio1 em cada opcional e0 na raiz.

Cobrir todos os pares de blocos exige ativar ao menos3 dos4 blocos: se dois estivessem vazios, a demanda correspondente não seria coberta. Portanto rho=3 e RHS4.

Ponto racional: yi=5/8 para cada opcional, y0=1; x de raiz a cada opcional=1/4; x interno de cada bloco=5/8; x entre vértices de blocos opcionais diferentes=1/16. A soma entre blocos é7/2, violando4.

O vetor y pertence ao casco convexo de TODAS as seleções válidas: escolher uniformemente3 de4 blocos; em cada bloco ativado escolher ambos os vértices com probabilidade2/3, só o primeiro com1/6 e só o segundo com1/6. Cada vértice tem marginal (3/4)*(5/6)=5/8. Toda seleção cobre as seis demandas e satisfaz quota. Assim, nenhum corte puramente em y, mesmo desconhecido, elimina esse ponto. A desigualdade envolve realmente a estrutura de roteamento agregado.

O verificador exato em experiments/pcctp_cuts/verify_rank_witness.py usa somente Fraction e stdlib. Resultado em rank_witness_result.json. Verifica graus, linking, cobertura, quota,2048 SEC,134 cortes forçados por cobertura,2 cortes por prêmio, todos os154732800 subconjuntos ímpares de arestas de dentes via programação dinâmica de paridade,1 knapsack cover,189 seleções válidas e54624 circuitos distintos até reversão. Todos passam as famílias antigas; todos os circuitos válidos satisfazem o corte novo. O teste dos blossoms usa a forma forte correta x(delta(H))>=2x(F)-|F|+1, incluindo F não matching; isso domina o teste apenas de dentes disjuntos e evita depender de detalhes tipográficos da fórmula IJISE.

Nota: para partição em singletons, o corte se reduz a cardinalidade de seleção e identidade global de graus. Não é o caso desta testemunha: massa opcional y=5, mas a deficiência está nos blocos ativados.

### Realização métrica e objetivo

Fixe custo1 entre blocos de V e custo epsilon=1/10 dentro de cada Bi. Distâncias de W a V são1 quando há cobertura,2 caso contrário (raiz não cobre W); entre demandas W, distância2. Com D=1, são as vizinhanças requeridas. Essa matriz é métrica: dois vértices de um mesmo Bi têm idêntica pertinência em todos os Rw, e as únicas arestas menores que1 são internas desses blocos. Qualquer percurso viável visita raiz+3 blocos e custa pelo menos4; visitar um vértice de cada um de3 blocos dá custo4. O ponto fracionário tem custo7/2+(5/2)*(1/10)=15/4. A testemunha dá separação objetiva métrica, sem afirmar realização Euclidiana plana.

### Separação e riscos

Subfamílias por grafo de obrigações: para a partição fixa, criar vértice para cada bloco não raiz. Acrescentar aresta ij quando algum Rw está contido em Bi união Bj (e não intersecta P0), ou quando p(V\(Bi união Bj))<Q. A primeira condição exige que pelo menos um dos dois blocos seja visitado para cobrir w; a segunda exige visita a pelo menos um para atingir quota. Blocos com obrigatórios, ou com um Rw inteiro, ou com complemento de prêmio insuficiente, são obrigações unárias. Logo, o conjunto de blocos ativados é vertex cover desse grafo que respeita as obrigações unárias. Um bound certificado l para essa cobertura de vértices fornece rho>=l e, quando l>=1, x(Eentre)>=l+1.

Se esse grafo contém clique Kq, q>=2, ao menos q-1 blocos são ativados e x(Eentre)>=q. Se contém ciclo ímpar de comprimento2k+1, ao menos k+1 blocos são ativados e x(Eentre)>=k+2. A testemunha usa exatamente a clique K4. Tais consequências têm provas imediatas (independent set máximo1 na clique; no ciclo, somar z_i+z_j>=1 e arredondar para cima). Encontrar clique de tamanho arbitrário é difícil em geral; para q fixo pequeno, enumerar q-tuplas é polinomial. Em grafos de obrigações bipartidos, tau é computável por matching; em poucos blocos, resolver vertex cover exato é viável. Não confundir dificuldade de obter a partição com dificuldade do bound após fixá-la.

Pode-se fortalecer o oráculo com capacidade de prêmio por bloco p(Bi): soma p(Bi)zi+p(P0)>=Q; combinar esse knapsack com cobertura/obrigações pode aumentar rho. Emissão de clique/ciclo constitui especialização candidata de rank-partition, com novidade bibliográfica pendente.

Para uma partição fixada: obter rho certificado e testar custo interblocos é direto. Separação global sobre partições não foi caracterizada aqui; não alegar polinomialidade. Heurística inicial: contrair arestas de x alto ou construir comunidades, limitar a4–10 blocos, resolver o pequeno set cover/quota em blocos e emitir somente violações certificadas. Ranking por ganho/tempo é experimental, não parte do teorema.

## 3. Comb ativado por conjuntos obrigatoriamente visitados

### Adendo: testemunha genuinamente mista de cobertura e quota

O arquivo experiments/pcctp_cuts/verify_mixed_partition_witness.py mantém x,y e os quatro blocos da testemunha anterior, mas usa SOMENTE cinco demandas: todos os pares exceto B3 união B4. Prêmio por vértice de B1 ou B2 é1/2; de B3 ou B4 é1; quota3. Cobertura sozinha admite visitar apenas blocos1,2; quota sozinha admite blocos3,4. Conjuntamente, duas escolhas de blocos nunca bastam: a dupla1,2 tem prêmio máximo2<3; todas as demais duplas omitem alguma demanda dos cinco pares presentes. Portanto ranks separados são2 e2, mas rank conjunto é3. O verificador confirma esses três valores por enumeração exata de todas as seleções.

O vetor y pertence mesmo ao casco inteiro da seleção CONJUNTA: escolher3 de4 blocos uniformemente e, dentre seus6 vértices, omitir1 uniformemente. Os24 termos têm peso1/24 e marginal (3/4)*(5/6)=5/8. Cada seleção cobre todas as demandas, e seu prêmio é ao menos3: qualquer tripla de blocos tem prêmio total pelo menos4 e a exclusão retira no máximo1.

Verificação exata atual:2048 SEC,116 cortes por cobertura,212 cortes por prêmio,106 cover cuts,106 extended cover cuts da fórmula IJISE,106 residual-knapsack-cover cuts (truncados pela quota residual),1878 CI,154732800 odd-F blossoms certificados por DP,141 seleções e54336 circuitos viáveis. Todos os testes passam. A matriz métrica de14 locais (9 visitáveis+5 demandas) satisfaz2744 triângulos. Custo fracionário15/4; ótimo inteiro enumerado41/10. O corte de partição é7/2>=4 e elimina o ponto; não se afirma que esse único corte feche todo o gap objetivo, pois o ótimo é4.1 e o argumento de cruzadas fornece4.

Os resultados reproduzíveis estão em mixed_rank_witness_result.json. Essa variante demonstra o benefício de combinar obrigações de cobertura e prêmio em BLOCOS, apesar do resultado negativo para simples cortes de um só conjunto S.

Natureza poliédrica precisa no exemplo K4: o corte é uma desigualdade de Chvátal–Gomory de primeira rodada. Seja A a soma de arestas raiz–blocos e B a soma entre blocos opcionais. Somar as seis desigualdades delta(Bi união Bj)>=2 produz3A+4B>=12. A equação de grau da raiz fornece A=2. Multiplicar aquela soma por1/4 e adicionar A/4=1/2 resulta A+B>=7/2; como A+B é inteiro nas soluções, arredondar fornece A+B>=4. Portanto o corte NÃO é combinação não negativa sem arredondamento dos cortes velhos, mas deriva deles com uma rodada de integralização. Isso orienta uma contribuição honesta como separação estruturada de CG/rank de partição, em vez de sugerir que toda a força vem de informação lógica ausente do modelo. Novidade conceitual genérica de CG não existe; possível novidade está na família especializada, caracterização, separação e evidência computacional.

Escolha handle H e número ímpar t>=3 de dentes disjuntos Ti. Para cada dente, Ai=Ti interseção H e Bi=Ti\H devem ambos conter algum vértice visitado em TODA solução válida. Um certificado suficiente é conter um obrigatório, um Rw completo ou um conjunto cujo complemento não atinge a quota. Então

    x(delta(H))+soma_i x(delta(Ti)) >=3t+1.

Prova autônoma: como cada Ai e Bi é visitado, x(delta(Ai)),x(delta(Bi))>=2; como existem outros dentes visitados, x(delta(Ti))>=2. A identidade delta(Ai)+delta(Bi)=delta(Ti)+2x(Ai,Bi) implica, somando metade das três desigualdades, delta(Ti)+x(Ai,Bi)>=3. Os conjuntos de arestas Ai–Bi são disjuntos e contidos em delta(H). Somando, LHS>=3t. Todo corte de um ciclo tem cardinalidade par; LHS é par, enquanto3t é ímpar, logo LHS>=3t+1.

É uma aplicação/lift de comb conhecido, com prova de arredondamento por paridade; não é soma não negativa sem arredondamento das desigualdades antigas. Contudo, o exemplo canônico de6 vértices com dentes de2 nós já é 2matching e NÃO demonstra avanço sobre IJISE. Não obtivemos aqui testemunha estrita contra todas as 2matching para dentes maiores. Prioridade teórica menor que o corte de partição, cuja testemunha está verificada.

## 4. Fontes primárias e auditoria de novidade

- [Grötschel e Padberg, 1978, comb/generalized comb](https://www.zib.de/userpage/groetschel/pubnew/paper/groetschelpadberg1978.pdf): antecedente clássico; impede apresentar comb padrão como descoberta.
- [Fischetti, Salazar González e Toth, 1997, GTSP branch-and-cut](https://doi.org/10.1287/opre.45.3.378): referência obrigatória para cortes por clusters e famílias generalizadas.
- [Maziero, Usberti e Cavellucci, covering salesman, preprint2021](https://arxiv.org/abs/2104.01173), [texto integral](https://arxiv.org/html/2104.01173): traduz cortes GTSP para cobertura sobreposta e propõe Cover Intersection inequalities (CI), envolvendo união de dois cutsets. O verificador foi ampliado para testar exatamente a transposição de todas as CI usando os seis conjuntos Rw e o singleton obrigatório raiz. A testemunha passa também essas CI; isso mostra independência para essa transposição, NÃO dominância entre famílias em geral.
- [Prize Collecting TSP, primitive comb, 2016](https://www.scitepress.org/Papers/2016/57581/57581.pdf): evidencia antecedentes de generalized SEC e primitive comb em coleta de prêmios; verificar detalhes antes de nomear a família.

Conclusão epistemológica: há teorema válido, ponto racional de independência contra as famílias PCCTP listadas e protocolo reproduzível. Ainda faltam equivalência com famílias GTSP/CSP conhecidas, análise de facets, complexidade de separação, teste de separação nas instâncias reais e impacto computacional. Nenhuma melhora computacional é afirmada.
