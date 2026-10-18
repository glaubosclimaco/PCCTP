# Partições revelam obrigações coletivas no PCCTP

**Encontramos uma desigualdade de partição válida, uma testemunha fracionária que ela separa mesmo quando toda a seleção de vértices já pertence ao seu casco convexo, e uma realização métrica com diferença estrita de custo.** Isso estabelece um fortalecimento matemático relativo a uma relaxação explicitamente testada. Também certificamos que o corte define faceta nas duas instâncias pequenas examinadas. Ainda não estabelecemos novidade bibliográfica, condições gerais de faceta, separação global eficiente ou melhoria nas instâncias dos artigos. A direção mais forte é estudar obrigações coletivas de visita a blocos, em vez de tentar aumentar arbitrariamente o número de cruzamentos de um único corte. O código do solver não foi modificado. Os verificadores usam aritmética racional e enumeram instâncias pequenas; sua função é certificar exemplos e eliminar conjecturas erradas, não estimar desempenho.

## O modelo precisa distinguir validade, força e novidade

Considere um grafo não orientado completo sobre os vértices visitáveis \(V=T\cup R\), com uma raiz obrigatória \(r\in T\), prêmios não negativos \(p_v\), quota total \(Q\), e demandas \(w\in W\). Cada demanda tem conjunto de cobertores \(N_w\subseteq V\). Uma solução é um ciclo simples com pelo menos três vértices, contendo \(T\), satisfazendo \(\sum_v p_vy_v\ge Q\) e \(\sum_{v\in N_w}y_v\ge1\). A variável \(x_e\) indica aresta do circuito e \(y_v\) indica visita, com \(y_t=1\) e \(x(\delta(v))=2y_v\). Quando só opcionais aparecem na expressão de prêmio, a quota deve ser residual, após subtrair os obrigatórios.

O artigo IJISE já contém conectividade generalizada, cortes forçados por cobertura e prêmio, covers e covers estendidas de prêmio, além de 2-matching. A auditoria não encontrou invalidade nas duas fórmulas examinadas: a cover estendida é válida por complementaridade \(u=1-y\); usar \(|T|\), em vez de \(|H\cap T|\), na constante da 2-matching apenas a enfraquece. A comparação abaixo adota a versão forte correta, evitando obter uma falsa contribuição por corrigir uma constante fraca. O callback encontrado em `hib.cpp` permanece comentado nesse fluxo, portanto a existência de um separador no repositório não demonstra sua ativação. ([Artigo PCCTP](https://www.inderscience.com/info/inarticle.php?artid=133526), [Referência sobre blossom e b-matching](https://www.lancaster.ac.uk/staff/letchfoa/articles/2008-bmatching.pdf)).

Há quatro perguntas diferentes. Validade exige prova para toda solução inteira; fortalecimento exige um ponto aceito pela relaxação anterior e rejeitado pela candidata; faceta exige dimensão e independência afim suficientes; novidade exige confronto com famílias anteriores e suas projeções. Nesta etapa temos prova de validade, testemunhas de fortalecimento relativo e certificados de faceta em duas instâncias de nove vértices. **Não provamos condições gerais de faceta nem originalidade da família.**

Um resultado negativo orienta a busca. Se \(O=V\setminus S\) contém todos os obrigatórios, tem pelo menos três vértices e não há restrições superiores de visitas, custo ou capacidade, então existe solução inteiramente em \(O\) exatamente quando \(O\) cobre todas as demandas e \(p(O)\ge Q\). Basta selecionar todos os vértices de \(O\): monotonicidade preserva cobertura e prêmio, e completude permite ordená-los em ciclo. Logo, combinar esses dois testes não descobre uma razão adicional para \(x(\delta(S))\ge2\) quando ambos já passam. Essa conclusão deixa de valer com esparsidade, limites superiores ou outras restrições. Também é falso deduzir \(x(\delta(S))\ge2k\) da necessidade de visitar \(k\) vértices de \(S\): o circuito pode entrar uma vez, visitar todos e sair uma vez.

## O posto de uma partição gera um corte válido

Seja \(\mathcal P=(P_0,P_1,\ldots,P_k)\) uma partição de \(V\), com \(r\in P_0\), e seja \(E_{\mathcal P}\) o conjunto de arestas entre partes diferentes. Defina \(\rho(\mathcal P)\) como o menor número de blocos não raiz ativados por uma seleção que satisfaça obrigatórios, cobertura e prêmio. A definição considera seleção monotônica, sem exigir que ela forme circuito; isso fornece uma obrigação válida, eventualmente mais fraca que a do problema de roteamento. Então, para \(\rho(\mathcal P)\ge1\), vale

\[
\boxed{x(E_{\mathcal P})\ge \rho(\mathcal P)+1.}
\]

**Prova.** Contraia cada bloco visitado do circuito em um vértice, apague laços e preserve multiplicidades. Se \(b\ge2\) blocos foram visitados, o multigrafo resultante é conexo e euleriano. Cada vértice tem grau par pelo menos dois, logo seu número de arestas é pelo menos \(b\). O bloco da raiz é visitado, portanto \(b\ge\rho+1\), e essas arestas são precisamente os cruzamentos da partição. Reentradas em blocos não invalidam o argumento, pois são preservadas como multiplicidades. Se \(\rho=0\), a desigualdade universal segura tem lado direito zero: um circuito pode permanecer em \(P_0\).

Para uma partição fixa, podemos calcular o posto com variáveis binárias \(z_i\) de ativação dos blocos não raiz. Blocos contendo obrigatórios têm \(z_i=1\). Para cada demanda não coberta por \(P_0\), impomos \(\sum_{i:P_i\cap N_w\ne\varnothing}z_i\ge1\). Acrescentamos \(p(P_0)+\sum_i p(P_i)z_i\ge Q\), minimizando \(\sum_i z_i\). Disponibilizar todos os vértices de um bloco ativado torna esse oráculo exato para a seleção monotônica. Um **limite inferior certificado** do oráculo também fornece lado direito válido; o valor de uma solução heurística é limite superior e não serve. Arredondar para cima um limite inferior LP é permitido pela integralidade do posto.

Isso ainda não é um separador geral eficiente: mesmo o oráculo contém cobertura mínima de conjuntos como caso particular. Encontrar a melhor partição é um problema adicional, cuja complexidade não foi caracterizada aqui. Uma implementação inicial pode contrair arestas de alto valor fracionário, gerar partições pequenas e resolver oráculos de quatro a dez blocos; só se emitem cortes com certificado. Essa é uma estratégia experimental, não uma conclusão teórica de polinomialidade.

Quando os blocos são singletons, a expressão pode reduzir-se a cardinalidade via equações de grau. Quando todos os blocos necessários são individualmente obrigatórios, somar seus cutsets pode bastar. A oportunidade está em **blocos não triviais cuja obrigação é coletiva**, com alternativas de visita e cobertura sobreposta.

## Uma testemunha separa roteamento de seleção

A testemunha persistente tem nove vértices: raiz \(0\) e blocos \(B_1=\{1,2\}\), \(B_2=\{3,4\}\), \(B_3=\{5,6\}\), \(B_4=\{7,8\}\). Apenas a raiz é obrigatória. Existem seis demandas, com cobertores \(B_i\cup B_j\), uma para cada par de blocos. Todo opcional tem prêmio um, a raiz tem prêmio zero e a quota é um. Se dois blocos ficassem vazios, a demanda correspondente não seria coberta; logo \(\rho=3\) e o corte exige quatro cruzamentos.

Considere o ponto racional

\[
y_0=1,\quad y_v=5/8\ (v\ne0),\qquad
x_{0v}=1/4,\quad x_{uv}=5/8\text{ dentro de }B_i,
\quad x_{uv}=1/16\text{ entre blocos opcionais}.
\]

Os graus e ligações \(x_{uv}\le\min(y_u,y_v)\) são satisfeitos. Entretanto,

\[
x(E_{\mathcal P})=8(1/4)+24(1/16)=7/2<4.
\]

Mais forte: **o vetor \(y\) já pertence ao casco convexo completo das seleções viáveis dessa instância**. Escolha uniformemente três dos quatro blocos; em cada um deles, selecione os dois vértices com probabilidade \(2/3\), só o primeiro com \(1/6\), ou só o segundo com \(1/6\). Cada seleção cobre as seis demandas e satisfaz quota. A marginal de cada opcional é \((3/4)(5/6)=5/8\). A decomposição tem 108 termos e pesos racionais. Assim, nenhuma desigualdade exclusivamente em \(y\) elimina o ponto. O déficit é estruturalmente de roteamento entre grupos.

O verificador testa 2.048 instâncias das SEC generalizadas, 134 cortes forçados por cobertura, dois por prêmio e a cover de prêmio aplicável. Uma programação dinâmica exata de paridade certifica todas as 154.732.800 escolhas ímpares de \(F\subseteq\delta(H)\) contabilizadas pelo script, sem enumerá-las individualmente, usando a forma forte

\[
x(\delta(H))\ge 2x(F)-|F|+1.
\]

O teste permite inclusive \(F\) sem a restrição de matching, portanto cobre a versão clássica com dentes disjuntos. Foram enumeradas 189 seleções e 54.624 circuitos viáveis distintos até reversão; todos satisfazem o corte. O programa também aceita 2.254 expressões obtidas pela transposição explícita de cover-intersection registrada no script. Isso não deve ser ampliado para uma alegação sobre toda a família do covering salesman: as hipóteses da fórmula impressa e sua versão publicada ainda exigem auditoria. ([Covering salesman, preprint dos autores](https://arxiv.org/abs/2104.01173), [Versão publicada](https://doi.org/10.1051/ro/2023055)).

A testemunha admite custos métricos: custo \(1/10\) dentro de cada bloco, um entre blocos, distância de uma demanda a seus cobertores igual a um e aos demais igual a dois, e dois entre demandas. Com raio um, as vizinhanças são exatamente as descritas. Foram verificadas as 3.375 desigualdades triangulares ordenadas nos 15 pontos. O custo fracionário é \(7/2+(5/2)(1/10)=15/4\); o ótimo inteiro enumerado vale quatro, pois visitar raiz e um vértice de cada um de três blocos atinge esse custo. O corte implica custo pelo menos quatro nesta instância. Isso estabelece diferença estrita numa métrica, **sem afirmar realização euclidiana plana**.

Essa instância prova força sobre a relaxação descrita, não dominância universal sobre cada família citada nem independência frente a todas as desigualdades GTSP/Steiner existentes. Sua quota tampouco é essencial. Por isso verificamos uma **segunda testemunha genuinamente conjunta**: retiramos a demanda de cobertores \(B_3\cup B_4\), atribuímos prêmio total um a cada um dos blocos \(B_1,B_2\) e dois a cada um de \(B_3,B_4\), distribuído igualmente entre seus vértices, com quota três. Cobertura isolada permite selecionar \(B_1\cup B_2\); prêmio isolado permite \(B_3\cup B_4\). Mas cobertura e prêmio juntos exigem três blocos: o único par que cobre as cinco demandas é \(B_1,B_2\), cujo prêmio total é apenas dois. Assim, os postos separados são dois e o conjunto é três.

O mesmo ponto fracionário passa pelos testes. Seu vetor \(y\) tem decomposição explícita em 24 seleções viáveis equiprováveis: escolher três blocos e omitir um dos seis vértices disponíveis. Foram verificados 2.048 SEC, 116 cortes de cobertura, 212 de prêmio, 106 covers simples, 106 estendidas e 106 knapsack covers residuais, além das mesmas 154.732.800 escolhas ímpares certificadas por programação dinâmica. Foram enumeradas 141 seleções e 54.336 circuitos viáveis. O custo fracionário continua \(15/4\), contra ótimo inteiro \(41/10\), e as 2.744 verificações triangulares passam. A diferença maior não deve ser atribuída inteiramente ao novo corte: ele impõe apenas custo interblocos de pelo menos quatro.

Há uma limitação conceitual importante: nesta configuração \(K_4\), o corte é obtido por **arredondamento Chvátal–Gomory de posto um** dos seis cortes forçados para uniões de dois blocos. Escrevendo \(a=x(\delta(\{0\}))=2\) e \(b\) para arestas entre blocos opcionais, a soma fornece \(3a+4b\ge12\). Da igualdade da raiz resulta \(a+b\ge7/2\), que se arredonda para quatro porque cruzamentos são inteiros. Na variante mista, cinco obrigações vêm da cobertura e a sexta do prêmio. Portanto a testemunha prova não implicação pela relaxação linear testada, mas **não** independência frente ao seu fecho CG de posto um. Uma contribuição precisaria caracterizar e separar eficientemente essa estrutura, ou demonstrar uma generalização não absorvida por famílias conhecidas.

Também examinamos a dimensão dos politopos nessas duas instâncias. Há 44 coordenadas livres, correspondentes a 36 arestas e oito visitas opcionais. As nove equações de grau independentes dão dimensão no máximo 35; circuitos viáveis com diferenças de posto 35 fornecem o limite inferior coincidente. Acrescentar a equação do corte dá dimensão facial no máximo 34, e circuitos justos fornecem posto 34. A independência foi certificada por eliminação modular no primo 1.000.003: um menor não nulo módulo primo é não nulo sobre os racionais. Isso certifica **faceta apenas nesses dois politopos concretos**. O script para ao atingir limites coincidentes, após 12.601 e 12.313 circuitos processados, respectivamente; não exige enumerar todos para provar dimensão. As rotas geradoras ficam no JSON. Não decorre daí que todo corte de posto de partição seja facetário.

## Obrigações de pares geram subclasses interpretáveis

Fixada a partição, construa um grafo cujos vértices são os blocos não raiz. Insira a aresta \(ij\) se existir demanda com \(N_w\subseteq P_i\cup P_j\), ou se \(p(V\setminus(P_i\cup P_j))<Q\). Em ambos os casos, ao menos um desses blocos precisa ser visitado. Obrigações unárias decorrem de obrigatórios, cobertores inteiramente num bloco ou complemento de prêmio insuficiente. Os blocos ativos formam uma cobertura de vértices desse grafo, respeitando as obrigações unárias.

Uma clique de tamanho \(q\ge2\) exige \(q-1\) blocos ativos, fornecendo \(x(E_{\mathcal P})\ge q\). Um ciclo ímpar de comprimento \(2h+1\) exige \(h+1\) blocos, fornecendo \(x(E_{\mathcal P})\ge h+2\). A prova é a cardinalidade mínima de vertex cover: na clique omite-se no máximo um vértice; no ciclo, somam-se as obrigações de pares e arredonda-se. A testemunha anterior corresponde a uma clique \(K_4\). Cliques pequenas podem ser enumeradas para tamanho fixo; clique de tamanho arbitrário não se torna fácil. Em grafos bipartidos, matching fornece a cobertura mínima. Esses são instrumentos para obter postos certificados, sem afirmar descoberta de clique ou odd-cycle inequalities.

Há também fortalecimento conjunto puramente em seleção, útil como controle. Considere opcionais \(a,b,c\), prêmios \((1,1,2)\), quota dois e demanda coberta por \(\{a,b\}\). Toda seleção viável usa pelo menos dois opcionais, logo \(y_a+y_b+y_c\ge2\). O ponto \((1/2,1/2,1/2)\) passa quota, cobertura, todas as knapsack covers e covers estendidas dessa instância, mas viola a cardinalidade. Ela é uma consequência **Chvátal–Gomory de posto um conhecida**: cobertura fornece \(y_a+y_b\ge1\), quota fornece \(y_a+y_c\ge1\) e \(y_b+y_c\ge1\); somar, dividir por dois e arredondar produz o resultado. Não constitui família genérica nova. ([Knapsack-cover e antecedentes](https://theoryofcomputing.org/articles/v014a014/)).

O segundo verificador estende esse ponto a quatro terminais obrigatórios: \(x=1/4\) em terminal–opcional, \(x=5/12\) entre terminais e zero entre opcionais. Executou 749 verificações racionais, incluindo 63 cutsets sem raiz, três knapsack covers, três covers estendidas e 672 máximos blossom por handle/cardinalidade ímpar; enumerou 540 ciclos viáveis. É uma testemunha relativa rigorosa, mas explica por que a proposta principal precisa ultrapassar apenas rank em \(y\).

## A contribuição exige confronto poliédrico e teste real

A prioridade bibliográfica é comparar o posto de partição com Steiner F-partition, generalized subtour e cortes de GTSP/group TSP com grupos sobrepostos. Set-cover rank é nomenclatura estabelecida; escrever um ótimo auxiliar no lado direito não prova uma família nova. GTSP já tem estudo poliédrico e branch-and-cut com famílias de facetas, e Steiner 2-edge-connected já oferece desigualdades de partição. Precisamos confrontar hipóteses, coeficientes e projeções, não apenas nomes. ([Rank de set cover](https://www.zib.de/userpage/groetschel/students/borndoerfer.thesis.pdf), [GTSP](https://doi.org/10.1287/opre.45.3.378), [Steiner F-partition](https://www.numdam.org/item/10.1051/ro:2008022.pdf)).

Uma comb ativada por conjuntos obrigatoriamente visitados também possui prova válida: para handle e número ímpar de dentes disjuntos, exigir visita nas duas metades de cada dente permite obter \(x(\delta(H))+\sum_i x(\delta(T_i))\ge3t+1\) por conectividade e paridade. Porém combs e liftings são antecedentes clássicos; o exemplo mínimo já cai em 2-matching, e não obtivemos testemunha estrita contra todas as blossom para dentes maiores. Sua prioridade é menor. ([Comb clássico](https://www.zib.de/userpage/groetschel/pubnew/paper/groetschelpadberg1978.pdf), [Generalized comb com opcionais](https://arxiv.org/abs/1508.01813)).

O próximo passo é gerar politopos pequenos não simétricos, enumerar facetas, remover equivalências por graus e reconhecer famílias conhecidas. Em paralelo, devemos medir violações do posto de partição nas relaxações das instâncias reais antes de desenvolver um separador complexo. Só depois vale comparar limite na raiz, tempo de separação, nós, gap global, tempo de prova e quantidade resolvida, com e sem cortes, incluindo baselines corrigidos e knapsack-cover residual. Nenhuma melhoria nesses indicadores foi testada nesta etapa.

A reprodução usa apenas Python padrão, a partir da raiz do repositório:

```sh
python3 experiments/pcctp_cuts/verify_rank_witness.py
python3 experiments/pcctp_cuts/verify_joint_witness.py
python3 experiments/pcctp_cuts/verify_mixed_partition_witness.py
python3 experiments/pcctp_cuts/verify_partition_face.py
```

Os certificados persistentes são `experiments/pcctp_cuts/rank_witness_result.json` e `experiments/pcctp_cuts/joint_witness_results.json`, além de `mixed_rank_witness_result.json` e `partition_face_results.json` na mesma pasta. A aritmética exata certifica o ponto específico e a enumeração finita; a validade para tamanho arbitrário vem das provas. Mineração pode sugerir partições ou suportes promissores, mas padrões observados nunca substituem prova de validade. O resultado atual justifica uma investigação teórica dirigida: **obrigações coletivas entre blocos geram força que seleção isolada não recupera**. A eventual publicação dependerá de estabelecer o que há de novo nessa estrutura e de demonstrar utilidade além da testemunha.
