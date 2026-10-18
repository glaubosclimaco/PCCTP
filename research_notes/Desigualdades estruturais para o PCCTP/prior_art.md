# Antecedentes poliédricos e filtro de novidade para o PCCTP

## Quais famílias já existem e o que está no IJISE?

### Takeaway
Subtours generalizados, knapsack covers, 2-matchings, combs com lifting e partições Steiner são famílias antigas. Uma especialização ao PCCTP pode ser útil, mas não deve ser chamada nova sem comparar fórmulas, hipóteses e projeções.

### Cited Findings
- **Fonte primária local:** `artigos/IJISE/sample (single col).tex`, linhas 229–395. Já contém: grau `x(δ(i))=2y_i`; obrigatório `y_t=1`; GSEC `x(δ(S))≥2y_i`; cobertura forçada `x(δ(S))≥2` se algum `R_w⊆S`; cover de quota `Σ_{i∈C}y_i≥1` quando `p(R\C)<Q`; extensão com itens exteriores de prêmio alto; 2matching; corte forçado por prêmio exterior insuficiente e argumento de incumbent. [Publicação](https://www.inderscience.com/info/inarticle.php?artid=133526).
- **CTP / strong 2matching:** Erdoğan, Cordeau e Laporte, The Attractive Traveling Salesman Problem, CIRRELT-2007-33, proposições 3–5, atribuem ao CTP de Gendreau–Laporte–Semet1997 as desigualdades `x_ij≤y_i,y_j`, cortes de conectividade e `x(E(H))+x(F)≤Σ_{i∈H}y_i+(|F|−1)/2`, onde F é conjunto ímpar ≥3 de arestas cruzando H, duas a duas sem vértices comuns. [PDF dos autores, pp9–10](https://www.cirrelt.ca/documentstravail/cirrelt-2007-33.pdf). Ao fixar obrigatórios, a constante correta é `|H∩T|`, não `|T|` como escrito no manuscrito local; a versão local é apenas mais fraca.
- **Prize-collecting TSP:** Balas1995, The prize collecting traveling salesman problem: II. Polyhedral results, DOI10.1002/net.3230250406, já apresenta transferência/lifting de famílias comb, odd CAT, SD, clique-tree e lifted-cycle do ATSP. [Editora](https://onlinelibrary.wiley.com/doi/abs/10.1002/net.3230250406). Texto completo das fórmulas não ficou acessível nesta busca; não usar coeficientes inventados.
- **GTSP:** Fischetti, Salazar González e Toth1997, Operations Research45(3):378–394, DOI10.1287/opre.45.3.378, fornece branch-and-cut com separação de famílias de facetas; remete ao estudo poliédrico Networks1995. Hipótese central: clusters disjuntos. [Editora](https://pubsonline.informs.org/doi/abs/10.1287/opre.45.3.378).
- **Generalized comb explícito:** Sundar–Rathinam, Generalized multiple depot traveling salesmen problem—polyhedral study and exact algorithm, preprint2015, prop3.6 (adaptada do GTSP1995). Para handle H e t dentes disjuntos, t≥3 ímpar, cada dente interceptando H e seu exterior, `x(E(H))+Σ_jx(E(T_j))+Σ_v β_v(1−y_v)≤|H|+Σ_j|T_j|−(3t+1)/2`. Coeficientes dependem de cada região conter cluster completo; sem cluster, escolhem-se representantes a_j∈H∩T_j e b_j∈T_j\H: β=1 no handle fora dos dentes, β=2 na interseção salvo β_a=1, β=1 fora do handle salvo β_b=0; zero fora da comb. [PDF dos autores, prop3.6, p11](https://arxiv.org/pdf/1508.01813). Transferência de condições de cluster ao PCCTP exige prova, e facetness não se transfere automaticamente.
- **Knapsack-cover residual:** para `Σp_i y_i≥Q`, A com `p(A)<Q`, U=Q−p(A)>0, vale `Σ_{i∉A}min(p_i,U)y_i≥U`. Esta é a família clássica de Carr–Fleischer–Leung–Phillips SODA2000, não uma descoberta PCCTP. Fórmula/atribuição e discussão de separação no artigo de Bazzi–Fiorini–Huang–Svensson, Theory of Computing2018, eq1.1. [PDF primário](https://theoryofcomputing.org/articles/v014a014/v014a014.pdf). É diferente do cover não ponderado no IJISE e vale testar como baseline mais forte.

### Inferences
- Usar comb de TSP com RHS constante sem corrigir ausência de vértices opcionais pode ser inválido. Usar família conhecida com y/lifting é baseline, não reivindicação de nova família.
- O prêmio Q deve estar normalizado: se variáveis existem só em R, Q é quota residual após obrigatórios; não misturar com PRIZE total.
- Não tratar corte condicionado ao custo de incumbent como desigualdade do politopo de TODAS as soluções viáveis: ele pode preservar apenas soluções melhores/ótimas sob limite objetivo declarado.

### Gaps
- Nenhuma fórmula multistar específica ao PCCTP foi confirmada; multistars de VRP dependem de capacidade/demandas e não podem ser transportadas só pela semelhança do desenho.
- Falta verificar texto completo Balas1995 e Fischetti1995 antes de reivindicar lifting novo.

## O que rank, partições e interseção de cobertura já cobrem?

### Takeaway
Um RHS definido pelo ótimo de um problema auxiliar de seleção é um rank por definição. Sua projeção sobre partições pode gerar cortes úteis, mas requer estudo de dominância e antecedentes Steiner para se tornar contribuição estrutural.

### Cited Findings
- **Set-cover rank:** para matriz binária A, `β(A)=min{1ᵀy:Ay≥1,y∈{0,1}ⁿ}`, a desigualdade `Σy_i≥β(A)` é chamada rank inequality. Borndörfer explicita definição e minors, closedness/nonseparability na tese, seção set covering, página impressa48 (PDF65). [PDF institucional](https://www.zib.de/userpage/groetschel/students/borndoerfer.thesis.pdf). Logo uma cardinalidade mínima de seleção ou sua reformulação via graus é antecedente conhecido.
- **Steiner 2-edge-connected:** Mahjoub–Pesneau2008, RAIRO42:259–283, DOI10.1051/ro:2008022, contém Steiner F-partition e generalizações com coeficientes0,1,2. Exemplo explícito válido na wheel com terminais u1,u2,u3: `2x(e1)+2x(e2)+x(e3)+x(e4)+x(f4)≥4`. Em partições específicas de wheel, família simplificada `x(Δ1)+2x(Δ2\F)≥q+1` com q ímpar, Δ1/Δ2 construídos segundo regiões sem/com terminais, F arestas do centro a regiões com terminais. [PDF primário, eq2.2–2.4](https://www.numdam.org/item/10.1051/ro:2008022.pdf). Não transportar coeficientes sem todas as condições topológicas.
- **Partições com ordem terminal2025:** The multiple steiner TSP with cyclic order on terminals: valid inequalities and polyhedra, DOI10.1007/s10878-025-01288-1, prop5 usa partição (V1,…,Vp) com cada parte contendo terminal do veículo; algumas partes contêm q_i≥2 terminais/sequências não sucessivos. Para S de nós Steiner, `x^k(δ_{G\S}(V1,…,Vp))≥p+Σ_{i=1}^r q_i−r−|S|`. Caso r=0 dá p−|S|. [Texto primário](https://link.springer.com/article/10.1007/s10878-025-01288-1). A ordem cíclica fixa é essencial ao fortalecimento q_i; PCCTP livre não a tem.
- **Cover intersection:** Maziero–Usberti–Cavellucci, preprint2021 / RAIRO2023 57:1149–1166 DOI10.1051/ro/2023055. Usa GSEC base `x(δ(S))≥2(y_i+y_j−1)`, i∈S,j∉S; especializações com conjuntos de cobertura forçando presença; e CI envolvendo união de dois cortes: com S_v=S∩C(v), preprint eq10 imprime `x(δ(S)∪δ(S_v))≥2` para S contendo um conjunto C(u). [Preprint primário, seção3](https://arxiv.org/pdf/2104.01173). **Cautela:** as hipóteses de validade dessa forma impressa exigem auditoria: um ciclo todo em S_v daria LHS zero. Não implemente a forma geral só pela citação; conferir versão publicada/condição ausente e contraexemplos pequenos. [Publicação](https://doi.org/10.1051/ro/2023055).

### Inferences
- Para partição P0 contendo raiz e P1,…,Pk, seja ρ mínimo número de blocos não raiz ocupados por seleção que atende T, cobertura e quota. Se ρ≥1, qualquer ciclo factível usa `x(E_entre_partes)≥ρ+1`: após contração, multigrafo conexo Euleriano tem ≥número de vértices ocupados. Se ρ=0, não impor RHS1, pois ciclo pode ficar dentro de P0. **Candidato de especialização rank/partition, novidade não estabelecida.**
- Se todos os blocos necessários são individualmente forçados, a desigualdade é metade da soma de seus cutsets; pouca perspectiva de novidade. Se singleton, equivale a rank de cardinalidade via grau. A oportunidade está em prova de não dominância para partições não triviais com obrigação coletiva, lifting e separação eficiente.
- Definir função auxiliar ρ por MIP e utilizá-la como RHS é válido, mas não produz por si uma família poliédrica explicitamente caracterizada nem separação escalável.
- Para produzir rigor: procurar fração satisfazendo todas famílias antigas e violando candidata; provar que não é combinação não negativa das antigas; testar dimensão/facetness nas instâncias pequenas; confrontar casos degenerados (sem quota, sem cobertura, todos obrigatórios).

### Gaps
- A condição precisa do CI publicado e seu efeito com raiz obrigatória continuam pendentes. No PCCTP, se raiz está fora de S e S já é forçado pela cobertura, a união dos cortes é dominada por `x(δ(S))≥2` e não fortalece nada.
- Não foi demonstrada novidade da função ρ quota-cobertura sobre quociente frente a generalized Steiner partition / group Steiner / GTSP.

## Quais propostas devem ser rejeitadas cedo e o que merece descoberta?

### Takeaway
RHS2 forçado por prêmio exterior insuficiente já existe; RHS2k por exigir k visitas dentro de um conjunto não vale em geral. É mais promissor descobrir coeficientes por facetas pequenas e depois interpretar o padrão combinatório do que nomear uma fórmula previsível como nova.

### Cited Findings
- O próprio IJISE seção cut-sets já cobre a impossibilidade de atender quota fora de S. [Artigo](https://www.inderscience.com/info/inarticle.php?artid=133526).
- Rank e lifting geral já têm teoria extensa, conforme fontes set-cover, GTSP e PCTSP acima. [Rank](https://www.zib.de/userpage/groetschel/students/borndoerfer.thesis.pdf), [Balas](https://onlinelibrary.wiley.com/doi/abs/10.1002/net.3230250406).

### Inferences
- **No-go monotônico:** em grafo completo sem upper bound em visitas, existe seleção factível completamente fora de S se e somente se V\S contém T, cobre todos W e seu prêmio total é ≥Q, além de permitir ciclo de pelo menos3 vértices. Selecionar TODOS fora de S simultaneamente maximiza cobertura e prêmio. Portanto uma suposta sinergia quota+cobertura para provar apenas entrada em S não acrescenta critério além dos já conhecidos.
- **Contraexemplo conceitual a 2k:** se k visitas são exigidas dentro de S, um ciclo pode entrar uma vez, percorrer um caminho pelos k vértices e sair uma vez. Exige-se apenas2 cruzamentos. Maior RHS precisa de desconexão interna, arestas proibidas, múltiplos componentes ou paridade coletiva demonstrada.
- **Condicionais padrão:** se desligar U força visita a S, a desigualdade `x(δ(S)) + 2Σ_{i∈U}y_i≥2` é válida com raiz fora de S, mas é disjuntiva elementar e pode ser dominada por GSEC+cover. Novidade requer coeficientes levantados e testemunha fracionária.
- **Programa de descoberta recomendado:** enumerar ciclos viáveis de modelos PCCTP pequenos não simétricos, computar facetas, remover equivalências por graus e permutações, reconhecer covers/comb/partition/KCI, estudar facetas remanescentes que misturem x e y e dependam simultaneamente de quota e cobertura. Uma família candidata deve vir acompanhada de prova, exemplar de separação, não dominância e caminho para separador prático.
- **Mineração como ferramenta:** pode sugerir suportes de cortes ou blocos para separação heurística, mas padrões observados nunca provam validade; separar geração de hipótese e certificação matemática.

### Gaps
- Nenhuma família nova certificada é reivindicada neste memorando. Falta enumeração poliédrica e confronto completo das referências próximas.
- A coordenação encontrou testemunha adicional para partição com raiz+4 blocos de2 vértices e seis demandas, cada uma coberta pela união de um par de blocos: rank de ocupação não raiz ρ=3, RHS4, fração com custo interblocos3,5. A fração de y pertence ao hull de seleções viáveis (escolher3 blocos uniformemente; em cada um, escolher ambos com probabilidade2/3 e cada singleton com1/6), e a coordenação verificou passagem pelas antigas famílias. Isto é evidência de ganho de roteamento, não apenas rank em y; estes cálculos não foram refeitos independentemente neste memorando. Com quota não ativa, a testemunha prova fenômeno de cobertura, não interação essencial entre quota e cobertura.
- Busca focada por `covering salesman/covering tour + partition inequalities` não encontrou fórmula exata equivalente a ρ+1. Ausência de resultado de busca não prova novidade. Os antecedentes prioritários restantes são generalized Steiner F-partition, rank em quocientes de set-cover/group-Steiner, e GTSP com grupos sobrepostos.
