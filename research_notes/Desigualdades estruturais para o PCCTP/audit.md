# Auditoria independente de desigualdades para o PCCTP

## 1. As suspeitas de invalidade do manuscrito procedem?

**Conclusão: não, nas duas fórmulas examinadas.** A cover estendida é válida; o termo constante da 2-matching torna a desigualdade mais fraca, não inválida. Não há contraexemplo de ciclo factível sob as definições corretas. A ausência de contraexemplo aqui decorre das provas abaixo, não de buscas inconclusivas.

Notação: Q é prêmio residual exigido dos opcionais R, depois de subtrair prêmios obrigatórios. Se o manuscrito usa PRIZE para a quota original, sua linha de prêmio em R precisa desse esclarecimento. Para evitar ambiguidade, z_i=1 em T e z_i=y_i em R.

### Cover estendida: prova

O manuscrito define C⊆R com p(R\C)<Q e B={i∈R\C:p_i≥max_{j∈C}p_j}. A desigualdade é y(C)+y(B)≥1+|B|.

Defina u_i=1−y_i e U=p(R)−Q. A quota implica p·u≤U, e p(C)>U. Qualquer conjunto de |C| itens em E=C∪B tem prêmio total pelo menos p(C), pois todo item de B é tão pesado quanto o maior de C. Portanto não se podem omitir |C| itens de E: u(E)≤|C|−1. Substituir u=1−y produz exatamente a fórmula publicada. A minimalidade de C afeta força e separação, não é necessária nesta prova.

### 2-matching: prova por paridade

Para H⊆V e F⊆δ(H), |F| ímpar, graus dão 2z(H)=2x(E(H))+x(F)+x(δ(H)\F). Como x é binário, x(F)−x(δ(H)\F) é inteiro par e no máximo |F|, logo no máximo |F|−1. Assim:

    x(E(H))+x(F) ≤ |H∩T|+y(H\T)+floor(|F|/2).

O manuscrito usa |T| em vez de |H∩T|; isso acrescenta |T\H| à RHS, preservando validade e perdendo força. Teeth sem extremidades comuns não são necessárias para esta prova de validade da desigualdade blossom; são relevantes para a apresentação clássica e propriedades de faceta. A prova exige F realmente contido no corte, sem contar a mesma aresta várias vezes, e |F| ímpar.

Fonte primária para o contexto de blossom e separação: [Letchford, Reinelt e Theis, Odd Minimum Cut Sets and b-Matchings Revisited](https://www.lancaster.ac.uk/staff/letchfoa/articles/2008-bmatching.pdf). Fórmulas auditadas: [manuscrito local](/Users/glaubosclimaco/Projects/papers/PCCTP/artigos/IJISE/sample%20(single%20col).tex:275).

### O código implementa a fórmula impressa?

- `include/Mycallback.h`, função `cover_cut`, linhas 995–1148: trabalha explicitamente com 1−y e capacidade `prize_V-prize`; constrói cover por p(C)>U e adiciona itens com prêmio≥máximo de C. É a transformação correta **se** os dois totais estiverem contabilizados consistentemente após reduções. Essa consistência não foi auditada integralmente.
- `kruskal2Matching`, linhas 479–525: RHS usa `handle.nos.size()+floor(nTeeths/2)`, isto é, |H|; não contém as variáveis opcionais y e não é a fórmula adaptada impressa. O corte clássico é válido porém pode ser mais fraco que a forma com z(H).
- `src/hib.cpp`, linhas 882–886: criação e registro do callback estão comentados. Portanto encontrar código de separação não significa que esteja ativo no executável deste fluxo.
- Esta auditoria não valida todas as decisões do separador nem substitui testes de candidatos de cortes contra soluções inteiras. Não editar solver nesta etapa.

## 2. Quais candidatos devem ser rejeitados como novos?

### KC condicional é dominada por KC convencional

Para A⊆R com p(A)<Q, r=Q−p(A)>0, a desigualdade conhecida

    sum_{i∉A} min(p_i,r)y_i ≥ r

vale globalmente, mesmo que A não esteja selecionado. De fato, p(R\A selecionado)≥Q−p(A)=r; para seleção binária, truncar individualmente por r preserva soma≥r. A variante com RHS r[1−sum_{i∈A}(1−y_i)] é mais fraca porque o fator é≤1 para todo y∈[0,1]. Não deve ser proposta como nova família. Ver [Bazzi et al., 2018](https://theoryofcomputing.org/articles/v014a014/).

### Certificado conjunto de inviabilidade do complemento: resultado de impossibilidade

Suponha grafo completo, sem orçamento de custo ou comprimento, e O=V\S contendo todos os obrigatórios, |O|≥3. Então O contém uma solução factível **se e somente se** p(O)≥quota original e O cobre todas as demandas. Necessidade é imediata. Para suficiência, selecionar TODOS os vértices de O preserva prêmio/cobertura por monotonicidade; qualquer ordem cíclica é um ciclo simples existente no grafo completo.

Consequência: nenhum teste puramente conjunto de cobertura+prêmio torna mais forte o certificado “não pode existir solução só no complemento” quando ambos testes individuais passam. Portanto não existe nova razão estrutural para δ(S)≥2 nesse cenário apenas por combinar as duas condições. Falha quando |O|<3, o grafo é esparso ou existem orçamento, capacidades ou outras restrições. Limites de custo relativos a incumbent podem preservar ótimo mas não são desigualdades válidas de todo o politopo factível sem corte de objetivo explícito.

## 3. Onde há fortalecimento real e como demonstrá-lo?

### Rank conjunto de seleção

Para U⊆R, permita todos os vértices fora de U. Subtraia seus prêmios da quota, e remova demandas que eles cobrem. Defina rho(U) como o mínimo número de vértices de U capaz de satisfazer simultaneamente quota residual e demandas restantes. Então y(U)≥rho(U) é válido: a parte em U de qualquer solução factível serve ao problema residual, pois disponibilizar todos os externos só facilita suas restrições. Um limite INFERIOR comprovado sobre rho basta; um valor de uma solução heurística do subproblema é limite superior e NÃO pode ser colocado na RHS.

Calcular rho exatamente é NP-difícil: com quota zero, resolve cobertura mínima de conjuntos. Isso descarta a expectativa de oráculo geral polinomial, mas não separadores parciais ou subclasses tratáveis. O nome rank não implica novidade: trata-se de consequências válidas do politopo conjunto de seleção.

### Testemunha estrita, racional e reproduzível

R={a,b,c}, prêmios (1,1,2), Q=2; uma demanda tem cobertores {a,b}. Toda seleção factível usa≥2 opcionais, logo y_a+y_b+y_c≥2. Entretanto y_a=y_b=y_c=1/2 passa quota, cobertura e **todas** as KC e covers estendidas de prêmio.

Extensão completa ao roteamento: quatro terminais obrigatórios de prêmio zero; x=1/4 em cada aresta terminal–opcional, x=5/12 em cada par de terminais, zero entre opcionais. Esse ponto passa os graus, contagem de arestas, ligações x≤y, todos os 63 cutsets sem raiz e todas as 2-matchings com RHS forte z(H)+floor(|F|/2). A rank conjunta corta com violação 1/2. A família de vizinhanças de cobertura é realizável por uma métrica com distâncias 1 ou 2 (distância 1 aos dois cobertores, 2 aos demais), com D=1.

Importante: a rank é uma consequência Chvátal–Gomory de posto 1 das três obrigações de pares:

    ya+yb≥1  (cobertura)
    ya+yc≥1  (quota/KC)
    yb+yc≥1  (quota/KC)

Somar, dividir por 2 e arredondar dá ya+yb+yc≥2. Portanto há fortalecimento genuíno sobre a interseção das famílias separadas testadas, mas **não uma nova família genérica**. A existência de uma testemunha não prova novidade no PCCTP nem força sobre o fecho CG das restrições existentes.

### Subclasse útil: ciclos ímpares de obrigações de pares

Construa grafo em R com aresta ij quando y_i+y_j≥1 é obrigação válida, por cobertura ou por quota (p(R\{i,j})<Q). Todo ciclo ímpar C gera y(C)≥(|C|+1)/2. Isso é a desigualdade conhecida de ciclo ímpar do politopo de vertex cover, agora aplicada ao grafo misto de obrigações. Separação polinomial: com pares já satisfeitos, usar pesos não negativos w_ij=y_i+y_j−1 e buscar circuito ímpar de peso<1 em grafo expandido por paridade. Um passeio ímpar contém ciclo ímpar; pesos não negativos permitem extrair um ciclo violado. Esta é uma adaptação justificável e implementável, mas a possível contribuição reside na caracterização específica, separação eficiente no PCCTP e evidências, não na alegação de descobrir odd-cycle cuts.

### Artefatos executados

- [Verificador](/Users/glaubosclimaco/Projects/papers/PCCTP/experiments/pcctp_cuts/verify_joint_witness.py)
- [Resultados JSON](/Users/glaubosclimaco/Projects/papers/PCCTP/experiments/pcctp_cuts/joint_witness_results.json)
- Execução: `python3 experiments/pcctp_cuts/verify_joint_witness.py`.
- Resultado: PASS; 749 verificações exatas; 63 cutsets; 3 KC; 3 covers estendidas; 672 máximos blossom por handle/cardinalidade ímpar; 540 ciclos inteiros factíveis enumerados. Para cada cardinalidade de teeth, escolher as maiores componentes x maximiza a violação, verificando implicitamente todos os F possíveis. Não depende de Gurobi, tolerâncias flutuantes nem licença.

O script não testa separador real nem desempenho. Sua função é estabelecer uma testemunha de não redundância relativa às famílias explícitas auditadas e validar a rank em todos os ciclos da instância diminuta.

### Auditoria cruzada adicional da testemunha por blocos

Foi lido independentemente `experiments/pcctp_cuts/verify_rank_witness.py`, de outro pesquisador: o DP de paridade maximiza corretamente sum(2x−1) sobre todos F ímpares, e a identidade de graus torna seu teste equivalente à blossom forte. Quatro blocos opcionais cujas uniões aos pares são demandas exigem pelo menos três blocos opcionais visitados; com o bloco obrigatório, qualquer ciclo cruza a partição pelo menos quatro vezes.

Também executei independentemente uma extensão com três terminais e oito opcionais (11 vértices): bloco obrigatório K3 com x=2/3 nas suas três arestas, x=1/12 nas arestas terminal–opcional, quatro blocos opcionais de tamanho dois com x interno 5/8 e x entre blocos 1/16; y opcional=5/8. Resultado racional exato: graus corretos, todos os 2046 cortes não vazios próprios passam SEC/cobertura/obrigatórios, e todos os máximos blossom por DP passam; massa das arestas entre blocos=7/2<4. Isso confirma que a testemunha não depende de haver apenas um terminal. Essa checagem extra foi executada em comando Python autônomo; o artefato persistente principal de blocos é o script do outro pesquisador.

### Certificado exato de faceta em duas instâncias finitas

Criados [verify_partition_face.py](/Users/glaubosclimaco/Projects/papers/PCCTP/experiments/pcctp_cuts/verify_partition_face.py) e [partition_face_results.json](/Users/glaubosclimaco/Projects/papers/PCCTP/experiments/pcctp_cuts/partition_face_results.json). São 44 coordenadas (36 arestas e 8 y opcionais), com nove equações de grau independentes, portanto dimensão≤35; adicionar igualdade da partição fornece dez equações independentes, portanto face≤34.

O programa encontra rank afim 35 no politopo inteiro e rank afim 34 nos ciclos que cruzam a partição exatamente quatro vezes, tanto na instância original de coberturas aos seis pares de blocos quanto na instância mista. Esta última remove a demanda formada pelos blocos 3 e 4; cada vértice dos blocos 1 e 2 tem prêmio 1, cada vértice dos blocos 3 e 4 tem prêmio 2, quota=6. Visitar apenas blocos 1 e 2 cobriria as cinco demandas mas daria prêmio≤4, de modo que três blocos opcionais continuam necessários.

Todos os ranks são calculados por eliminação modular no primo 1000003. Um menor não nulo nesse corpo é um menor inteiro não nulo, certificando limite inferior sobre rank racional; coincidir com os limites superiores torna as dimensões exatas. Os JSONs guardam as rotas geradoras que permitem reproduzir os menores/ranks. A validade da desigualdade é checada em todas as seleções factíveis antes de interromper a coleta de rank. Conclusão: **a desigualdade define faceta para estas duas instâncias específicas**. Não demonstra condição geral de faceta nem novidade bibliográfica; a certificação é forte evidência para orientar uma prova geral futura.
