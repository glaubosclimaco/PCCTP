# Matheuristics e mecanismos de mineração para PCCTP

## O que já existe e limita a alegação de novidade?

### Takeaway
Combinar soluções heurísticas, restringir componentes e resolver MIPs já é uma família estabelecida. Acrescentar aprendizado ao ciclo também já foi publicado; MDM+MIP, sozinho, seria uma adaptação cuja força dependeria dos resultados e do mecanismo específico.

### Cited Findings
- CMSA (Blum, Pinacho, López-Ibáñez e Lozano, 2016) gera soluções probabilisticamente, reúne componentes numa subinstância, chama um solver e adapta a subinstância com o retorno. Isso é diretamente relevante para unir arestas de soluções do PCCTP e resolver o subgrafo. [Artigo original](https://doi.org/10.1016/j.cor.2015.10.014); [manuscrito dos autores](https://lopez-ibanez.eu/doc/BluPinLopLoz2015cor.pdf).
- RL-CMSA (Reixach e Blum, publicado em setembro de 2024) já treina online usando dados da busca. Valores de componentes sobem ou descem conforme sua presença na solução do solver; esses valores orientam novas construções. Portanto, realimentar mineração/construção pelo MIP também precisa ser distinguido desse antecedente. [Artigo original](https://link.springer.com/article/10.1007/s10479-024-06243-7).
- Local branching (Fischetti e Lodi, 2003) usa desigualdades de distância de Hamming para definir vizinhanças exploradas pelo MIP. Escolher adaptativamente o raio é uma extensão de um mecanismo clássico, não uma nova classe de hibridização. [Tutorial de um autor, Universidade de Bologna](https://site.unibo.it/operations-research/en/research/library-of-codes-and-instances-1/mic2003-lb.pdf/@@download/file/mic2003-lb.pdf); [artigo](https://doi.org/10.1007/s10107-003-0395-5).
- Kernel Search já organiza variáveis promissoras em um núcleo e expande com grupos de variáveis; a aplicação de Guastaroba e Speranza (2012) documenta essa linha. Restringir arestas promissoras com expansão, por si só, possui antecedente. [Artigo dos autores](https://doi.org/10.1016/j.ejor.2011.09.004).
- CL-LNS (Huang et al., 2023) aprende seleção de vizinhanças de ILP com exemplos positivos e negativos de uma heurística especialista e redes de atenção em grafos. Avalia qualidade ao longo do tempo, incluindo integral primal. Logo, usar rede neural ou aprendizado para selecionar variáveis livres também não basta como novidade. [Preprint original dos autores](https://arxiv.org/abs/2302.01578).

### Inferences
- A comparação conceitual mais próxima à proposta inicial do usuário é CMSA, além da própria MPO_S. Se o novo método só usar união de arestas com envelhecimento, estará reproduzindo sua estrutura.
- O diferencial científico deve demonstrar informação adicional dos padrões conjuntos além da frequência marginal de cada aresta, explorando a interação específica entre visita, prêmio, cobertura e conectividade.

### Gaps
- Esta consulta focada não demonstra ausência de trabalho idêntico especificamente para PCCTP; é necessário cruzar referências e citações da literatura do problema.

## Qual mecanismo seria uma hipótese científica concreta?

### Takeaway
A direção mais forte é uma vizinhança de substituição de cobertura orientada por padrões, que possa abandonar padrões frequentes quando eles bloqueiam mudanças coordenadas. É uma proposta de pesquisa, sem alegação de novidade já provada.

### Cited Findings
- Os antecedentes CMSA e RL-CMSA trabalham com componentes e utilidade observada na busca; CL-LNS ilustra o uso de sinais positivos/negativos para escolher vizinhanças. [CMSA](https://doi.org/10.1016/j.cor.2015.10.014), [RL-CMSA](https://link.springer.com/article/10.1007/s10479-024-06243-7), [CL-LNS](https://arxiv.org/abs/2302.01578).

### Inferences
- Proposta: representar não somente padrões de arestas, mas blocos estruturados com (i) caminho de visita, (ii) prêmio obtido, (iii) assinatura de clientes cobertos e (iv) alternativas de substituição de vértices opcionais. Mineração identifica coocorrência; o MIP testa se a relação ajuda a melhorar, em vez de apenas perpetuar incumbentes.
- Operador: escolher um bloco de vértices opcionais; liberar simultaneamente seus indicadores de visita, arestas incidentes e vértices alternativos que cobrem os mesmos clientes de W; manter as demais partes sob restrições locais. Essa liberação conjunta permite trocas que movimentos de aresta não alcançam facilmente.
- Manter padrões de soluções distintas que representam coberturas semelhantes; não fixar simultaneamente todos os padrões frequentes, pois podem ser incompatíveis em grau/conectividade ou descartar alternativa superior.
- Critério adaptativo: estimar melhoria por segundo para decisões de preservar versus destruir um padrão, condicionada a folga de prêmio, redundância de cobertura e diversidade do conjunto elite. Usar a utilidade de um padrão inteiro, comparada à soma das arestas, é a hipótese central, e não apenas colocar um bandit na busca.
- Hipótese falsificável: em instâncias com forte sobreposição de cobertura e muitas escolhas de visita, padrões com informação de cobertura proporcionam melhor integral primal que frequência de arestas, sob mesmo tempo total e mesmo número de variáveis livres. Se não vencer essas comparações, a hipótese não foi sustentada.
- Experimento diagnóstico anterior à implementação completa: nos casos pequenos com ótimo conhecido, medir cobertura de arestas ótimas, conservação de alternativas de vértices e gap ótimo do subproblema em função do tamanho do núcleo. Isso distingue economia computacional de mero viés em favor do incumbente.

### Gaps
- Não há garantia de que mineração estruturada seja superior. É preciso testar marginalidade, número de padrões, overhead e viés do conjunto elite.

## Como avaliar e preservar alegações corretas de exatidão?

### Takeaway
O trabalho deve separar melhoria primal de prova global e desenhar ablações que isolem padrões conjuntos, cobertura e feedback do MIP.

### Cited Findings
- Local branching descreve partições por uma restrição de distância e seu complemento; a estrutura completa de busca diferencia exploração local de certificação global. [Tutorial de Lodi](https://site.unibo.it/operations-research/en/research/library-of-codes-and-instances-1/mic2003-lb.pdf/@@download/file/mic2003-lb.pdf).
- CL-LNS explicita avaliação anytime usando primal gap e primal integral; são métricas pertinentes quando se disputa qualidade por orçamento de tempo. [Artigo](https://arxiv.org/abs/2302.01578).

### Inferences
- Para minimização, resolver com x_e=0 fora de E' restringe o conjunto viável. Uma solução viável é um limite superior válido para o PCCTP original, mas o dual bound da subinstância NÃO é em geral limite inferior global. Ótimo da subinstância não prova ótimo original. Infeasibilidade de E' tampouco prova inviabilidade original.
- Preservar um incumbente inteiro em E' garante viabilidade do subproblema somente se suas demais restrições forem mantidas corretamente. Não é suficiente preservar arestas sem visitas, cobertura, prêmio e conectividade.
- Certificação: executar separação e branch-and-cut no modelo integral completo, usando mineração apenas para starts, prioridades e heurísticas; alternativamente permitir exclusão de arestas apenas com prova válida (por exemplo, reduced-cost fixing com limite inferior global e incumbente, sob hipóteses adequadas).
- Baselines mínimos: MDM-GRASP original; MPO_S; solver completo com mesmo MIP start; CMSA adaptado; local branching; núcleo por frequência de arestas. Comparações opcionais mais caras: RL-CMSA e CL-LNS adaptados, sobretudo se a alegação for aprendizagem geral.
- Ablações: sem mineração; frequências marginais; padrões conjuntos; padrões conjuntos com cobertura; sem feedback; sem exploração; mesmo tamanho de núcleo e mesmo tempo total. Incluir mineração, construção e solver no tempo total, treinar/ajustar em famílias separadas, usar sementes pareadas e comparar por instância.
- Variar separadamente: quota de prêmio, proporção de opcionais/obrigatórios, raio e sobreposição de cobertura, geometria e tamanho. Buscar explicação de quando funciona, além de média de ganho.
- Se for proposta uma arquitetura geral, validar numa segunda família de cobertura com roteamento. Se o escopo for exclusivamente PCCTP, sustentar contribuição na estrutura e nos resultados desse problema, sem prometer generalidade não testada.

### Gaps
- Desempenho atual dos três algoritmos locais e validade de seus resultados precisam de reprodução. Não cabe estimar ganho percentual antes dos experimentos.
