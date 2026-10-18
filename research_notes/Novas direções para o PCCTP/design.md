# Agenda científica e desenho experimental para o PCCTP

## O que já existe e qual hipótese principal vale testar?

### Takeaway
A contribuição deve explicar quando a memória de soluções ajuda, como evitar seus vieses e por que as propriedades de cobertura e prêmio acrescentam informação à frequência. Apenas reunir MDM-GRASP, solver e local branching oferece novidade limitada.

### Cited Findings
- O manuscrito local minera padrões maximais de arestas com FP-max*, ordena-os pelo tamanho e usa seleção round-robin para iniciar GENIUS; MDM repete mineração quando o conjunto elite estabiliza. É a base existente, não uma proposta nova. [Manuscrito local, linhas 237–378](/Users/glaubosclimaco/Projects/papers/PCCTP/artigos/mdm/main_review.tex:237).
- A tese de setembro de 2020 já apresenta MPO_S, local branching, DM-GRASP e MDM-GRASP. MPO_S remove arestas nunca observadas e interrompe a heurística quando esse conjunto estabiliza; local branching permite explorar arestas externas ao conjunto observado, motivado pela ausência de uma aresta ótima. O texto inspecionado descreve estratégias separadas; o resumo não prova uma combinação operacional mineração+PLI. [Tese, §§3.5.3–3.5.5, páginas impressas 36–43](https://pesc.coppe.ufrj.br/uploadfile/publicacao/2964.pdf).
- Aprender tamanho de vizinhança e tempo de exploração de local branching já foi estudado, incluindo generalização por tamanho e entre instâncias; adaptar esses parâmetros por aprendizado não é novidade geral. [Liu, Fischetti e Lodi, versão de 2022](https://arxiv.org/abs/2112.02195).

### Inferences
- **H1, hipótese prioritária:** padrões que representam decisões conjuntas de visita, cobertura e prêmio obtêm maior melhoria por segundo que padrões escolhidos por suporte ou tamanho, especialmente quando há conflito entre cobertura e prêmio. Medir utilidade marginal sobre o estado atual: cobertura adicional sem duplicação, prêmio adicional até o déficit e custo de inserção/reparação. Prêmio e cobertura pertencem aos vértices, não diretamente às arestas: usar a união dos vértices de cada padrão e evitar dupla contagem.
- Primeiro fazer experiência controlada com o mesmo banco elite e mesmos padrões candidatos: seleção aleatória, round-robin por tamanho, frequência de arestas isoladas, frequência conjunta, utilidade proposta. Fixar orçamento e cardinalidade/tempo de construção; caso contrário padrões maiores parecem melhores porque simplesmente pulam mais inserções.
- Depois integrar o melhor seletor em uma matheurística que reabra decisões com baixa utilidade observada e preserve as de alta utilidade. Comparar explicitamente ao local branching da tese e a uma combinação direta MDM+MIP. O diferencial seria utilidade condicionada ao estado e mecanismo de revisão, não a presença do solver.
- Registrar melhoria de objetivo, tempo de reparação, viabilidade, sobrevivência de arestas após RVND e repetição de ótimos locais. Frequência e sobrevivência são descritores; o efeito causal do padrão exige intervenção controlada ou ablação.
- Não chamar suporte de probabilidade de uma aresta pertencer ao ótimo: ele é frequência dentro de amostra produzida pela própria busca.

### Gaps
- Não há evidência experimental ainda para H1. Novidade específica precisa de confronto com mineração de alta utilidade, memória adaptativa e matheurísticas de roteamento; esta nota não esgota tais literaturas.

## Como testar generalização e evitar conclusões restritas às instâncias antigas?

### Takeaway
Transformar a avaliação em uma pergunta sobre regimes de dificuldade: quando predominam conectividade, cobertura ou prêmio? A evidência mais forte combina base clássica, famílias novas controladas e teste fora da distribuição usada para ajustar parâmetros.

### Cited Findings
- O artigo MDM usa 144 problemas derivados da TSPLIB; 37 são usados para calibragem e 107 para avaliação. As variantes alteram conjuntos R/T/W e limiar de prêmio. [Manuscrito, linhas 398–415](/Users/glaubosclimaco/Projects/papers/PCCTP/artigos/mdm/main_review.tex:398).
- Recomendações de benchmarking abrangem objetivos claros, problemas bem definidos, algoritmos adequados, medidas de desempenho, análise, desenho experimental e reprodutibilidade. [Bartz-Beielstein et al., 2020](https://arxiv.org/abs/2007.03488).

### Inferences
- **H2:** ganho da memória depende da escassez de alternativas de cobertura e da tensão prêmio/custo; há regimes em que frequência induz convergência prematura. Criar fatores controlados: densidade e sobreposição das coberturas, número de cobertores por demanda, proporção de obrigatórios, déficit de prêmio após obrigatórios, correlação prêmio/custo e geometria uniforme versus agrupada.
- Dividir treino/validação/teste por geometria TSPLIB original, não aleatoriamente entre variantes da mesma geometria. Reservar tamanhos maiores e geometrias novas como teste final. Instâncias irmãs não são evidência independente de generalização geográfica.
- Para métodos aprendidos, transferir descritores normalizados e políticas, não identificadores de arestas. Um padrão de arestas identificado pelos números dos vértices não é diretamente transferível entre grafos.
- Baselines: GRASP reativo, DM-GRASP, MDM-GRASP, MPO_S, local branching da tese, solver completo com mesmo warm start, MDM+MIP direto, e uma heurística moderna simples de destruição/reparo com mesma infraestrutura. Ajustar todos sob orçamento de tuning semelhante.
- Experiência inicial de baixo custo: repetir mesmos bancos elite, variando apenas seletor e orçamento. Só avançar para aprendizado pesado se a informação de cobertura/prêmio superar frequência simples em teste separado.

### Gaps
- A separação 37/107 do artigo, sozinha, não demonstra separação por geometria; verificar lista das instâncias antes de criticar ou comparar.
- Não há quantidade universal correta de repetições. Um piloto deve estimar variabilidade e custo; 30 sementes por instância é ponto de partida proposto, não requisito estatístico nem garantia de poder.

## Que alternativa exata e que protocolo tornam a contribuição convincente?

### Takeaway
Uma contribuição exata pode ser mais forte que acrescentar outro componente heurístico, se produzir novas desigualdades comprovadas e limites melhores. Em qualquer caminho, distinguir ótimo do modelo restrito de ótimo do problema completo e contabilizar todo o custo computacional.

### Cited Findings
- O artigo de abordagens anteriores já inclui branch-and-cut, desigualdades e reduções; esses componentes devem ser referências de comparação. [Fonte local, artigo IJISE](/Users/glaubosclimaco/Projects/papers/PCCTP/artigos/IJISE/sample%20(single%20col).tex:86).
- A literatura metodológica exige que desenho, métricas e reprodutibilidade correspondam à pergunta científica. [Benchmarking in Optimization](https://arxiv.org/abs/2007.03488).

### Inferences
- **H3, alternativa:** combinar condições de cobertura, prêmio residual e conectividade gera cortes mais fortes do que tratá-las isoladamente. Antes de implementar um catálogo de cortes, procurar conjuntos de vértices cuja ausência inviabiliza cobertura ou prêmio e estudar como essa obrigação implica travessias de cortes. Provar validade, buscar exemplos de dominância ou independência das desigualdades existentes e avaliar separação. Não reivindicar família nova antes da revisão poliédrica.
- Para algoritmo exato, reportar relaxação na raiz, limite dual global, gap primal-dual global, nós, tempo de prova, instâncias resolvidas e custo da separação. Um limite inferior de um modelo com arestas eliminadas heuristicamente não certifica o problema original. Fixações por mineração exigem reabertura ou prova para preservar exatidão.
- Para heurísticas, orçamento principal igual de tempo total, incluindo construção do elite, mineração, reparo, montagem do modelo, solver e atualizações. Comparação com mesmas iterações pode permanecer como diagnóstico do mecanismo, mas não comprova vantagem em qualidade sob mesmo tempo.
- Mesma máquina, compilador, versão do solver e número de threads. Reexecutar baselines: não comparar tempos novos aos publicados em outro equipamento. Registrar uso de múltiplos núcleos e tempo de relógio; custo de treinamento offline deve ser declarado e avaliado por amortização quando aplicável.
- Métricas primárias pré-definidas: gap em tempos fixos e área sob curva de gap com normalização explícita; métricas complementares: tempo até vários alvos, taxa de sucesso e consumo de memória. Censurar execuções que não atingem o alvo no timeout, sem descartá-las. BKS é melhor conhecido; somente certificado recebe nome de ótimo.
- Usar validador independente sobre instância original: obrigatórios, cobertura, prêmio, graus, conectividade e custo. Soluções após reduções precisam ser reconstruídas e validadas no original.
- Resultados por instância e família, mediana/média/quantis por sementes e intervalos de confiança. Análise pareada por instância; reamostragem hierárquica por geometria/família e depois sementes. Não tratar todas as sementes de variantes irmãs como observações independentes. Wilcoxon/Friedman com ajuste de Holm pode complementar, mas tamanho do efeito e intervalo são mais informativos que contagem de p<0,05.
- Ablations mínimas: sem mineração; frequência simples; padrões conjuntos; utilidade cobertura; utilidade prêmio; ambas; seleção fixa versus adaptativa; elite de qualidade apenas versus elite com diversidade; reabertura aleatória versus informada, com mesma capacidade e orçamento.
- Publicar gerador, instâncias, sementes, configurações, logs com timestamps, rotas, bounds e hashes de versão. Definir antecipadamente comparações principais e limiar de relevância prática para impedir escolher apenas alvos favoráveis.

### Gaps
- A auditoria funcional do código e reprodução dos três artigos ainda não foram executadas nesta pesquisa. Inspeções anteriores do chat não substituem esses experimentos.
- H1–H3 são propostas falsificáveis, não descobertas estabelecidas, nem garantia de publicação.
