# PCCTP: publicações diretas e lacunas defensáveis

## O que os três trabalhos já cobrem?

### Takeaway
A combinação genérica de GRASP, mineração e programação inteira já foi motivada nos próprios trabalhos. Uma continuação precisa introduzir um mecanismo verificável e explicar quando ele funciona, em vez de apenas reunir os componentes.

### Cited Findings
- O Reactive GRASP é de 2021, RAIRO 55(3), 1441–1457, DOI 10.1051/ro/2021057. Usa construção GENIUS, RVND e adaptação reativa; compara formulações anteriores e heurísticas. Fonte primária: [artigo](https://www.rairo-ro.org/articles/ro/pdf/2021/04/ro200488.pdf). Texto local: `artigos/rairo/main.tex`, especialmente 752–773.
- New approaches saiu em IJISE 45(1), 101–134, **2023**, DOI 10.1504/IJISE.2023.133526 (não confundir com registro antecipado de 2021). Apresenta branch-and-cut, desigualdades, reduções e MPO_S. [Editora](https://www.inderscience.com/info/inarticle.php?artid=133526). O texto local `artigos/IJISE/sample (single col).tex`, 580–593, descreve remoção das arestas nunca vistas em ótimos locais; 1072 registra que parada precoce pode eliminar arestas promissoras; 1248 já propõe mineração + ILP como continuidade.
- Enhancing a GRASP heuristic ... through Data Mining Techniques saiu em IJLSM 53(1), 136–162, **2026**, DOI 10.1504/IJLSM.2026.150960, online 6 janeiro 2026, aceito em 23 abril 2023. [Editora](https://www.inderscience.com/info/inarticle.php?artid=150960). Texto local `artigos/mdm/main_review.tex`, 265 em diante, contém DM-GRASP e MDM-GRASP; conclusão 882–894 já sugere warm starts e seleção de variáveis para branching via padrões.

### Inferences
- Melhor reivindicação candidata: selecionar, reparar e expandir dinamicamente uma vizinhança exata usando estrutura de cobertura, prêmio e conectividade, com mineração como uma fonte de informação; isso ataca a perda de arestas promissoras documentada em MPO_S.
- Frequência individual de arestas é controle obrigatório: só há evidência de valor da mineração se coocorrência contribuir além de popularidade e custo geométrico.
- Limites de um MIP restrito não certificam o problema original. Manter um bound válido do modelo completo e distinguir explicitamente redução segura de restrição heurística.

### Gaps
- Não foi feita comparação linha a linha de todas as regras e formulações; novidade poliédrica requer essa análise e provas.
- A publicação em 2026 não significa que seus experimentos foram produzidos em 2026; a aceitação é de 2023.

## O que há de posterior ou próximo na literatura?

### Takeaway
Há trabalho recente diretamente relacionado que enfraquece propostas genéricas de adicionar veículos e recombinar rotas com solver. Existe também uma linha de modelos exatos de cobertura com coleta que deve entrar na revisão, sem tratar problemas diferentes como equivalentes.

### Cited Findings
- Assanvi e Roboredo, **preprint SSRN postado em 20 maio 2026**, Exact and Heuristic Approaches for the Multi-Vehicle Covering Tour Problem with Prize Collecting, DOI 10.2139/ssrn.6804412. O resumo anuncia branch-cut-and-price com formulação em rotas e HBR com ruin-and-recreate e set-partitioning sobre pool de rotas; cita os trabalhos PCCTP de 2021 e 2023. [Página dos autores no SSRN](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=6804412). O conteúdo indexado foi acessível, mas fetch direto retornou 403; classificação como preprint e resultados apenas conforme resumo, sem confirmação independente.
- Morandi, Leus e Yaman: Exact algorithms for budgeted prize-collecting covering subgraph problems, COR 144 (2022), 105798, DOI 10.1016/j.cor.2022.105798. Maximiza prêmio sob orçamento, permite cobertura capacitada, compara branch-and-cut e Benders para ciclo e árvore; portanto difere do PCCTP de custo mínimo com prêmio mínimo. [Editora](https://www.sciencedirect.com/science/article/abs/pii/S0305054822000867), [preprint dos autores](https://arxiv.org/abs/2101.03546).
- O benchmark público PCCTP de 2021 contém 144 instâncias TSPLIB e explicita prêmio mínimo não satisfeito somente pelos obrigatórios. [Dataset do autor](https://data.mendeley.com/datasets/gw76znr8bm/1). O depósito anterior contém resultados completos de formulações e heurísticas. [Resultados 2019](https://data.mendeley.com/datasets/gpnmh9kdr9/1).
- A tese do próprio autor já anuncia quatro heurísticas híbridas para PCCTP, combinando ILP e mineração; deve ser auditada antes de reivindicar uma combinação inédita. [Tese UFRJ](https://pesc.coppe.ufrj.br/uploadfile/publicacao/2964.pdf).

### Inferences
- Multi-veículo por si só não é uma proposta de novidade segura em outubro de 2026. Uma extensão deve definir necessidade operacional adicional e comparar com o preprint de 2026.
- Benders e branch-and-cut genéricos também não constituem novidade por si. A contribuição precisaria ser uma decomposição nova específica, cortes fortes, teoria e vantagem computacional.

### Gaps
- A busca focada por título e variantes PCCTP não encontrou sucessor single-tour direto que substitua claramente os três trabalhos; isso **não prova** inexistência nem autoriza chamar o código atual de estado da arte de 2026.
- O texto completo SSRN precisa ser obtido para comparar definições (quota, cobertura, obrigatórios, capacidade) antes de adaptar benchmarks.
- Fontes consultadas: editora, repositórios dos autores e datasets primários; agregadores usados apenas para descoberta.

## Qual pergunta vale testar agora?

### Takeaway
A direção de menor distância entre código e ciência é estudar o risco de confinamento em subgrafos heurísticos e um mecanismo reversível de expansão informado pela estrutura do PCCTP.

### Cited Findings
- Motivação empírica direta: MPO_S documenta perda de qualidade após eliminação precoce de arestas e planeja melhorar E0 com mineração. [Artigo publicado](https://www.inderscience.com/info/inarticle.php?artid=133526), texto local 1072 e 1248.
- Motivação de informação: MDM utiliza recorrência conjunta em soluções boas e reconhece dependência da qualidade/quantidade da amostra. [Artigo publicado](https://www.inderscience.com/info/inarticle.php?artid=150960), texto local conclusão.

### Inferences
- Hipótese falsificável: em igual orçamento de tempo, mineração conjunta + expansão guiada por feedback melhora a curva anytime e reduz confinamento, frente a frequência marginal, união de elites e MPO_S estático.
- Pipeline candidato: MDM gera incumbente/padrões; formar subgrafo mantendo incumbente; resolver vizinhança temporariamente; medir estagnação, cobertura crítica e gargalos; reabrir arestas; repetir. Relaxação completa fornece informação dual e bound global válido, quando computacionalmente viável.
- Ablations: mineração desligada; somente frequência; conjunto aleatório com mesmo tamanho; expansão aleatória; expansão determinística por custo reduzido; subgrafo estático; sem sinais de cobertura/prêmio. Separar ganho de melhor incumbente do ganho do mecanismo de busca.
- Dados: reproduzir 144 instâncias, criar gerador com tamanho, densidade de cobertura, correlação prêmio-custo e razão de obrigatórios controlados; separar tuning e teste por instância-base TSPLIB, não por variantes da mesma geometria.
- Medidas: gap primal global, tempo até alvo, curvas anytime, custo da mineração e do MIP, robustez por sementes; em pequenas instâncias com ótimo conhecido, medir retenção de arestas de pelo menos um ótimo e perda causada pela restrição. NÃO inferir ausência de ótimo só porque uma solução ótima específica perdeu arestas.
- Critério de abandono: se frequência + expansão simples empatar com padrões em instâncias novas, não forçar narrativa favorável à mineração; o resultado pode ser que adaptação e recuperação importam mais que padrões.

### Gaps
- Novidade externa do mecanismo candidato ainda depende de confrontar kernel search, local branching, RINS/RENS e mineração em matheurísticas; este memorando cobre PCCTP direto, não toda a literatura de matheurísticas.
