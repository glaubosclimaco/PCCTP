# Aprendizado de vizinhanças estruturadas para PCCTP

Direção escolhida: aprender quais decisões opcionais liberar conjuntamente,
considerando alternativas de cobertura, quota de prêmio e custo de circuito.
O MIP reotimiza a vizinhança; a política aprendida somente escolhe onde procurar.
Trata-se de uma hipótese de pesquisa, não de contribuição inédita já estabelecida.

## Pergunta e critérios científicos

**A estrutura conjunta de cobertura, prêmio e reparação do circuito permite
selecionar trocas melhores que escolhas aleatórias, custo/frequência simples e
aprendizado genérico, sob orçamento total comparável e em novas geometrias?**

A contribuição pretendida exige mais que acrescentar ML ao MIP:

1. Identificar complementaridade entre decisões de visita e suas alternativas de
   cobertura/prêmio, além da importância marginal de cada vértice ou aresta.
2. Propor uma representação e um operador que explorem essa complementaridade.
3. Isolar o efeito da estrutura por ablações e explicar em quais regimes ajuda.
4. Demonstrar generalização e benefício líquido, incluindo preparação dos dados,
   treinamento, inferência e resolução dos subproblemas.

Antecedentes obrigatórios:

- CL-LNS, ICML 2023: https://proceedings.mlr.press/v202/huang23g.html
- Supervised LNS, preprint 2025: https://arxiv.org/abs/2501.10778
- CMSA: https://doi.org/10.1016/j.cor.2015.10.014
- RL-CMSA: https://doi.org/10.1007/s10479-024-06243-7
- Os três artigos locais e o local branching da tese do autor.

Esses trabalhos impedem reivindicar novidade apenas por aprender o destroy
operator, usar feedback do solver ou combinar MDM-GRASP com MIP.

## O que está implementado

`pilot.py` cria uma solução inicial por vizinho mais próximo, remoções opcionais
viáveis e 2-opt. Não reproduz o GRASP histórico. Constrói candidatos aleatórios,
por sobreposição de cobertura e por prêmio; todos os seletores recebem o mesmo
pool e cada candidato libera exatamente k variáveis opcionais y e as arestas
incidentes nesses k vértices. Os demais x/y permanecem na solução inicial.
Assim, o número de variáveis x/y livres é igual dentro de cada instância.

O reparo usa uma formulação inteira com fluxo de uma mercadoria para assegurar
conectividade, implementada com scipy.milp/HiGHS. A formulação é correta para
soluções inteiras, mas sua relaxação não é o branch-and-cut dos artigos. Fluxos
continuam livres. Todos os candidatos são avaliados a partir do MESMO estado;
não existe ainda uma trajetória iterativa completa de LNS.

O decodificador verifica integridade, graus, um único ciclo, obrigatórios,
cobertura, quota, objetivo e respeito às variáveis fixadas. Um limite dual do
reparo nunca é apresentado como limite global do PCCTP. O tempo de construção,
solução e avaliação é registrado; o limite solicitado de solver não inclui todo
esse trabalho. scipy.milp não recebe MIP start: a solução corrente é preservada
externamente caso o reparo não devolva melhoria.

`evaluate_ranker.py` treina um rankeador linear ridge de diferenças entre pares
de candidatos do MESMO estado. Exclui da padronização e do treinamento todos os
estados da geometria avaliada. O modelo genérico usa visita, economia de remoção
e dispersão; o estrutural acrescenta prêmio e características de cobertura.
Alpha=1 foi fixado antes da avaliação e não ajustado nos resultados de teste.
Isto não equivale a teste confirmatório: o desenho é exploratório e somente
nove geometrias pequenas foram examinadas.

`verify_repair.py` compara o ótimo do reparo com enumeração independente de
circuitos admissíveis nas mesmas vizinhanças de uma instância pequena.

## Piloto de 18 de outubro de 2026

Seleção: primeira variante em ordem lexicográfica de cada geometria cuja
instância original tem N<=76. São nove geometrias diferentes. Configuração:
seed=42, k=5, 12 candidatos distintos e 0,5 segundo de solver por candidato.
O pool contém tipos em proporções variáveis devido à deduplicação; o controle
aleatório da avaliação é uniforme sobre esse pool compartilhado.

- 108 subproblemas; 108 ótimos locais certificados pelo solver.
- 27 candidatos melhoraram a solução inicial; sete dos nove estados tinham
  pelo menos uma melhoria entre os candidatos amostrados.
- A enumeração independente verificou dois reparos: 13 circuitos e ótimo 36;
  278 circuitos e ótimo 25. Ambos coincidiram com o MIP.

Redução relativa média do custo inicial ao escolher UM candidato por estado:

| Seletor | Redução média |
|---|---:|
| Uniforme no pool, valor esperado | 0,371% |
| Regra de economia de remoção | 0,383% |
| Rankeador genérico | 0,636% |
| Rankeador com cobertura/prêmio | 0,155% |

**O piloto NÃO demonstra vantagem das características estruturais.** Não há
intervalos de confiança ou afirmação de superioridade estatística, de tempo
total ou sobre os métodos publicados. Uma geometria influencia bastante a média.
O melhor candidato do pool é um oráculo AMOSTRADO, não a melhor vizinhança possível.
Melhorias sobre a inicialização simples não implicam melhorias sobre MDM-GRASP.

## Limitação estrutural que o próximo operador precisa enfrentar

Liberar apenas arestas incidentes aos opcionais escolhidos pode impedir o
atalho necessário para remover um vértice: a nova aresta entre predecessor e
sucessor pode ter ambos os extremos fixos fora do conjunto. O operador atual
permanece válido, mas pode ser demasiadamente restritivo. Custos e ganho de
cobertura, isoladamente, não representam essa dificuldade de reconexão.

A próxima hipótese será uma vizinhança com contexto de rota: liberar as decisões
opcionais e permitir a reconexão dos fragmentos afetados, mantendo um orçamento
explícito e comparável de variáveis livres. Devem-se estudar vizinhanças que
troquem visitados por não visitados e sinais de reparabilidade conjunta.
Esse operador será comparado ao incidente atual antes de atribuir ao aprendizado
qualquer ganho. Não aumentar a rede apenas para compensar um operador inadequado.

## Reprodução

Dependências: Python com numpy, scipy>=1.9 e networkx. O Python `python3` desta
máquina tinha scipy 1.13.0/numpy 1.24.3; o Python do runtime empacotado não tinha
scipy. Nenhum pacote foi instalado e o solver C++ não foi modificado.

```sh
python3 experiments/ml_lns/verify_repair.py
python3 experiments/ml_lns/pilot.py \
  data/g_att48_VT16_T16_W16_75 data/g_berlin52_VT16_T18_W18_50 \
  data/g_brazil58_VT11_T12_W35_75 data/g_eil51_VT17_T17_W17_25 \
  data/g_eil76_VT14_T16_W46_75 data/g_gr48_VT16_T16_W16_25 \
  data/g_hk48_VT16_T16_W16_50 data/g_pr76_VT14_T16_W46_25 \
  data/g_st70_VT14_T14_W42_25 \
  --seconds 0.5 --k 5 --candidates 12 \
  --output experiments/ml_lns/geometry_pilot.json
python3 experiments/ml_lns/evaluate_ranker.py \
  experiments/ml_lns/geometry_pilot.json \
  --output experiments/ml_lns/ranker_pilot.json
```

## Condições para avançar a uma avaliação de publicação

- Congelar protocolo após o desenvolvimento exploratório; reservar novas
  geometrias. A base atual tem 144 variantes de 48 geometrias: variantes irmãs
  não podem ser distribuídas entre treinamento e teste.
- Obter estados de soluções fortes e de diferentes fases da busca, incluindo
  MDM-GRASP; tratar trajetórias/geometrias como grupos, não cada candidato como
  observação estatística independente.
- Comparar LNS aleatória, heurísticas estruturadas sem ML, modelo genérico,
  seleção adaptativa de operadores, MDM-GRASP, MPO_S, local branching e métodos
  modernos de LNS aprendida adaptados. Não limitar o confronto a métodos antigos.
- Ablacionar cobertura, prêmio, contexto de rota, padrões conjuntos e modelo
  aprendido, com orçamento de ajuste comparável.
- Testar quotas, sobreposição/escassez de cobertores, tamanhos e geometrias novos;
  medir evolução do gap e custo total, além da qualidade final.
- Se o aprendizado estrutural não melhorar controles simples fora das famílias
  de desenvolvimento, rejeitar ou reformular a hipótese. Não prometer periódico
  nem originalidade com base neste piloto.

A sondagem anterior em `experiments/ml_branching/` é secundária. Sua relaxação
omite famílias de cortes dos artigos e não sustenta comparação de algoritmos
de branching. Serve apenas para diagnosticar oportunidades e reutilizar leitura
de instâncias, validada nas 144 instâncias locais.

## Continuação: contexto de rota e controles pareados

Implementado um operador que, além das arestas incidentes aos k opcionais,
libera até 2k arestas adicionais: atalhos entre os extremos de trechos removíveis,
arestas atuais cuja abertura permite inserções e um preenchimento por proximidade.
O número de arestas adicionais é truncado pelo número disponível.

Quatro expansões recebem exatamente o mesmo orçamento de x/y livres:

- `route_context`: atalhos e inserções informados pela rota;
- `random_extra`: expansão aleatória no conjunto de arestas não incidentes;
- `geometric_extra`: expansão por proximidade geométrica aos opcionais;
- `route_random`: expansão aleatória priorizando arestas atuais e atalhos,
  um controle mais forte que escolher entre todas as arestas.

O operador `incident` original libera menos variáveis; sua comparação serve como
diagnóstico, não como demonstração sob o mesmo tamanho de vizinhança.
Todos preservam a solução inicial como factível. A expansão contém a vizinhança
original; o programa verifica que um ótimo expandido não é pior que o original.

Desenho exploratório: as mesmas nove geometrias de desenvolvimento, sementes
42/43/44, k=5, 12 candidatos por pool e 0,5 segundo por reparo. São 27 pools
de candidatos de **nove estados iniciais distintos**, não 27 estados independentes.
A ordem de execução dos operadores foi embaralhada sem usar os objetivos.
O conjunto de controles fortes foi acrescentado após observar o primeiro
contraste com expansão uniforme; portanto esta análise é de desenvolvimento.

Foram resolvidos 1.620 reparos (324 por operador), todos com ótimo local provado.
Ganho relativo médio por vizinhança, agregando primeiro dentro de geometria:

| Operador | Redução média do custo inicial |
|---|---:|
| Incidente original, orçamento menor | 0,264% |
| Expansão uniforme | 0,272% |
| Expansão geométrica | 0,348% |
| Expansão aleatória orientada pela rota | 0,456% |
| Contexto de rota | 0,710% |

Nos 324 pares contra `route_random`, contexto ganhou 64, perdeu um e empatou
259. Contra `geometric_extra`, ganhou 78 e empatou 246. Contra expansão
uniforme, ganhou 91 e empatou 233. Isso mede a qualidade dessas vizinhanças
amostradas, não qualidade final nem velocidade de uma LNS completa. O fato de
liberar o mesmo número de variáveis não iguala o esforço efetivo de solução;
os tempos observados estão nos JSONs, e aqui todos os reparos terminaram ótimos.

Mudanças da seleção opcional ocorreram em 117 reparos de contexto, contra 84 de
`route_random`; outros ganhos podem vir somente da ordem do circuito. Ainda
precisamos isolar esses dois efeitos e avaliar incumbentes do MDM-GRASP.

O rankeador foi reavaliado sobre os pools de contexto, usando agora informações
de atalho, segmentos liberados e inserção. O controle genérico também recebeu
essas informações de rota, para não atribuir à cobertura/prêmio um ganho que
vem da reconexão. Todos os pools da geometria de teste foram excluídos de
treinamento e padronização. O treinamento mantém alpha=1 e continua exploratório.

| Seletor de uma vizinhança de contexto por pool | Redução média |
|---|---:|
| Uniforme no pool, valor esperado | 0,710% |
| Economia de remoção | 0,874% |
| Rankeador genérico com contexto de rota | 1,123% |
| Rankeador acrescido de cobertura e prêmio | 0,793% |

**O novo operador melhorou este diagnóstico, mas ainda não há vantagem demonstrada
do aprendizado específico de cobertura/prêmio sobre aprendizado genérico.** Os
324 pares não são 324 instâncias independentes. Não ajustar a arquitetura para
selecionar apenas resultados positivos; o próximo protocolo precisa reservar
novas geometrias e estados de soluções fortes. A melhora sobre uma inicialização
simples não estabelece melhora sobre os trabalhos publicados ou novidade do operador.

Reprodução dos controles fortes e avaliação do rankeador:

```sh
python3 experiments/ml_lns/verify_repair.py
python3 experiments/ml_lns/context_study.py --strong-controls \
  --output experiments/ml_lns/context_strong_controls.json \
  --ranker-input experiments/ml_lns/context_strong_ranker_input.json
python3 experiments/ml_lns/evaluate_ranker.py \
  experiments/ml_lns/context_strong_ranker_input.json \
  --output experiments/ml_lns/context_ranker_results.json
```

O verificador compara por enumeração três tamanhos de conjunto e cinco
operadores, e inclui um caso em que liberar a aresta de atalho permite uma
melhoria que a vizinhança incidente não consegue realizar.

O usuário pediu para guardar uma extensão futura ao CTP, depois de consolidados
os resultados do PCCTP. Essa decisão está em `CLAUDE.md` e
`memory/projects/pcctp.md`; nenhum trabalho de CTP foi iniciado nesta etapa.
