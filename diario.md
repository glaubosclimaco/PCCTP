# Diário da pesquisa PCCTP

Este arquivo permite continuar a investigação em outra IA ou ferramenta sem
depender do histórico do chat. Criado por solicitação do usuário em
20/10/2026, às 12h42 no horário de Fortaleza (UTC−3). Os registros anteriores
abaixo são uma reconstrução a partir de arquivos e auditorias existentes,
identificada como tal; não são anotações contemporâneas inventadas.

## Leia isto para retomar

**Objetivo:** melhorar a resolução do Prize-Collecting Covering Tour Problem
(PCCTP) com evidência suficiente para um artigo científico. O usuário aceita
aplicação, adaptação e combinação de métodos conhecidos; não exige um método
inteiramente novo nem uso de ML. O artigo final deve ter aproximadamente
20–30 páginas, formulação, metodologia ilustrada, comparações justas, análises
pertinentes, referências verificadas e PDF compilado e inspecionado.

**Local e branch:** `/Users/glaubosclimaco/Projects/papers/PCCTP`,
`codex/pcctp-learning-submission`. Os caminhos relativos neste diário são
relativos à raiz do repositório. Confira `git status` antes de editar ou
selecionar arquivos para commit: há muitas alterações e artefatos históricos.
`include/Utils.h` contém uma alteração preexistente do usuário: preserve-a.

**Estado científico atual:** há ganhos exploratórios em TRAIN, mas ainda
não foi demonstrada superioridade sobre o MDM completo sob o protocolo de
70 iterações/10 sementes da literatura. Não há contribuição final ou artigo
pronto. ML e várias reparações inteiras tiveram controles negativos; não
apague esses resultados. Não iniciar CTP e não abrir as12 geometrias reservadas.

**Estado operacional confirmado em 20/10, após a retomada das 16h09:** um único
otimizador, `experiments/lkh/awake300_train.py`, retoma o lote de 81 execuções
de 300 s. Última auditoria: 69/81, 559 rotas e 3003 eventos, sete geometrias
completas. Um bloco de três execuções foi preservado por falha ambiental e
será repetido integralmente. Sessão atual: 13438. A fila durável
`published70_after_awake.py`, sessão 26890, aguarda esse lote e suas auditorias
para iniciar a calibração MDM70, seed80, nove casos TRAIN.
Esses IDs só valem na sessão da ferramenta atual: outra IA deve verificar os
processos reais e locks, e não presumir que um ID antigo possa ser reutilizado.

**Primeiro próximo passo:** acompanhar as 12 execuções restantes sem duplicar
o otimizador. O estudo isolado de cache já tem controlador e auditoria
verificados em dados sintéticos; não foi lançado nem enfileirado. Não o lançar
antes de terminar e auditar tanto o lote81 quanto as nove referências de
calibração. Ainda não existe resultado real de qualidade de solução do cache.

## Arquivos que contêm as decisões e evidências

- `AGENTS.md`: instruções portáveis para qualquer assistente, incluindo a
  obrigação de continuar este diário.
- `CLAUDE.md` e `memory/projects/pcctp.md`: decisões e restrições do usuário.
- `research_notes/EXECUTION_PLAN.md`: critérios científicos e objetivo final.
- `research_notes/CONTINUATION.md`: estado operacional detalhado, histórico
  de retomadas, sessões, artefatos e próximos passos. Use junto com este diário.
- `research_notes/DECISION300_AWAKE_REPORT.md`: diagnóstico do lote300s novo;
  parcial enquanto as81 execuções não terminarem.
- `research_notes/paired70_ten_seed_design.md`: comparação prospectiva com
  MDM70 completo, dez sementes pareadas. Código preparado, não lançado.
- `research_notes/cached30_conditional_design.md`: estudo de cache prospectivo,
  separado da comparação principal.
- `submission/main.tex` e `submission/build/main.pdf`: rascunho e PDF atual.
  Consulte também os arquivos de auditoria e proveniência em `submission/`.

## Como verificar a retomada

1. Leia o resumo acima e as entradas recentes de CONTINUATION. Confira branch,
   alterações locais, processos e arquivos de lock. Não remova um lock apenas
   porque seu arquivo existe: o bloqueio do kernel e o processo proprietário
   são o que importa.
2. Examine `experiments/lkh/awake300_train.queue.lock` e
   `experiments/lkh/published70_awake_calibration.queue.lock`. Não crie outra
   fila ou outro otimizador se os atuais estiverem ativos.
3. Para auditar resultados já gravados do lote300s, execute, sem concorrência
   com outro observador, `/usr/bin/python3 experiments/lkh/observe_awake300.py`.
   Esse comando não lança otimizador, mas atualiza relatórios. Só arquivos
   completos e blocos ambientalmente válidos sustentam a análise.
4. Não altere `screen.py`, `ruin_repair_absolute.py`, `awake300_train.py`, fontes
   de calibração ou binários congelados durante a coleta. Compare hashes
   contra o protocolo antes de afirmar continuidade do mesmo experimento.
5. Ao terminar um lote, confira custos e viabilidade, eventos dentro do prazo,
   suspensão do computador, cobertura de todas as sementes e geometrias,
   hashes e critérios definidos antes da coleta. Um arquivo completion sozinho
   não comprova validade científica.
6. Atualize este diário e CONTINUATION com o que realmente ocorreu. Não
   transforme previsões, testes sintéticos ou resultados parciais em conclusão.

## Histórico reconstruído:18–20/10/2026

### Escolha do problema e do tipo de contribuição

O usuário escolheu PCCTP, inicialmente com ML para selecionar vizinhanças por
cobertura/prêmio e reparação inteira. Depois esclareceu que melhorar soluções
por métodos conhecidos ou híbridos também é suficiente para investigar uma
publicação. A prioridade passou a ser desempenho demonstrado, não fabricar
novidade metodológica. A decisão editorial final depende dos experimentos.

O código legado MDM fornecido foi priorizado sobre uma reconstrução anterior.
Foram isoladas correções de índice não inicializado no decodificador e de
contagem de elite. Essas correções devem ser declaradas: o executável usado
não é uma reprodução byte a byte do programa publicado. O legado foi mantido.

### O que os experimentos completos já mostraram

- Híbrido com inicialização MDM e destruição/reparação/LKH,60s totais:
  54 execuções TRAIN, nove geometrias  × três sementes  × dois métodos;
  redução média por geometria de 1,886205%, oito ganhos e uma perda.
  Auditoria de 334 rotas. É seleção exploratória, não confirmação independente.
- Repetição em sete geometrias de validação já reutilizadas, cinco sementes:
  redução média de 0,907347%, dois ganhos, uma perda e quatro empates.
  O critério de investimento definido para essa etapa falhou. Não esconder
  isso nem chamar essa validação de teste independente.
- Reparo iniciado de estados comuns, comparado com apenas reordenar o mesmo
  suporte: redução média de 2,014272%, nove ganhos em nove geometrias.
  A geração dos estados está excluída do tempo: é evidência condicional de
  componente, não desempenho do algoritmo completo.
- Inserção conjunta por MIP, 324 comparações condicionais: variante de três
  lacunas próximas ficou em média 0,180897% pior que o reparo guloso, com apenas
  três geometrias favoráveis. O critério falhou. Não integrar por preferência
  pelo MIP quando a comparação observada não o apoia.
- Diversas políticas de ML, vizinhanças MIP, regras de redução e variantes
  tiveram resultados negativos. O inventário detalhado e os caminhos dos
  artefatos estão em CONTINUATION. Não repetir uma política apenas trocando
  sementes para procurar um resultado favorável.

Em quatro das nove instâncias TRAIN, os vértices obrigatórios já cobrem todas
as demandas. Logo, ganhos nelas não podem ser atribuídos à necessidade de
restaurar cobertura. Todas as nove ainda têm déficit de prêmio. Evidência em
`experiments/lkh/hard_train_instance_profile.json`.

### Comparabilidade com os artigos e suspensão do computador

A literatura heurística usa70 iterações e10 sementes;300s e3sementes são um
diagnóstico complementar, não o mesmo experimento. Tempos médios publicados
não são limites de execução e não devem ser convertidos por uma regra de GHz.
O plano pareado70 usa a duração real da referência completa de cada par como
orçamento de suas candidatas. Censura não equivale a completar70 iterações.

O primeiro lote300s foi afetado por suspensão do computador. Seus arquivos
foram preservados e marcados como inadequados à comparação temporal. Todas
as81 execuções foram reiniciadas em uma pasta separada, com monitor de relógios,
registro de suspensão e proteção de vigília. Não mesclar resultados bons dos
dois lotes. Uma rota viável não torna justo um experimento com tempo interrompido.

### Ferramentas, solvers e escrita

Os experimentos reais atuais usam MDM compilado e LKH 3.0.6, não Gurobi. As
investigações de MIP também usaram o backend disponível, explicitado nos
respectivos artefatos. Gurobi não tem licença ativa. Por pedido do usuário,
CPLEX foi instalado isoladamente em `.venv-cplex`; os testes indicaram licença
Community com limite de tamanho, não uma licença acadêmica irrestrita.

O rascunho LaTeX tem 18 páginas, formulação, metodologia, duas figuras TikZ,
gráfico e tabelas de estudos concluídos. A compilação e inspeção visual estão
registradas; abstract e conclusão continuam pendentes. Não preencher esses
trechos com superioridade ainda não demonstrada. As nove entradas TRAIN foram
comparadas por SHA256 com os arquivos públicos do depósito Mendeley V2; isso
não representa reprodução de todas as 107 instâncias do artigo ou de todo o
depósito. Referências reais e limitações estão nos ledgers de `submission/`.

## 20/10/2026,12h26–12h42 Fortaleza — estudo isolado de cache

**Ideia e motivo:** a serialização repetida da matriz de distâncias consome
parte da pequena janela de 0,2 s de cada chamada LKH. Guardar os textos inteiros
das distâncias pode reduzir essa preparação. O microbenchmark mostrou montagem
de texto mais rápida; isso não garante melhor solução. Construção e memória
do cache precisam ser cobradas, e uma busca real com controle é necessária.

**O que foi feito:** `cached_polish.py` implementa cache lazy por matriz e
preserva o contrato do polidor original. `cached30_phase.py` mantém o reparador
original e aceita somente melhorias até o prazo externo. O auditor escalar
independente confere rotas, custos, timestamps e contadores. Há marcadores de
tentativa para impedir repetição silenciosa depois de interrupção. Os
27 estados anteriores de TRAIN foram ligados aos hashes de resultados auditados
em `cached30_starts.py`/JSON. `cached30_pairs.py` controla ordem, duas fases,
ambiente e retomada dos pares.

**Verificação concluída:** 12 verificações do polidor,21 da fase e 8 da orquestração,
todos sintéticos, com certificados de fontes. Sem novos otimizadores ou leitura
de teste reservado. Commits relevantes: `a5e2276` e `743e4e6`.

**Avanço em andamento nesta retomada:** `cached30_batch_audit.py` foi escrito
para reconstruir independentemente os27 pares/54 fases, preservar o denominador,
calcular efeitos por geometria com frações exatas e conferir ambiente. Seu
fixture `verify_cached30_batch.py` passou 16 verificações sintéticos. O novo
`cached30_train.py` foi escrito para adquirir o mesmo lock do paired70,
exigir os dois estudos prioritários completos/auditados e verificar certificados.
**Esse controlador ainda não foi verificado e não foi lançado nem enfileirado.**
Se outra ferramenta retomar agora, não o execute como se estivesse pronto.

**Próximo passo concreto:** criar o fixture do controlador e conferir os certificados
após as últimas edições, testar bloqueio por pré-requisito incompleto e uso do
lock compartilhado, revisar documentação e salvar mudanças seletivamente.
O lançamento real do cache só poderá ocorrer após 81 execuções e nove calibrações70
completas/auditadas. O estudo de cache é condicional, 30 s por fase, não substitui
as duas candidatas principais do protocolo paired70.

## 20/10/2026,12h42 Fortaleza — solicitação de diário portável

O usuário pediu um minidiário durante todo o projeto, para poder transferir a
investigação a outra IA. Criei este arquivo com resumo de retomada, contexto,
evidências, limites e próximos passos; criei `AGENTS.md` e registrei a regra em
CLAUDE/memória. As próximas retomadas devem acrescentar entradas materiais,
incluindo fracassos e correções, e atualizar o resumo operacional. Este diário
não substitui os dados brutos nem as auditorias: aponta para eles e explica
por que cada investigação foi escolhida.

## 20/10/2026, 16h09–16h10 Fortaleza — falha ambiental revista e coleta retomada

Ao retomar, os IDs antigos de sessão já não existiam e uma consulta aos
processos confirmou ausência de MDM/LKH e das duas filas antigas. Havia 72
resultados finais, mas o bloco ts225/seed72 estava marcado ENVIRONMENT_FAILURE:
gap máximo do monitor 9,407546s (>5s) e diferença de incrementos dos relógios
3,608818s (>3s), sem evento de suspensão registrado. O coletor parou sozinho
nesse bloco; a ausência de suspensão não autoriza relaxar os outros limites.
Não é possível atribuir uma causa específica apenas a esses registros.

Conferi os 14 hashes congelados de fontes/dependências/binários e os hashes de
todos os resultados dos 24 blocos: íntegros; 23 blocos ambientalmente válidos.
O observador inicialmente recusou corretamente o bloco inválido. Preservei
os 16 artefatos desse bloco, incluindo todos os três métodos, em
`experiments/lkh/ruin_total300_awake_interrupted_blocks/20261020T190941Z/`.
`review.json` registra motivos, política de repetição e SHA256 de cada arquivo.
Não selecionei a repetição pelo custo encontrado e não aproveitei métodos
favoráveis do bloco. As 69 execuções válidas anteriores foram mantidas.

Reiniciei o coletor congelado, sem editar suas fontes, com
`/usr/bin/caffeinate -i /usr/bin/python3 -u experiments/lkh/awake300_train.py`.
Nova sessão 13438; faltam 12 execuções: repetir os três métodos de ts225/seed72
e executar os nove restantes de tsp225. A ordem e sementes são as do protocolo.
Implementei a fila durável `published70_after_awake.py`, sessão 26890, que usa
o lock já previsto, espera os 81 resultados, roda as auditorias e só então
inicia as nove referências seed80. Não foi lançado outro otimizador.

Nova auditoria dos 69 registros aceitos: 559 rotas e 3003 eventos; sete
geometrias completas. H10 tem ganho médio de −1,650450% (cinco ganhos, duas
perdas); H60, +1,253701% (sete ganhos, nenhuma perda). A regra H10 já não pode
ser satisfeita com duas perdas, mas ambos os braços serão concluídos, sem
interrupção seletiva ou ajuste. São resultados parciais TRAIN, não demonstração
de superioridade sobre MDM70. Próximo passo: fechar/auditar o lote e executar
a calibração completa; se outro bloco falhar, preservar e revisar antes de
qualquer nova tentativa.

### Consolidação do estudo de cache nesta retomada

O controlador `cached30_train.py` e seus testes já foram implementados no
trabalho iniciado anteriormente: 11 verificações sintéticas do controlador,
16 do auditor de lote, nove da orquestração de pares, 21 da fase e 12 do polidor.
Os certificados foram agora reconferidos contra os arquivos atuais: PASS.
O controlador usa o mesmo lock do paired70 e bloqueia enquanto os estudos
prioritários não estiverem completos e auditados. Nada foi lançado ou
enfileirado para cache. Durante a verificação houve falha de importação a
partir da raiz, pois o guard legado usa import local; a correção ficou no
adaptador novo, sem alterar a fonte congelada. Também foi reforçada a
checagem de todo o cronograma antes da primeira chamada: um par incompleto
no fim da ordem impede iniciar até os pares anteriores ainda não executados.
Esses testes validam contratos, não desempenho. As alterações e certificados
serão salvos seletivamente, mantendo `include/Utils.h` fora dos commits.

## Modelo para próximas entradas

### AAAA-MM-DD,HH:MM Fortaleza — assunto

O que fiz e por quê. Hipótese ou decisão que está sendo examinada.

Arquivos e comandos necessários para reproduzir ou retomar; processos ativos
e artefatos congelados, quando houver. Resultados com escopo, tamanho do lote,
auditoria e limitações. Falhas e correções sem apagar o resultado anterior.

Próximo passo concreto e condições que precisam ser satisfeitas antes dele.
