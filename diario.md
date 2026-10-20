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

Estado operacional confirmado em 20/10, 16h52: o otimizador 68208 e a fila
68316 seguem ativos. tsp225 tem cinco braços fechados. Semente 70, ambiente
PASS: MDM 2447, hybrid60 2449, hybrid10 2449. Semente 71: hybrid60 2449 e
hybrid10 2449, ainda sem MDM fechado. O MDM da semente 71 está em escrita
nativa desde 16h51. Faltam esse MDM e os três braços da semente 72. A
auditoria publicada continua 72/81. Não lance outro otimizador.

**Primeiro próximo passo:** deixar o otimizador 68208 e a fila 68316 terminarem.
Não lançar outro otimizador. O cache já tem controlador e auditoria sintéticos,
não foi lançado, e só pode rodar depois do lote de 81 e das nove calibrações
auditados. Ainda não existe resultado real de qualidade de solução do cache.

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

## 20/10/2026, 16h20 Fortaleza — diagnóstico passivo da perda de H10

Confirmei as sessões 13438 (único otimizador) e 26890 (fila aguardando). A
repetição de ts225/seed72 já gravou H60 e H10, ambos com custo108805; o bloco
de três ainda precisa terminar e passar pela auditoria ambiental. O resumo
aceito continua sendo69/81: não promover linhas isoladas a um bloco validado.

Investiguei por que o ganho médio de H10 ficou negativo quando rd400 entrou
no conjunto completo. Novo script passivo
`experiments/lkh/diagnose_awake300_initialization.py` lê somente os resultados
completos e ligados aos hashes das auditorias; conferiu escalarmente custos e
rotas iniciais/finais, reconciliação de contadores e efeitos por geometria.
Foram63resultados dos três braços, dos quais42híbridos, sete geometrias.
Nenhum otimizador ou dado reservado foi usado. Relatório legível em
`research_notes/awake300_initialization_diagnostics.md`, dados/proveniência em
`experiments/lkh/awake300_initialization_diagnostics.json`. As três auditorias
de entrada foram preservadas em `awake300_initialization_snapshots/`, pois os
observadores as atualizarão quando o lote avançar.

Em rd400, H10 perde19,479595% na média das três sementes. Nas sementes70e72,
não houve incumbente MDM tempestivo nos10s: a fase começou do fallback,
custos15499e14143, e terminou em13975em ambas. Os contadores dessas fases são,
respectivamente,1065chamadas/0exit0/1032timeouts/33other e
1187chamadas/0exit0/1186timeouts/1other. H60 melhora0,959129% em rd400.
Exit0 significa apenas retorno normal do processo, não uma rota melhor e
validada; contadores não fornecem decomposição detalhada do tempo por chamada.

Contraponto importante: pr226/seed70/H10 também não teve incumbente MDM no
prazo, mas recuperou o fallback88536para67286, com3189registros exit0. Além
disso, rd400/seed72/H60 melhora11486para10750mesmo sem exit0. Portanto, os
dados sugerem interação entre qualidade inicial e esforço de reparação/
reordenação; não isolam causalidade nem provam que cache resolva o problema.
Não eliminar rd400, trocar sementes, alterar limiares ou modificar candidatos.
Próximo passo permanece concluir os dois braços e a calibração MDM70.

## 20/10/2026, 16h26 Fortaleza — retomada sem duplicar o lote

O usuário pediu para continuar o trabalho no mesmo protocolo e manter este
diário. Confirmei os processos reais, sem usar os IDs de sessão anteriores.
O otimizador 68208 e a fila 68316 estão vivos e seguram os dois locks. Não
iniciei cópia, não alterei fonte congelada e não lancei o cache.

No disco há 72 resultados finais, oito geometrias com os três métodos e as
três sementes. O bloco repetido de ts225/semente 72 tem os três resultados e
o arquivo ambiental PASS. Isso ainda não substitui a auditoria de 69/81: o
observador não foi executado de novo, para não disputar a escrita com a fila,
que rodará `observe_awake300.py` quando existir `completion.json`. tsp225 está
em andamento no MDM da semente 70, binário
`build_timed_elite_bounded/original_mdm`, processo 76796.

Próximo passo: acompanhar as nove execuções de tsp225. Se o coletor parar sem
`completion.json`, preservar o bloco e revisar antes de qualquer repetição. A
calibração das nove referências de 70 iterações só começa pela fila já ativa,
depois das auditorias. O cache continua condicionado a esses dois lotes.

## 20/10/2026, 16h29 Fortaleza — bloco repetido aprovado e contabilidade revista

O observador terminou sem erro: 72/81 execuções, 582 rotas e 3063 eventos
auditados, 24 blocos ambientais válidos e oito geometrias completas. O bloco
repetido ts225/seed72 foi aprovado. Ganhos parciais: H10 −1,429044% (seis
geometrias com ganho, duas com perda); H60 +1,112088% (oito ganhos, zero perdas).
Falta tsp225 nas três sementes e três métodos. Não finalizar a comparação
antes de 81 resultados válidos nem chamar esse diagnóstico de comparação MDM70.

Revisei as fronteiras de tempo no estudo de cache e acrescentei uma tabela a
`research_notes/cached30_conditional_design.md`. Leitura, validação, preparação
e cache lazy consomem os 30s de busca; geração do estado comum está excluída.
Auditoria independente, proveniência, gravação e enumeração dos objetos do
cache ocorrem depois da busca e não fornecem novas soluções. `wall_seconds`
termina no retorno do núcleo: não representa tempo total de um processo
completo. Memória reportada é estimativa de objetos Python, não RSS/pico.
Essa revisão é de código, sem medições novas e sem alterar fontes ativas.

Motivo: evitar apresentar ganho de montagem de texto ou aumento de chamadas
como melhora de solução, ou chamar um estudo condicional de algoritmo
completo. Próximo passo: terminar o lote e deixar a fila iniciar a calibração
completa. Cache continua bloqueado e não enfileirado. Relatórios passivos
anteriores conservam seus snapshots e escopos originais. Uma primeira edição
deste resumo não encontrou o texto esperado porque a entrada de 16h26 já
havia atualizado o parágrafo; reli o arquivo e preservei aquela entrada.

## 20/10/2026, 16h38 Fortaleza: cópia enviada ao GitHub

O usuário pediu para colocar o projeto em outro PC. O remoto já era
`https://github.com/glaubosclimaco/PCCTP` e estava público. A branch
`codex/pcctp-learning-submission` não tinha upstream. O índice local, cerca
de mil arquivos de artigos, experimentos, notas e o rascunho em `submission/`,
foi commitado e enviado nessa branch.

Ficaram de fora do commit: `include/Utils.h`, alteração preexistente do
usuário, e quatro locks de fila (`awake300_train.queue.lock`,
`published70_awake_calibration.queue.lock`, `published70_calibration.queue.lock`,
`ruin_total300_train.queue.lock`). Um clone não herda esses locks. Os
processos desta máquina, se ainda estiverem ativos, continuam aqui. No outro
PC, confira processos e locks antes de lançar otimizador. Não trate o push
como fim do lote 300s nem como resultado científico novo.

No outro computador:

```
git clone https://github.com/glaubosclimaco/PCCTP.git
cd PCCTP
git checkout codex/pcctp-learning-submission
```

## 20/10/2026, 16h39 Fortaleza — revisão da calibração e correção preventiva

O único coletor continua na sessão13438 e a fila26890 espera o lote completo.
Já foram gravados MDM/seed70/tsp225 (2447) e H60/seed70/tsp225 (2449), mas ainda
não há um bloco completo novo para agregar. O snapshot publicado permanece
72/81. Não lancei outro otimizador ou ensaio de fidelidade concorrente.

Conferi o protocolo e o código da calibração: seed80, orçamento positivo
desativado, limite de segurança43200s tratado como censura, prioridade gil262
e rd400. A conclusão exige payload completo,70iterações, snapshot70, exit0,
ambiente válido e auditoria escalar. Resultados censurados permanecem no
denominador das nove referências e impedem chamar a calibração de completa.
Nenhuma dessas configurações foi alterada.

A revisão encontrou um erro no controlador de cache ainda não executado:
ele procurava LKH na lista de binários da calibração, que contém apenas MDM
e os auxiliares de mineração. Isso bloquearia o controlador por chave ausente
mesmo após referências válidas. Corrigi somente `cached30_train.py` para
conferir o LKH no protocolo awake300, que de fato o usa e congela. Acrescentei
dois testes sintéticos ao fixture: calibração sem LKH e recusa de binário
alterado. Agora são13verificações do controlador, todas PASS. Conferi também
o hash do LKH real contra awake300 e todos os certificados atuais: PASS,
sem executar busca. Não havia dados de cache para invalidar ou selecionar.

Próximo passo: deixar concluir os sete resultados de tsp225 ainda pendentes
no momento dessa conferência, auditar o lote e acompanhar as referências70
pela fila. Registrar tempos reais e censuras antes de estimar o custo do
estudo principal. O artigo permanece incompleto e o cache não foi enfileirado.

## 20/10/2026, 16h42 Fortaleza: tsp225 semente 70 fechada, monitoramento de 30 minutos

Conferência passiva, sem lançar observador nem outro otimizador. O lote
300s segue no processo 68208. A fila 68316 segue no lock de calibração.
Oito geometrias permanecem completas no disco. tsp225 semente 70 gravou os
três braços e o ambiente PASS. Custos fechados: MDM 2447, hybrid60 2449,
hybrid10 2449. Os dois híbridos perderam 2 unidades para o MDM nesta semente.
Não entra em média auditada e não muda a leitura de H10/H60 já registrada
em 72/81. Faltam seis execuções (sementes 71 e 72). O arquivo nativo
hybrid60 da semente 71 e o LKH 88149 indicam busca aberta, não resultado.

O usuário pediu acompanhamento e um texto de situação pelo menos a cada
30 minutos, diário nos fatos relevantes e commit quando houver avanço de
método, resultado ou decisão de experimento. A regra está em
`.cursor/rules/commitar-avanco-relevante.mdc` e em `AGENTS.md`. Um aviso
sem mudança de resultado não gera commit.

Próximo passo: deixar 68208 terminar tsp225 e a fila auditar o lote antes
da calibração. Cache continua bloqueado.

## 20/10/2026, 16h45 Fortaleza — tsp225 semente 70 fechou o bloco

O coletor não foi duplicado. Os três resultados de tsp225/semente 70 estão no
disco e o arquivo ambiental está PASS, com folga de monitor abaixo de 2 s e
sem evento de suspensão. A semente 71 já tem o artefato nativo de hybrid60; o
processo 68208 segue em CPU, então o intervalo entre o nativo e o JSON final
é processamento, não parada. Ainda faltam seis execuções. Nenhum custo novo
entra na comparação antes do observador. A fila 68316 continua aguardando.

## 20/10/2026, 16h46 Fortaleza: hybrid60 da tsp225 semente 71 fechou em 2449

O braço hybrid60 da semente 71 gravou JSON viável, custo 2449, igual ao
pré-custo e ao nativo, orçamento nativo de 60 s, relógio do braço cerca de
300 s. Ainda não há MDM desta semente, então 2449 não se compara com o 2447
da semente 70. O bloco da semente 70 permanece o único completo em tsp225:
os dois híbridos ficaram 2 acima do MDM. Isso é uma semente, não uma
conclusão de método.

O nativo de hybrid10 da semente 71 começou às 16h45 e continua aberto. Não
foi commitado. Faltam esse braço, o MDM da semente 71 e os três braços da
semente 72. Otimizador 68208 e fila 68316 seguem únicos. Cache bloqueado.

## 20/10/2026, 16h50 Fortaleza — preparação do fechamento, sem análise parcial nova

Conferência passiva:77arquivos finais,25blocos ambientais PASS, nenhum bloco
com falha e quatro resultados ainda ausentes. A fila continua esperando a
conclusão. Não atualizei os agregados publicados de72/81 apenas por haver
arquivos de um bloco ainda aberto; a fila fará o observador completo antes
da calibração. Nenhum arquivo nativo em escrita foi selecionado para commit.

Preparei `experiments/lkh/close_awake300.py` para executar depois do término e
das auditorias finais. Ele exige81chaves esperadas, nove geometrias,27blocos
ambientais validados e consistência das auditorias/curvas. Reconfere hashes e
custos/viabilidade das81rotas finais, recompõe os dois ganhos por geometria
com frações exatas e verifica concordância com os critérios já definidos.
Preserva uma cópia das auditorias/protocolo e calcula sensibilidade descritiva
ao retirar uma geometria por vez, mantendo todos os dados e sem criar outro
critério de aprovação ou teste de significância. Isso permitirá mostrar se
o efeito médio depende muito de uma geometria, sem eliminá-la da conclusão.

Verificação realizada agora: o script recusou corretamente o lote incompleto,
antes de criar análise de fechamento, snapshot ou otimizador. A execução
completa ainda não ocorreu e precisa ser conferida após a fila publicar as
auditorias81. Comando futuro: `/usr/bin/python3 experiments/lkh/close_awake300.py`.
Não o incluir no processo ativo nem substituir o protocolo300s pelo de70.

Próximo passo: terminar os quatro resultados, acompanhar a calibração já
enfileirada e fechar a análise descritiva quando os pré-requisitos forem reais.
O artigo e as inferências permanecem condicionados às evidências completas.

## 20/10/2026, 16h52 Fortaleza: hybrid10 da tsp225 semente 71 fechou em 2449

O JSON fechou viável, custo 2449. O nativo de 10 s devolveu 2483 e o
pós-processamento chegou a 2449, o mesmo custo do hybrid60 desta semente.
Ainda não há MDM da semente 71, então os dois híbridos não têm par nesta
semente. Na semente 70, os dois híbridos continuam 2 acima do MDM 2447.
Não atualizei a média de 72/81.

O MDM da semente 71 está aberto: log nativo às 16h51 e checkpoint às 16h52.
Não foi commitado. Faltam esse braço e os três da semente 72. Otimizador
68208 e fila 68316 seguem únicos. O fechamento descritivo continua recusado
enquanto o lote não tiver 81 resultados.

## Modelo para próximas entradas

### AAAA-MM-DD,HH:MM Fortaleza — assunto

O que fiz e por quê. Hipótese ou decisão que está sendo examinada.

Arquivos e comandos necessários para reproduzir ou retomar; processos ativos
e artefatos congelados, quando houver. Resultados com escopo, tamanho do lote,
auditoria e limitações. Falhas e correções sem apagar o resultado anterior.

Próximo passo concreto e condições que precisam ser satisfeitas antes dele.
