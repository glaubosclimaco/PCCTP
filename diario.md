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

Estado operacional confirmado em 22/10, 00h58: quatro referências aceitas,
semente 80. gil262 2001, rd400 10577, a280 2001 e gr229 96884. Auditoria
parcial PASS, 4 de 9, 60 verificações. A coleta atual é kroA200, semente 80,
MDM 14745, cerca de 3 minutos de CPU. O checkpoint tem 20 iterações, custo
provisório 23495, busca incompleta e 0 minerações. O log tem 17 marcações
e a última cabeçalho lida foi a iteração 16. O custo 23495 permanece
provisório. Não lançar outra calibração. O cache continua parado.

**Primeiro próximo passo:** acompanhar a calibração sem duplicá-la. Ela não
substitui a avaliação final de dez sementes.

**Rascunho conferido às17h17:**20páginas, incluindo resultados300s completos,
tabela de nove geometrias e figura dos dois braços. Ainda faltam abstract,
conclusão e evidências das referências/comparações70; não está pronto para
submissão. Compilação e inspeção estão em `submission/compilation_audit.json`.

## Arquivos que contêm as decisões e evidências

- `AGENTS.md`: instruções portáveis para qualquer assistente, incluindo a
  obrigação de continuar este diário.
- `CLAUDE.md` e `memory/projects/pcctp.md`: decisões e restrições do usuário.
- `research_notes/EXECUTION_PLAN.md`: critérios científicos e objetivo final.
- `research_notes/CONTINUATION.md`: estado operacional detalhado, histórico
  de retomadas, sessões, artefatos e próximos passos. Use junto com este diário.
- `research_notes/DECISION300_AWAKE_REPORT.md`: diagnóstico do lote de 300 s.
  Em 20/10, 17h10, está completo: 81/81, hybrid60 passou e hybrid10 não passou.
  Ainda não é a comparação com MDM de 70 iterações.
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

## 20/10/2026, 16h57 Fortaleza: tsp225 semente 71 fechou o bloco

O MDM da semente 71 gravou 2447, viável, orçamento nativo de 300 s. Com os
híbridos já em 2449, o bloco repete a semente 70: os dois híbridos ficaram
2 acima do MDM. O ambiente está PASS, folga máxima do monitor 2,158 s,
diferença de relógio 0,003 s, sem suspensão. A ordem desta semente foi
hybrid60, hybrid10, MDM.

Isso são duas sementes da mesma geometria, não a geometria completa e não
a média auditada de 72/81. A semente 72 está aberta no nativo de hybrid10
desde 16h55. Esse arquivo não entra no commit. Faltam os três braços dela.
Otimizador 68208 e fila 68316 seguem únicos. O fechamento descritivo continua
bloqueado até as 81 execuções.

## 20/10/2026, 17h00 Fortaleza — figura final preparada, ainda sem renderização

Conferi os78resultados e26blocos ambientais PASS observados no início desta
retomada; não havia falha ambiental. Depois, o coletor gravou H10/seed72/tsp225
com custo2449. Restam dois resultados do último bloco; não usar a semente
incompleta no agregado. O otimizador e a fila permanecem os mesmos. Não rodei
observador concorrente com o fechamento que a fila fará ao concluir81runs.

Li a habilidade `data:create-viz` e preparei
`experiments/lkh/plot_awake300.py`. A figura só aceitará o arquivo de
fechamento auditado completo: nove geometrias, dois braços e três sementes.
Ela terá18barras de médias e54pontos pareados, com o mesmo eixo para ambos
os métodos, perdas mantidas, cores acessíveis e valores rotulados. Os pontos
não serão apresentados como intervalos de confiança. O script recompõe os
efeitos a partir dos custos com frações e exige concordância com o fechamento.
Saídas previstas PNG/SVG e proveniência com dados/hashes; inspeção visual
continua obrigatória antes de incorporar ao artigo.

Checagem realizada: com o lote incompleto, a obtenção dos dados para o gráfico
foi recusada antes de criar imagem ou carregar renderização. Nenhuma figura
de resultados finais foi produzida ainda. Comandos futuros, em sequência,
depois das auditorias completas da fila: `close_awake300.py` e
`plot_awake300.py`, ambos em `experiments/lkh/` com `/usr/bin/python3`.
Próximo passo é analisar o lote completo, não preencher o artigo com previsão.

## 20/10/2026, 17h02 Fortaleza: hybrid10 da tsp225 semente 72 fechou em 2449

O JSON está viável. O pré-custo, o nativo de 10 s e o custo final são 2449,
o mesmo valor dos hybrid10 das sementes 70 e 71. O bloco da semente 72
continua aberto: o nativo de hybrid60 começou às 17h01 e não entra no commit.
Falta também o MDM desta semente. Sem esses dois, a geometria não fecha e a
média de 72/81 não muda.

Otimizador 68208 e fila 68316 seguem únicos. O gráfico final continua
recusado até o fechamento auditado das 81 execuções.

## 20/10/2026, 17h07 Fortaleza: hybrid60 da tsp225 semente 72 fechou em 2449

O JSON está viável. Pré-custo, nativo de 60 s e custo final são 2449, o
mesmo dos outros cinco híbridos já fechados nesta geometria. Nas sementes
70 e 71 o MDM ficou em 2447, 2 abaixo. O MDM da semente 72 ainda não fechou:
o log nativo e o checkpoint são de 17h07. Esse arquivo não entra no commit.
Sem ele, o bloco e a geometria continuam abertos, e a média de 72/81 não muda.

Otimizador 68208 e fila 68316 seguem únicos. Fechamento descritivo e figura
continuam bloqueados até as 81 execuções.

## 20/10/2026, 17h12 Fortaleza — lote de 81 auditado; calibração iniciada

Acompanhei o otimizador 68208 sem lançar outro. tsp225, sementes 70, 71 e 72,
gravou os três métodos. O coletor encerrou e escreveu
`experiments/lkh/ruin_total300_awake_train/completion.json` com 81 corridas.
A fila 68316, que já esperava esse arquivo, executou o observador. Não rodei
um segundo observador.

Auditoria de rotas: PASS, 81 corridas, 624 rotas, lote completo. Auditoria
ambiental: PASS, 27 blocos, elegível. Relatório
`research_notes/DECISION300_AWAKE_REPORT.md`. O critério prévio, em
`research_notes/awake300_repetition_protocol.md`, pede ganho médio de pelo
menos 0,2%, ao menos seis vitórias e no máximo uma perda em nove geometrias.
hybrid10 ficou em −1,279343%, com seis vitórias e três perdas: não passou.
hybrid60 ficou em 0,979441%, com oito vitórias e uma perda: passou. As duas
variantes permanecem no relatório. tsp225 piorou nos dois híbridos nas três
sementes, 2449 contra a média MDM 2447. Isso entra na média e não foi removido.

A fila em seguida iniciou `published70_calibration.py`. Às 17h12 o binário
MDM corre gil262, semente 80, 70 iterações, processo 65637, com caffeinate no
65629. O limite de segurança é 43200 s por caso. Nove geometrias de treino,
nenhuma reservada. Essa calibração mede a duração real da referência; não é
a comparação final de dez sementes e não autoriza escrever superioridade no
artigo.

Próximo passo: deixar a fila terminar as nove referências e a auditoria
`audit_published70_calibration.py`. Se um bloco ambiental falhar, o próprio
script preserva a tentativa e para. Não repetir em silêncio e não lançar o
cache antes dessa auditoria.

## 20/10/2026, 17h17 Fortaleza — fechamento81 e incorporação auditada ao rascunho

O coletor terminou com81execuções elegíveis e a fila realizou suas auditorias:
27blocos ambientais válidos,624rotas na auditoria de endpoint e3371registros
de rotas/eventos na reconstrução de horizontes. Os contadores se referem a
passagens distintas de auditoria, não a novas repetições experimentais.
H60 passou no critério exploratório: redução0,979441%, oito ganhos e uma
perda. H10 falhou:−1,279343%, seis ganhos e três perdas. tsp225 manteve os
seis resultados híbridos2449contra os três MDM2447. Todas as perdas e o
bloco ambiental inválido arquivado foram preservados.

Executei `close_awake300.py` com o lote completo: hashes e81rotas finais
conferidos, custos e critérios recompostos por frações, upstream preservado
em `experiments/lkh/awake300_completed_snapshots/`. Sensibilidade descritiva:
H60 varia0,840469–1,112088% ao omitir cada geometria; H10 varia−1,682921a
0,995688%. O extremo positivo de H10 exclui rd400 somente nesse diagnóstico,
sem excluir a instância da média principal ou alterar o critério. Conferi
independentemente18efeitos e18omissões por Decimal com45dígitos.

Gerei a figura por `plot_awake300.py`, com18barras e54marcadores. A primeira
inspeção mostrou rótulos próximos dos pontos; movi os valores para uma coluna
e ajustei o tamanho para leitura na página. Mantive o eixo comum, inclusive
a perda de rd400. Aritmética dos72valores conferida independentemente. Saídas
PNG/SVG e proveniência em `submission/figures/awake300_effects.*`.

Após esboçar a seção, acrescentei resultados e tabela300s, sensibilidade e
descrição da repetição ambiental ao LaTeX. Uma primeira aplicação de patch
falhou por contexto textual impreciso, sem alterar as fontes; reli e apliquei
o contexto correto. A redação preserva a validação60negativa, a reparação MIP
negativa e a distinção entre este screen TRAIN e MDM70. Não introduzi nova
referência ou alegação de novidade/BKS. Histórico das proveniências anteriores
foi copiado para `submission/provenance_history/`, e os vínculos atuais foram
reconferidos. O marcador da habilidade PDF foi executado uma vez antes da
autoria; duas compilações seriais pdfLaTeX terminaram sem avisos finais.

O PDF agora tem20páginas. Renderizei70dpi: páginas1–10 têm pixels idênticos
às revisadas antes; inspecionei individualmente11–20, incluindo as quatro
figuras e quatro tabelas. Não vi cortes ou sobreposição. Uma referência
isolada na última página será reequilibrada quando forem escritas as partes
pendentes, sem preencher páginas artificialmente. Abstract e conclusão
continuam pendentes; o número de páginas não significa artigo concluído.

A fila já iniciou a calibração seed80 por gil262. Conferência às17h17:
zero referências completas e sem marcador de falha ambiental; não iniciar
cache ou uma segunda calibração. Próximo passo: acompanhar checkpoints
apenas como observações parciais, preservar eventuais censuras e só usar
referências70completas/auditadas para dimensionar a comparação principal.

## 20/10/2026, 17h20 Fortaleza — acompanhamento passivo das nove referências

A calibração segue no único processo enfileirado. Ainda há zero registros
finais de referência e nenhum marcador de falha ambiental. O checkpoint de
gil262/seed80 observado nesta retomada registrava20iterações, custo2017 e
nenhuma chamada de mineração exportada. Isso é observação parcial, não um
resultado MDM70 nem evidência de ausência de mineração durante todo o resto
da execução. Não estimei o término por regra de três das iterações.

Criei `experiments/lkh/calibration_progress.py`, observador que não inicia
otimizador nem atualiza a auditoria oficial da fila. Ele mantém os nove slots
TRAIN do protocolo, distingue ainda não iniciado, checkpoint, conclusão pelo
runner, censura/falha e falha ambiental. Confere escalarmente a rota exibida
do checkpoint e só marca cobertura pela auditoria independente quando o hash
do registro final coincide. Um checkpoint com flags de conclusão, sem o
registro final e a auditoria, não basta para promover a referência.

O comando `/usr/bin/python3 experiments/lkh/calibration_progress.py` gerou
`experiments/lkh/published70_progress.json` e
`research_notes/published70_progress.md`, com instante UTC de observação.
Resultado atual:PROGRESS_ONLY, nove referências previstas, zero concluídas
e zero cobertas por auditoria. São arquivos passivos de situação que serão
atualizados; não substituirão os artefatos finais nem a decisão de calibração.
Uma projeção grosseira do custo das dez sementes/três métodos só será
calculada depois de nove referências completas e auditadas, com limitações
explícitas. Nenhum dado reservado, nova busca ou mudança de protocolo.

Próximo passo: observar a execução sem duplicá-la, relatar censuras se
ocorrerem e esperar as referências completas antes de dimensionar o lote
principal ou liberar o estudo de cache. O manuscrito de20páginas não mudou.

## 20/10/2026, 17h32 Fortaleza — guia de reprodução e manifesto de integridade

Atualizei o observador passivo: ainda não há referência completa. O checkpoint
de gil262/seed80 agora registra35iterações, custo2003 e zero chamadas de
mineração exportadas, após827,711s observados. Esses dados continuam
provisórios; não comparei esse custo com os endpoints nem estimei uma duração
para70iterações. Uma leitura paralela do relatório ainda mostrou a versão
anterior; reli depois do término do observador e usei o timestamp UTC novo.
Ao retomar, leia o relatório depois que seu comando terminar.

Criei `submission/REPRODUCIBILITY.md`, que registra denominadores, arquivos,
ambiente local, correções declaradas do baseline, comandos passivos e comandos
que lançam otimização. A separação evita que outra IA execute o coletor como
se fosse uma auditoria ou duplique a calibração. O guia inclui as limitações
de renderização entre máquinas e os pré-requisitos dos dois estudos futuros.
Versões locais conferidas:Python3.9.6, NumPy1.24.3, SciPy1.13.0,
Matplotlib3.8.4 e highspy1.13.1. O lote300s usa MDM/FPMax/LKH, não exige CPLEX
nem Gurobi para seus81resultados.

Gerei `submission/reproducibility_manifest.json` com183artefatos do repositório
e o hash do LKH externo. A lista deriva dos81resultados elegíveis,27blocos,
inputs das nove instâncias, fechamento/snapshots, fontes e binários congelados,
fontes/figura/PDF atuais e registros de revisão ambiental. Não inclui
checkpoints ativos da calibração, teste reservado ou a alteração do usuário
em`include/Utils.h`. Isso é integridade do escopo declarado, não reprodução
independente ou aprovação científica de todos os estudos anteriores.

Novo comando passivo: `/usr/bin/python3 experiments/lkh/verify_submission_manifest.py`.
Executado com183hashes de artefatos PASS e LKH externo PASS. Em outra máquina,
a falta do binário externo aparece separadamente da integridade dos arquivos.
Não há instalação, download ou otimizador nesse verificador. Mudanças futuras
no rascunho devem produzir uma nova versão revisada do manifesto, preservando
a anterior; não atualizar hashes apenas para esconder divergências.

Próximo passo: manter a calibração única em andamento e usar apenas referências
completas/auditadas para os próximos experimentos. O PDF de20páginas não mudou.

## Modelo para próximas entradas

### AAAA-MM-DD,HH:MM Fortaleza — assunto

O que fiz e por quê. Hipótese ou decisão que está sendo examinada.

Arquivos e comandos necessários para reproduzir ou retomar; processos ativos
e artefatos congelados, quando houver. Resultados com escopo, tamanho do lote,
auditoria e limitações. Falhas e correções sem apagar o resultado anterior.

Próximo passo concreto e condições que precisam ser satisfeitas antes dele.

### 20/10/2026, 23:13 Fortaleza — correção do manuscrito e retomada ambiental da calibração

Corrigi trechos que ainda descreviam o estudo300s como incompleto: introdução, protocolo, método e sensibilidade agora distinguem81ensaios elegíveis concluídos da comparação70iterações ainda pendente. Outline em research_notes/awake300_status_consistency_outline.md. A segunda passagem TeX foi iniciada antes da confirmação da primeira; corrigi executando duas passagens explicitamente seriais, logs consistency-serial-pass-1/2.txt. PDF20p sem warnings/overfull/underfull. Páginas1,2,3,8,10,16,17 inspecionadas diretamente; outras13 pixel-iguais ao render já revisado. Abstract/conclusão continuam pendentes. Proveniências anteriores arquivadas em submission/provenance_history; manifesto183artefatos+LKH novamente PASS. Nenhum resultado científico mudou.

A leitura atual do relógio foi23h09, bem posterior ao horário nominal do heartbeat. Não assumo trabalho contínuo nesse intervalo. O monitor mostrou ENVIRONMENT_FAILURE para gil262seed80: lacuna6,178795s excedeu limiar fixo5s, diferença de relógios0,326508s, nenhum Sleep. Nenhum processo otimizador estava ativo. Checkpoint70/custo2001 não constitui referência válida. Preservei todos os seis arquivos, com hashes, em experiments/lkh/published70_calibration_interrupted/20261021T021307Z/, review.json; conferi todas as fontes/binários congelados. Mantive protocolo e repeti a mesma semente/configuração, sem seleção por qualidade nem alteração dos limiares.

Fila única reiniciada: /usr/bin/python3 -u experiments/lkh/published70_after_awake.py, sessão71426, PID71865. Aguarda/conferirá awake81 e executa calibração, depois auditor oficial. Não iniciar cache/paired70 até9referências completas/auditadas. Próximo passo: verificar nova execução e ambiente, preservar qualquer nova falha; completar comparação científica antes de declarar artigo pronto. Utils.h e teste12 preservados.

### 20/10/2026 23:14 Fortaleza — confirmação da retomada e integridade

Conferi processos: fila71865, coletor71966 e único MDM72003 em gil262seed80, com caffeinate71967. Observador passivo registra somente checkpoint1, sem referência completa/auditada; não é resultado científico final. verify_submission_manifest.py reconfirmou183artefatos e LKH PASS após a correção textual. Nenhuma nova busca lançada. Consolidação seletiva das fontes/proveniências e arquivo da tentativa ambiental; Utils.h excluído. Próximo: deixar coleta única avançar e auditar referências completas; cache e paired70 continuam bloqueados.

### 20/10/2026, 23:24 Fortaleza — portabilidade do ensaio excluído

Calibração retomada avança: checkpoint10, custo provisório2017; zero referências finais/auditadas. Não iniciar outro otimizador. Conferi os seis hashes do arquivo da tentativa ambiental contra review.json e identidade do protocolo da repetição; preservation_audit.json PASS. A preservação local já estava completa, mas regras .gitignore excluíam raw/checkpoint/log do commit anterior. Incluí explicitamente esses três arquivos pequenos no arquivo histórico para outra IA também poder recuperar as evidências pelo Git. Não os reclassifiquei como elegíveis. Próximo: aguardar calibração completa e ambiente válido; não extrapolar prazo ou qualidade do checkpoint.

Nota de integridade: git diff --check aponta espaços finais emitidos pelo log nativo arquivado. Mantidos byte a byte para preservar os hashes; não são erro de fonte editada.

### 20/10/2026, 23h29 Fortaleza: retomada viva na iteração 17

Conferência passiva, sem novo otimizador. Os processos da retomada continuam
os mesmos da entrada de 23h13: fila 71865, coletor 71966, caffeinate 71967 e
MDM 72003. O log nativo de gil262, semente 80, está na iteração 17 de 70.
O relatório `published70_progress.md`, observado às 23h23, ainda dizia
checkpoint 10 e custo provisório 2017. Esse custo não é resultado de 70
iterações. `completion.json` da calibração não existe. Zero referências
fechadas. O filtro de treino das 81 execuções não mudou.

O aviso de 30 minutos não disparou entre 17h44 e 23h28. A falha ambiental
já registrada explica a interrupção; esta entrada só confirma que a repetição
está em andamento. Próximo passo: deixar esta única calibração terminar as
nove referências e auditá-las. Não extrapolar prazo a partir da iteração 17.

### 20/10/2026, 23:34 Fortaleza — verificação de congelamento durante coleta

Observador passivo atualizado: gil262seed80 checkpoint35/70, custo2003 provisório, zero referências completas/auditadas. Conferi novamente todos os hashes de fontes e binários no protocolo ativo: PASS. Manifesto do estudo concluído e manuscrito183artefatos+LKH também PASS. Nenhum parâmetro, fonte congelada ou otimizador adicional alterado/iniciado. Entradas de outro executor às23h29 preservadas. Não há conclusão científica nova; aguardar referências completas e auditoria antes de cache/paired70.

### 21/10/2026, 09:36 Fortaleza — calibração ainda sem referência final

Leitura passiva: checkpoint35/70, custo2003 e elapsed exportado1166,892s; zero referências finais/auditadas, nenhum environment_failure.json atual. ps confirma fila71865/coletor71966/MDM72003 único, MDM em estadoR e processo com elapsed cerca10h23. Não inferir tempo efetivo ou travamento a partir do checkpoint; monitor ambiental e limite43200s do protocolo continuam responsáveis pela elegibilidade/censura. Não reiniciei nem dupliquei a coleta. A primeira conferência ps expirou na revisão automática; a repetição autorizada somente de leitura funcionou, sem bloqueio restante. O timestamp real observado foi posterior ao heartbeat nominal; não reconstruo trabalho no intervalo. Próximo: revisar saída final/timeout e ambiente; cache/paired70 seguem bloqueados até9referências válidas.

### 21/10/2026, 09:37 Fortaleza — distinção entre checkpoint e log nativo

Conferi implementação congelada do coletor: timeout de segurança aplicado por subprocess.run; ClockGuard verifica lacunas e relógios e consulta pmset após retorno. Não alterei fonte ativa ou limiares. Checkpoint permanece35/70, enquanto log chegou49 com mineração; saídas têm cadências diferentes, não há novo endpoint. Nenhum arquivo final ou falha ambiental atual. Atualizei CONTINUATION, que ainda destacava iteração17, para evitar que outra IA reinicie uma coleta viva por interpretar relatório antigo como travamento. Próximo: observar conclusão/timeout e auditar ambiente; nenhuma busca adicional.

### 21/10/2026, 09h41 Fortaleza: log na iteração 62, ainda sem referência

O aviso de 30 minutos não disparou durante a madrugada. Conferência passiva.
Os processos continuam os da retomada de 23h13. O MDM 72003 está em execução,
98% de CPU, cerca de 33 minutos de tempo de processador e 10h28 de relógio
de parede. O log nativo tem 63 marcações de iteração; a última cabeçalho
completa lida foi a 62, com a linha seguinte ainda aberta. O checkpoint
registrado às 09h36 permanece 35/70, custo provisório 2003. Os dois canais
não são o resultado de 70 iterações. Não há `completion.json` nem falha
ambiental nova. Zero de nove referências. Não reiniciei a coleta.

O filtro de treino das 81 execuções não mudou. Próximo passo: esperar o
retorno desta instância e a auditoria ambiental do coletor. Cache e a
comparação de dez sementes continuam bloqueados.

### 21/10/2026, 09:47 Fortaleza — suspensão da tampa invalida repetição

Coletor retornou e rejeitou gil262seed80 por Clamshell Sleep às23h34 e vários Maintenance Sleep. Gap máximo3619,117s, diferença relógios3617,782s; duração processo37720,267s não é referência elegível. Custo2001 não entra em comparação. ps confirmou nenhum otimizador ativo. Preservei seis arquivos byte a byte em experiments/lkh/published70_calibration_interrupted/20261021T124719Z com review/hashes; fontes/binários congelados PASS. Não relaxei critérios. Retomar mesma configuração/semente com fila única; usuário informado de que computador precisa permanecer acordado. A coleta válida não está concluída. Próximo: acompanhar repetição e auditar9referências antes dos estudos posteriores.

Retomada efetivamente lançada: sessão95614, fila PID76634, published70_after_awake.py; não iniciar outra. Artefatos da tentativa excluída permanecem locais íntegros, pendentes de commit seletivo.

### 21/10/2026, 09:57 Fortaleza — arquivo ambiental conferido

Nova tentativa gil262seed80 em checkpoint10, zero endpoints finais/auditados; nenhuma busca adicional. Conferi os seis hashes da tentativa excluída20261021T124719Z e identidade do protocolo da repetição: preservation_audit.json PASS. Preparado commit seletivo com logs/raw/checkpoint históricos explicitamente incluídos, espaços nativos preservados; Utils.h e diretório ativo excluídos. Próximo: continuar coleta única sessão95614 e aceitar somente endpoints com ambiente válido.

### 21/10/2026, 10h07 Fortaleza — conferência da coleta e instruções atuais

Observador atualizado: checkpoint35 gil262seed80, custo2003 provisório; zero referências finais/auditadas. Log parcial está sendo gravado, com marcação34, o que não altera o checkpoint atômico nem representa regressão comprovada. Todos os hashes das dependências/fontes/binários do protocolo ativo conferidos PASS. Atualizei resumo de CONTINUATION com fila atual76634/sessão95614 para evitar uso dos PIDs históricos. Arquivo ambiental já consolidado no commit a3b7f5b. Nenhuma nova otimização, alteração experimental ou resultado científico. Próximo: aguardar retorno válido e auditoria, depois avançar estudos condicionados.

### 21/10/2026, 10h15 Fortaleza: a terceira tentativa continua em gil262

Conferência passiva. A execução que o resumo de 09h41 ainda descrevia foi
encerrada: a máquina dormiu, a lacuna chegou a 3619 s e o coletor rejeitou
o ensaio. Isso já está arquivado em `20261021T124719Z` e no commit `a3b7f5b`.
A tentativa agora ativa é outra, iniciada perto de 09h47. Às 10h15 seguem
vivos a fila 76634, o coletor 76702, o caffeinate 76703 e o MDM 76712, com
cerca de 28 minutos. O log nativo está na iteração 34. O checkpoint continua
35/70, custo provisório 2003, e não é referência de 70 iterações. Zero de
nove resultados fechados. Não lancei outro otimizador.

O filtro de treino das 81 execuções não mudou. Para esta referência valer,
o computador precisa permanecer acordado até o coletor gravar o ambiente.

### 21/10/2026, 10h18 Fortaleza — revisão do critério de conclusão70

Enquanto a coleta única continua, revi o auditor congelado sem executar auditoria final prematura: exige export final completed=true/iterations70, snapshot70, exit0/sem timeout, rotas/custos escalares, hashes e ambiente PASS (semSleep, gap≤5s, diferença≤3s). Só completion+9casos permite all_references_completed; calibração não é comparação de dez sementes. Nenhuma fonte alterada. Observador atual ainda35/70, zero referências finais. Não iniciei outra busca nem interpretei checkpoints como vantagem. Próximo: aceitar somente saída final auditada e avançar o protocolo condicionado.

### 21/10/2026, 10h28 Fortaleza — auditoria escalar do checkpoint

Log atual chegou62 com mineração, mas checkpoint atômico segue35/70. Auditei independentemente a rota e todos os snapshots desse checkpoint contra o inputTRAIN gil262 com raw_route_audit.load/check, sem chamar otimizador. Resultado em experiments/lkh/published70_checkpoint_scalar_audit.json: PASS_CHECKPOINT_ONLY; hash vincula a observação, não endpoint final. Essa verificação não aprova tempo/ambiente nem substitui o auditor final; zero referências completas. Nenhuma fonte congelada alterada ou busca duplicada. Próximo: aguardar final+ambiente e só então avaliar calibração.

### 21/10/2026, 10:38 Fortaleza — nova suspensão observada, tentativa excluída

Gil262seed80 retornou70, custo2001, mas ambiente rejeitado: Clamshell Sleep10h07 por155s e Maintenance Sleep78s; gap154,477s e diferença153,468s. Nenhum endpoint aceito. ps confirmou ausência de otimizadores. Preservei seis arquivos e review/hashes em experiments/lkh/published70_calibration_interrupted/20261021T133833Z; integridade e fontes/binários congelados PASS. Não afrouxar o guard nem apresentar tempo interrompido como comparação válida. Nova tentativa da mesma semente/protocolo, sem seleção por custo, depende de host acordado. Usuário informado. Cache/paired70 continuam bloqueados.

Nova fila efetivamente iniciada: PID79387, sessão75056, published70_after_awake.py. Não duplicar. Arquivo desta falha preservado localmente; commit seletivo pendente.

### 21/10/2026, 10h45 Fortaleza: quarta tentativa no começo

Conferência passiva. A rejeição das 10h38 está confirmada: 70 iterações e
custo 2001 não entram, porque a tampa dormiu 155 s e houve manutenção de
78 s. Às 10h45 a fila 79387, o coletor 79433, o caffeinate 79434 e o MDM
79449 estão vivos há cerca de 6 minutos, outra vez em gil262, semente 80.
O log nativo ainda tem 0 marcações de iteração. Não há referência fechada.
Não lancei outro processo. O custo 2001 da tentativa excluída não é baseline.

### 21/10/2026, 11h15 Fortaleza: primeira referência de 70 iterações aceita

A quarta tentativa de gil262, semente 80, gravou `COMPLETED70`. Custo 2001,
70 iterações, 14 rotas conferidas e 3 chamadas de mineração. O ambiente está
PASS: nenhum evento de suspensão, lacuna máxima 1,040193 s e diferença de
relógios 0,001453 s. O wrapper levou cerca de 2135 s. Os arquivos fechados
são `experiments/lkh/published70_calibration/g_gil262_VT51_T158_W53_75.seed80.json`
e o ambiente correspondente. O auditor independente das nove referências
ainda não foi chamado. Uma referência isolada não completa a calibração.

A mesma fila, processo 79387, iniciou rd400, semente 80, às 11h14. O MDM
85124 está em execução e o log dessa instância ainda está vazio. Não
dupliquei a coleta. O custo 2001 da tentativa rejeitada por suspensão não
foi promovido. O valor aceito é o desta corrida com ambiente PASS.

O filtro de treino das 81 execuções de 300 s não mudou e não se compara,
neste registro, com a referência de 70 iterações. Próximo passo: deixar
rd400 e as sete instâncias seguintes terminarem com o computador acordado.

### 21/10/2026, 11:17 Fortaleza — primeira referência70 com auditoria independente

Executei audit_published70_calibration.py sobre o único registro final atômico existente, sem interferir na coleta de rd400. PASS parcial:1run COMPLETED70,15verificações de rotas/custos (inclui reconferência do endpoint),70iterações,3chamadas mineração, custo2001, processo2129,847022s e wrapper2135,259169s. Fontes/binários/exports/input/ambiente vinculados por hash e critérios ambientais conferidos. all_references_completed=false, batch_complete=false: não é aprovação das nove nem comparação de dez sementes. Atualizei observador:1concluída/1auditada. Não reutilizei tentativas excluídas. rd400 já iniciado pela fila única; não duplicar. Próximo: auditar novos endpoints atômicos e só liberar estudos condicionados após9válidos.

### 21/10/2026, 11h28 Fortaleza — preservação da primeira referência elegível

rd400 iniciou e tem checkpoint1, custo10794 provisório; não comparei esse checkpoint. gil262 continua único endpoint70 concluído/auditado, SHA confere com auditoria parcial. Preparada preservação seletiva no Git do endpoint, ambiente e exports/log de gil262, protocolo e auditoria parcial; não incluir arquivos mutáveis de rd400, locks ou Utils.h. Isso permite recuperar a primeira referência válida sem promover o lote incompleto. Nenhum novo otimizador iniciado. Próximo: aguardar rd400 e demais referências; manter9casos como denominador.

### 21/10/2026, 11h45 Fortaleza: rd400 segue no checkpoint 10

Conferência passiva. gil262 permanece a única referência fechada e auditada,
1 de 9. O MDM de rd400, processo 85124, está em execução há cerca de 31
minutos. O checkpoint das 11h45 tem 10 iterações, custo provisório 10691,
`completed` falso e 0 chamadas de mineração. O log nativo continua com 0
bytes. O relatório das 11h37 ainda dizia iteração 5. Nenhum desses números
é resultado de 70 iterações. Não há falha ambiental nova nem `completion.json`.
Não reiniciei a coleta. Próximo passo: esperar o endpoint de rd400 com o
computador acordado.

### 21/10/2026, 11h58 Fortaleza — rd400 continua consumindo CPU

Checkpoint10 permanece provisório, sem endpoint novo. Para não interpretar ausência de export como travamento, conferi ps do MDM85124: estadoR,43min37s de parede e42min56,72s de CPU. Processo correto rd400seed80, coletor79433; não reiniciar nem duplicar. Essa observação operacional não valida ambiente final ou estima término. Permanece1/9auditada; próximoaguardar endpoint/timeout e auditoria ambiental.

### 21/10/2026, 12h15 Fortaleza: rd400 no checkpoint 20

Conferência passiva. O processo 85124 segue o mesmo, cerca de 61 minutos de
CPU e 98% de uso. O checkpoint das 12h05 tem 20 iterações, custo provisório
10691, `completed` falso e 0 minerações. O tempo interno desse arquivo é
3091 s. O log nativo, atualizado às 11h56, tem 18 marcações. A nota das
11h58 ainda descrevia o checkpoint 10. Nenhum canal fechou a referência.
gil262 continua 1 de 9. Não reiniciei a coleta e não estimei o horário de
término. Próximo passo: esperar o endpoint de rd400 com o computador acordado.

### 21/10/2026, 12h38 Fortaleza — rd400 ativo apesar do checkpoint estável

Checkpoint20 permanece parcial. Conferência ps: MDM85124 estadoR, parede1h24min09s, CPU82min57,93s, filho do coletor79433. CPU avançou desde11h58; não há evidência operacional para reiniciar. Sem endpoint novo ou environment_failure; ambiente final ainda não validado. Nenhum processo duplicado e nenhuma alteração de fonte. Próximo: aguardar endpoint/timeout, preservando critérios e denominador9.

### 21/10/2026, 12h45 Fortaleza: rd400 no checkpoint 35

Conferência passiva. O processo 85124 segue o mesmo, cerca de 91 minutos de
CPU. O checkpoint das 12h42 tem 35 iterações, custo provisório 10644,
`completed` falso e 1 mineração. O tempo interno desse arquivo é 5305 s.
O log das 12h34 tem 31 marcações. A nota das 12h38 ainda descrevia o
checkpoint 20 e o custo 10691. O custo provisório caiu 47 unidades entre
esses dois checkpoints. Isso não é a referência de 70 iterações. gil262
continua 1 de 9. Não reiniciei a coleta. Próximo passo: esperar o endpoint
de rd400 com o computador acordado.

### 21/10/2026, 13h15 Fortaleza: o log de rd400 passou o checkpoint

Conferência passiva. O MDM 85124 segue o mesmo processo, cerca de 119 minutos
de CPU e 98% de uso. O checkpoint não mudou desde 12h42: iteração 35, custo
provisório 10644, busca incompleta. O log foi atualizado às 13h01 e tem 45
marcações. A última cabeçalho lida foi a iteração 44, com a linha seguinte
ainda aberta. Candidatos visíveis nesse trecho, como 10727, não substituem
o custo do checkpoint. gil262 continua 1 de 9. Não há falha ambiental nem
`completion.json`. Não reiniciei a coleta. Próximo passo: esperar o endpoint
de rd400 com o computador acordado.

### 21/10/2026, 13:38 Fortaleza — rd400 excluído pelo guard de monitor

rd400 terminou70iterações/custo10577, mas gap8,329185s excedeu5s; clockdiff0,075641s, nenhumSleep. Não promover resultado após observar qualidade. ps confirmou coletor encerrado. Preservei protocolo e cinco arquivos da tentativa em experiments/lkh/published70_calibration_interrupted/20261021T163849Z com review/hashes; fontes/binários e arquivo íntegros. gil262 e auditoria parcial válidos mantidos intactos. Retomar somente rd400 mesma semente80/protocolo; sem afrouxar limiar ou reiniciar gil262. Usuário informado. Próximo: obter ambiente válido e nove referências auditadas.

Retomada iniciada: fila90912/sessão68012; wrapper pulará gil262 válido e repetirá rd400. Não duplicar. Arquivo desta tentativa local preservado, commit seletivo pendente.

### 21/10/2026, 13h45 Fortaleza: rd400 recomeçou na iteração 1

Conferência passiva. A rejeição das 13h38 está confirmada no arquivo
`20261021T163849Z`: checkpoint com 70 iterações, custo 10577 e `completed`
falso; `native_cost_after_return` 10577; saída 0; lacuna 8,329185 s;
diferença de relógios 0,075641 s; nenhum evento de suspensão. Não promovo
esse custo. gil262 permanece a única referência válida.

Às 13h45 a nova coleta está viva há cerca de 6 minutos: fila 90912, coletor
90968, caffeinate 90969 e MDM 90998. O checkpoint tem 1 iteração, custo
provisório 10794 e 0 minerações. O log nativo tem 0 bytes. Não lancei outro
processo. Próximo passo: deixar esta repetição de rd400 terminar com o
computador acordado e sem lacuna de monitor acima de 5 s.

### 21/10/2026, 14h31 Fortaleza: repetição de rd400 no checkpoint 5

Conferência passiva. Os processos continuam os das 13h39. O MDM 90998 está
em execução, cerca de 28 minutos de CPU, 52 minutos de parede e 98% de uso
agora. O checkpoint das 14h00 tem 5 iterações, custo provisório 10691,
`completed` falso, 0 minerações e 1286 s internos. O log segue com 0 bytes.
Não há falha ambiental nova. O custo 10577 da tentativa rejeitada continua
de fora. gil262 segue 1 de 9. Não reiniciei a coleta. Próximo passo: esperar
o endpoint desta repetição.

### 21/10/2026, 15h01 Fortaleza: repetição de rd400 no checkpoint 10

Conferência passiva. Os processos continuam os das 13h39. O MDM 90998 está
em execução, cerca de 54 minutos de CPU, 82 minutos de parede e 97% de uso
agora. O checkpoint das 14h34 tem 10 iterações, custo provisório 10691,
`completed` falso, 0 minerações e 3329 s internos. O log das 14h54 tem 18
marcações. O custo provisório não mudou em relação ao checkpoint 5. gil262
segue 1 de 9. Não há falha ambiental nova. Não reiniciei a coleta. Próximo
passo: esperar o endpoint desta repetição.

### 21/10/2026, 16:08 Fortaleza — rd400 ainda executa na repetição

Observador mantém checkpoint20/custo10691 provisório;1/9final auditada. ps confirma MDM90998, pai90968, estadoR,2h29min43s de parede e81min46,11s de CPU. CPU avançou em relação à nota15h01; diferença parede/CPU não prova suspensão nem aprova ambiente. Nenhum endpoint novo ou environment_failure atual. Sem reinício/duplicação. Próximo: aguardar export final e guard ambiental, preservar qualquer rejeição.

### 21/10/2026, 16h09 Fortaleza: checkpoint 20 ainda é o mais novo

Conferência passiva. Os processos continuam os das 13h39. O MDM 90998 está
em execução, cerca de 82 minutos de CPU e 2h30 de parede. O arquivo de
checkpoint mais recente é o das 15h04: 20 iterações, custo provisório 10691,
`completed` falso, 0 minerações e 5128 s internos. O log não cresceu desde
14h54 e permanece com 18 marcações. gil262 segue 1 de 9. Não há falha
ambiental nova. Não reiniciei a coleta. Próximo passo: esperar o endpoint
desta repetição.

### 21/10/2026, 16h39 Fortaleza: repetição de rd400 no checkpoint 35

Conferência passiva. Os processos continuam os das 13h39. O MDM 90998 está
em execução, cerca de 108 minutos de CPU, 3 horas de parede e 98% de uso
agora. O checkpoint das 16h25 tem 35 iterações, custo provisório 10644,
`completed` falso, 1 mineração e 10007 s internos. O log das 16h15 tem 31
marcações. O custo provisório caiu 47 unidades em relação ao checkpoint 20.
O checkpoint 35 da tentativa rejeitada às 13h38 também estava em 10644. Essa
coincidência não recupera o custo 10577 daquela tentativa nem fecha esta.
gil262 segue 1 de 9. Não há falha ambiental nova. Não reiniciei a coleta.

### 21/10/2026, 17h09 Fortaleza: o log de rd400 passou o checkpoint 35

Conferência passiva. Os processos continuam os das 13h39. O MDM 90998 está
em execução, cerca de 134 minutos de CPU, 3h30 de parede e 86% de uso agora.
O checkpoint não mudou desde 16h25: iteração 35, custo provisório 10644.
O log foi atualizado às 16h51 e tem 45 marcações. A última cabeçalho lida
foi a iteração 44, com a linha seguinte ainda aberta. Um candidato visível
nesse trecho, 10727, não substitui o custo do checkpoint. gil262 segue 1 de
9. Não há falha ambiental nova. Não reiniciei a coleta.

### 21/10/2026, 19:14 Fortaleza — confirmação do worker após intervalo de retomada

O relógio de observação indica22h13UTC, posterior ao heartbeat nominal21h20UTC; não infiro trabalho no intervalo. rd400 permanece checkpoint35/custo10644 provisório, sem endpoint ou environment_failure. ps: MDM90998 estadoR, pai90968, parede5h34min44s, CPU146min04,20s. Há avanço de CPU desde nota17h09, mas parede não equivale a duração elegível; guard ambiental decidirá ao retornar. Não reiniciei nem dupliquei;1/9auditada. Próximo: endpoint/timeout e auditoria, preservando qualquer rejeição.

### 21/10/2026, 19:25 Fortaleza — repetição de rd400 excluída por suspensão

Retornou70/custo10577 mas19Sleep (inclui tampa), gap3367,213s/clockdiff3366,155s. Não aceitar custo/tempo como referência; gil262 permanece1/9válida. ps confirmou ausência de otimizador. Preservei seis arquivos em experiments/lkh/published70_calibration_interrupted/20261021T222541Z com review/hashes, fontes/binários íntegros. Repetir só rd400 mesma semente80/protocolo, sem alterar limiares. Usuário informado da necessidade de host continuamente acordado. Próximo: obter endpoint com ambiente válido; não liberar cache/paired70.

Nova fila iniciada PID161/sessão5542. PIDs antigos históricos; não duplicar. Commit seletivo do arquivo desta falha pendente.

### 21/10/2026, 19h28 Fortaleza: nova repetição de rd400 na iteração 1

Conferência passiva. A rejeição das 19h25 está confirmada: checkpoint com
70 iterações, custo 10577, `completed` falso, 3 minerações e saída 0. O
tempo de parede nativo foi cerca de 20583 s, mas 19 suspensões produziram
lacuna de 3367,213 s. Não promovo o custo. gil262 permanece a única
referência válida.

Às 19h28 a fila 161, o coletor 197, o caffeinate 199 e o MDM 204 estão vivos
há cerca de 2 minutos. O checkpoint tem 1 iteração, custo provisório 10794
e 0 minerações. O log tem 0 bytes. Não lancei outro processo. O aviso de
30 minutos não disparou enquanto a máquina dormia. Próximo passo: deixar
esta repetição terminar com o computador acordado.

### 21/10/2026, 19:34 Fortaleza — relatório descritivo prospectivo entre sementes

Preparei experiments/lkh/paired70_dispersion.py separado de qualquer fonte congelada/runner/gate. Completa a exigência já prevista no desenho de reportar dispersão: média, mediana, mínimo/máximo, SD amostral n−1, vitórias/empates/perdas e pares brutos, por geometria com dez sementes completas. Geometrias incompletas mantêm contagem sem resumo de custos; sem p/CI, sem combinar sementes como geometrias independentes, sem mudar efeito primário ou critérios. Verificador reproduzível verify_paired70_dispersion.py PASS_SYNTHETIC_ONLY com contas manuais (variância55/6), empates/ganhos constantes, falhas/ausências/duplicatas. Execução real do novo relatório retorna NO_RESULTS: nenhum experimento paired70 lançado. Isso é preparação, não resultado científico novo. Coleta ativa rd400 permanece única, gil2621/9auditada. Próximo: aguardar calibração válida; executar relatório somente sobre auditoria do estudo paired70 futuro.

### 21/10/2026, 19h58 Fortaleza: repetição de rd400 no checkpoint 10

Conferência passiva. Os processos continuam os das 19h26. O MDM 204 está
em execução, cerca de 30 minutos de CPU, 33 minutos de parede e 98% de uso
agora. O checkpoint das 19h56 tem 10 iterações, custo provisório 10691,
`completed` falso, 0 minerações e 1860 s internos. O log segue com 0 bytes.
Não há falha ambiental nova. O custo 10577 das tentativas rejeitadas
continua de fora. gil262 segue 1 de 9. Não reiniciei a coleta. Próximo
passo: esperar o endpoint com o computador acordado.

### 21/10/2026, 20h28 Fortaleza: repetição de rd400 no checkpoint 20

Conferência passiva. Os processos continuam os das 19h26. O MDM 204 está
em execução, cerca de 60 minutos de CPU, 63 minutos de parede e 99% de uso
agora. O checkpoint das 20h20 tem 20 iterações, custo provisório 10691,
`completed` falso, 0 minerações e 3286 s internos. O log das 20h11 tem 18
marcações. O custo provisório não mudou em relação ao checkpoint 10.
gil262 segue 1 de 9. Não há falha ambiental nova. Não reiniciei a coleta.
Próximo passo: esperar o endpoint com o computador acordado.

### 21/10/2026, 20h58 Fortaleza: repetição de rd400 no checkpoint 35

Conferência passiva. Os processos continuam os das 19h26. O MDM 204 está
em execução, cerca de 89 minutos de CPU, 93 minutos de parede e 98% de uso
agora. O checkpoint das 20h57 tem 35 iterações, custo provisório 10644,
`completed` falso, 1 mineração e 5502 s internos. O log das 20h49 tem 31
marcações. O custo provisório caiu 47 unidades em relação ao checkpoint 20.
O checkpoint 35 das tentativas rejeitadas também estava em 10644. Essa
coincidência não recupera o custo 10577 nem fecha esta busca. gil262 segue
1 de 9. Não há falha ambiental nova. Não reiniciei a coleta.

### 21/10/2026, 21h28 Fortaleza: o log de rd400 passou o checkpoint 35

Conferência passiva. Os processos continuam os das 19h26. O MDM 204 está
em execução, cerca de 119 minutos de CPU, 122 minutos de parede e 99% de
uso agora. O checkpoint não mudou desde 20h57: iteração 35, custo provisório
10644. O log foi atualizado às 21h16 e tem 45 marcações. A última cabeçalho
lida foi a iteração 44, com a linha seguinte ainda aberta. Um candidato
visível nesse trecho, 10727, não substitui o custo do checkpoint. gil262
segue 1 de 9. Não há falha ambiental nova. Não reiniciei a coleta.

### 21/10/2026, 21:56 Fortaleza — rejeição rd400 sem suspensão e diagnóstico necessário

70iterações/custo10577 novamente, mas gap10,333994s>5s, clockdiff0,285471s, zeroSleep. Tempo de processo não elegível pelo protocolo; não recuperar tentativas anteriores ou afrouxar limiar após custos. ps sem otimizador. Preservei seis arquivos em experiments/lkh/published70_calibration_interrupted/20261022T005659Z e review/hashes; fontes/binários íntegros. Gil262 válido mantido. Nova repetição mesma semente/protocolo; falhas recorrentes semSleep merecem diagnóstico do monitor/carga antes de especificar eventual protocolo futuro, sem retroclassificar dados. Próximo: coletar referência válida e registrar todas as interrupções.

Nova fila6691/sessão98125 iniciada; não duplicar. Commit seletivo da falha pendente.

### 21/10/2026, 21h58 Fortaleza: nova repetição de rd400 no começo

Conferência passiva. A rejeição das 21h56 está confirmada: checkpoint com
70 iterações, custo 10577, `completed` falso e 3 minerações. Saída 0 e
tempo de parede nativo de cerca de 8439 s. Não houve suspensão. A lacuna
foi 10,333994 s e a diferença de relógios 0,285471 s. Não promovo o custo.
gil262 permanece a única referência válida.

Às 21h58 a fila 6691, o coletor 6736, o caffeinate 6737 e o MDM 6746 estão
vivos há cerca de 1 minuto, outra vez em rd400, semente 80. O log nativo
tem 0 bytes. Não lancei outro processo. Três tentativas desta instância já
devolveram 10577 e foram rejeitadas. Duas foram por suspensão ou por lacuna
sem suspensão. O limiar de 5 s permanece. Próximo passo: deixar esta
repetição terminar e guardar qualquer rejeição nova.

### 21/10/2026, 22:06 Fortaleza — inventário das interrupções da calibração

Criei calibration_interruption_ledger.py, observador sem otimizador ou alteração de fontes ativas. Reconfere hashes de todos os arquivos nos seis archives de tentativas excluídas, liga environment e native por instância e registra duração excluída/CPU/causa observada. Execução real sobre arquivos históricos:6tentativas,3comSleep e3semSleep registrado; soma de tempos de processo excluídos79919,829414s (não tempo válido, nem desempenho). Lacuna semSleep não identifica causa de escalonamento/carga; não há inferência causal ou mudança de limiar. Output calibration_interruption_ledger.json preserva razões em vez de esconder custo operacional das repetições. Nova rd400segue checkpoint1;gil2621/9auditada. Próximo: investigar futuras lacunas com medições prospectivas separadas, sem promover archives; aguardar endpoint da coleta única.

### 21/10/2026, 22h28 Fortaleza: repetição de rd400 no checkpoint 10

Conferência passiva. Os processos continuam os das 21h57. O MDM 6746 está
em execução, cerca de 31 minutos de CPU, 31 minutos de parede e 98% de uso
agora. O checkpoint das 22h24 tem 10 iterações, custo provisório 10691,
`completed` falso, 0 minerações e 1651 s internos. O log segue com 0 bytes.
A nota das 22h06 ainda descrevia o checkpoint 1. Não há falha ambiental
nova. O custo 10577 continua de fora. gil262 segue 1 de 9. Não reiniciei
a coleta. Próximo passo: esperar o endpoint desta repetição.

### 21/10/2026, 22h58 Fortaleza: repetição de rd400 no checkpoint 20

Conferência passiva. Os processos continuam os das 21h57. O MDM 6746 está
em execução, cerca de 60 minutos de CPU, 61 minutos de parede e 98% de uso
agora. O checkpoint das 22h48 tem 20 iterações, custo provisório 10691,
`completed` falso, 0 minerações e 3078 s internos. O log das 22h39 tem 18
marcações. O custo provisório não mudou em relação ao checkpoint 10.
gil262 segue 1 de 9. Não há falha ambiental nova. Não reiniciei a coleta.

### 21/10/2026, 23h28 Fortaleza: repetição de rd400 no checkpoint 35

Conferência passiva. Os processos continuam os das 21h57. O MDM 6746 está
em execução, cerca de 90 minutos de CPU, 91 minutos de parede e 99% de uso
agora. O checkpoint das 23h25 tem 35 iterações, custo provisório 10644,
`completed` falso, 1 mineração e 5280 s internos. O log das 23h17 tem 31
marcações. O custo provisório caiu 47 unidades em relação ao checkpoint 20.
O checkpoint 35 das tentativas rejeitadas também estava em 10644. Essa
coincidência não recupera o custo 10577 nem fecha esta busca. gil262 segue
1 de 9. Não há falha ambiental nova. Não reiniciei a coleta.

### 21/10/2026, 23h58 Fortaleza: o log de rd400 passou o checkpoint 35

Conferência passiva. Os processos continuam os das 21h57. O MDM 6746 está
em execução, cerca de 120 minutos de CPU, 122 minutos de parede e 98% de
uso agora. O checkpoint não mudou desde 23h25: iteração 35, custo provisório
10644. O log foi atualizado às 23h44 e tem 45 marcações. A última cabeçalho
lida foi a iteração 44, com a linha seguinte ainda aberta. gil262 segue 1
de 9. Não há falha ambiental nova. Não reiniciei a coleta.

### 22/10/2026, 00:17 Fortaleza — rd40070 válido auditado;2/9referências

Coletor gravou rd400seed80 COMPLETED70 com ambientePASS. Auditor independente executado sobre dois finais:PASS parcial,30verificações rotas/custos,2/9COMPLETED70,semcensura nos registros elegíveis,all_references_completed=false. rd400custo10577,3minerações,processo8229,090270s,wrapper8232,830226s;SHAfinal d001cf5cc2ff45bb15a5dec8afbadeccbd4926b85f6aa4c0b744fa54ddee70af. Custo coincide com tentativas excluídas mas apenas este endpoint/ambiente entra. Não comparar diretamente com300s ou garantir vantagem híbrida. Fila única já inicioua280seed80,checkpoint5provisório. Próximo:auditar demais sete referências e só liberar estudos condicionados após9válidas.

### 22/10/2026, 00h28 Fortaleza: a280 no checkpoint 35, custo provisório 2007

Conferência passiva. A fila 6691 e o coletor 6736 seguem os mesmos desde
21h57. rd400 já está fechado e aceito. O MDM atual é o 12697, instância
a280, semente 80, em execução, cerca de 14 minutos de CPU e 98% de uso.
O checkpoint das 00h27 tem 35 iterações, custo provisório 2007, `completed`
falso, 1 mineração e 802 s internos. O log das 00h26 tem 32 marcações. A
última cabeçalho lida foi a iteração 31, ainda em escrita. O custo 2007
permanece provisório. gil262 e rd400 seguem 2 de 9. Não lancei outra coleta.

### 22/10/2026, 00:48 Fortaleza — terceira referência70 auditada, a280

Auditor independente PASS parcial3runs/45verificações de rotas. a280seed80 COMPLETED70 custo2001,2minerações,processo1399,835346s,wrapper1403,585073s,ambiente válido;SHA2751f1be208874d4beb6bcad0e190ff813af8aeef4d856bcba6cefa2f89f5666.3/9auditadas;all_references_completed=false. gr229já em execução,checkpoint20provisório; não duplicar. Nenhuma comparação com híbrido nem alegação de melhoria. Próximo:demais seis referências válidas e auditoria total. Preservar somente arquivos finaisa280,sem gr229ativo.

### 22/10/2026, 00:57 Fortaleza — gr22970 auditado,4/9

Auditor independente PASS parcial4runs/60verificações rotas. gr229seed80 COMPLETED70,custo96884,3minerações,processo1020,129882s,wrapper1024,087607s,ambiente válido;SHA23d5f18d5adb1b46000b772aad9d2a0e00491b859bb9eb8d58edc00e78a6461a.4/9auditadas,all_references_completed=false. kroA200já em execução checkpoint20provisório. Preservar arquivos finaisgr229,sem fontes/conteúdo ativo de kroA200. Próximo:cinco referências restantes e auditoria completa; nenhuma alegação de ganho híbrido.

### 22/10/2026, 00h58 Fortaleza: kroA200 no checkpoint 20, custo provisório 23495

Conferência passiva. A fila 6691 e o coletor 6736 seguem os mesmos. a280 e
gr229 já estão fechados e aceitos, com a280 no custo 2001 e gr229 no custo
96884. O MDM atual é o 14745, instância kroA200, semente 80, em execução,
cerca de 3 minutos de CPU e 98% de uso. O checkpoint tem 20 iterações,
custo provisório 23495, `completed` falso, 0 minerações e 149 s internos.
O log tem 17 marcações. A última cabeçalho lida foi a iteração 16. O custo
23495 permanece provisório. As quatro referências seguem 4 de 9. Não lancei
outra coleta.

### 22/10/2026, 01:07 Fortaleza — kroA20070 auditado,5/9

Auditor independente PASS parcial5runs/75verificações rotas. kroA200seed80 COMPLETED70,custo23495,5minerações,processo322,000953s,wrapper325,687910s,ambiente válido;SHA1dc4ea775335025c737e73092819641d7a50fcf51bbb2c188d6e9c73f69c6d9e.5/9auditadas,all_references_completed=false. pr226já em execução,checkpoint20provisório. Nenhuma comparação híbrida. Preservar só arquivos finais kroA200,sem pr226ativo. Próximo:quatro referências restantes e auditoria completa.
