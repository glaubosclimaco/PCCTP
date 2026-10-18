# Estado operacional para retomada autônoma

## Retomada às 7h11 — prazo comum, falha de relógio detectada e repetição

- Implementado --timed embuild.py, build_timed/ separado, originalsemwrapper
 detemporização; sóadaptive medeoperadores. Hooks acompanham construção e
 melhoriasglobais dentro da buscalocal, guardamanytime e deadline_cost somente
 se descoberta<=prazo; melhorias tardias separadas. Iterações livres sobB,
 estabilidade de mineração mantémnumIt70. Parada cooperativa entreoperadores;
 operaçõesatômicas podemultrapassar e issoéregistrado. Nenhumsoluçãoantesdoprazo
 viraNone,não usa solução tardia como fallback comparável.
- `verify_deadline.py` passou exclusãodemelhoratardia,débito startupenovo teste
 deépocanativaentreleiturasCLOCK_MONOTONICdoPython. Semorçamento,modooriginal
 reproduziu2650em70iteraçõesatt48. Parentstartdebited,parse/preparaçãocobradas.
- PRIMEIROLOTEccb595d69e76 terminou162runs(9DEV,3seeds,6políticas,1s),mas é
 INVALIDADO PORFALHA DE RELÓGIO,registrado eminvalid_timing_protocols.json.
 Auditordetectougr48/seed44/adaptive:nativo3.0515s vsparent0.9121s.
 ProbeconfirmouGCCsteady_clockretornaépocadewalltime~1.792bilhões, enquanto
 CLOCK_MONOTONICexplícito~421mil. Nãointerpretar esse loteou seuresumo como
 evidência de política. Todo lote/fontes/evidêncianativa preservados. Não
 selecionar sóoutlier:sumarizador recusa oPROTOCOLOINTEIRO,auditoria oexclui.
- FIX:research_hook.cpp ewrapperadaptiveTIMED usamclock_gettime(CLOCK_MONOTONIC)
 explícito,paitambémusaCLOCK_MONOTONICepassaPCCTP_STARTED_MONOTONIC. CLOCKSYS
 ePythonperf_counter não entram maisnocálculodeprazos/elapsed. `clock_source`
 salvo. .process_wallseconds inclui polling; .native_elapsedseconds já
 explicitamente salvo. Auditoria preservanative_records/ e calcula ambos.
- RepetiçãoINTEGRALagoraATIVA sessão36046:
 `python3 experiments/original_mdm/timed_experiment.py --budget 1 --limit 9
 --seeds 42 43 44`,log `/tmp/pcctp_timed_monotonic_development.log`.
 Protocolo d0df6114da43. Não duplicar. ÚltimaauditoriaPASS54arquivos novos,
 nenhumproblema derelógio. Antigoprocesso29568 encerrou. Builds39944/check47277/
64538encerrados. Programa novo deveacabar em~3min desdeinício.
- Novos scripts audit_timed.py/summarize_timed.py;sumário descritivo agrega por
 geometria,preserva pareamentos/missing/semviáveis.Nãohátestehypothesis.
 Rodarapósfim:python3 experiments/original_mdm/audit_timed.py e
 python3 experiments/original_mdm/summarize_timed.py d0df6114da43.
- IMPLICAÇÃOCIENTÍFICA:perfisantigosdeoperadores epath_dp.cpp usavamomesmo
 chrono. Custos/viabilidade/enumeração continuamválidos,mas médiastemporais e
 prazos precisamrevisão/recoleta. Pesosmodel52b8ainda congeladosANTIGOSpara
 diagnósticodarepetição; conclusãofinal exigeTRAIN comCLOCK_MONOTONIC explícito.
 Não misturar perfis/reutilizarmanifestantigosemverificarclock. Coletorhoje
 reusecompara build_manifest;novo manifesto de profiling deveincluir código
 instrumentado/hookseclock paraevitarreusoerrado. Usardiretóriodedatasetnovo.

Próximo: finalizarrepetição/auditar,semclaimantesdisso. Depois corrigirwrapper
 profiling(tambémdeveusarclockexplícito;hojesóTIMEDfoifixado), acrescentarhashes
 instrumentados/clockaomanifestoe recoletarTRAINemdiretórioseparado para
 retrainsemperderlegado. Fonteversõesloteinválidoem timed_source_versions/.
 Preservarresultadosnegativos,9DEVsomente; testesreservadosintactos,CTPfora.
 Nenhumacontribuiçãosuficiente/artigoprontoainda, nãoencerrarautomação.


## Retomada às 7h01 — política nativa e primeiro piloto completo

- `operator_state.cpp` + `verify_operator_state.py`:descritoresnativos vsPython
460estadosde23perfisPASS,erro6.66e-16,medianaobservada5.5us/estado viactypes,
 exclui preparaçãoúnicadografo. Outputoperator_state_verification.json.
- `scheduler.hpp` adaptaVertice/C/rota paraíndices/cobertura/prêmios ecarrega
 modelos52b8e9ae6089. Adaptadorvalidadoem60estados3DEV,erro2.22e-16. Produz
 features.jsonl apenasquandoPCCTP_SCHEDULE_VERIFYativado; pilotosemI/O extra.
- `build.py --initialize-decoder-index --scheduled` criaBUILDSEPARADO
 build_scheduled/ eheader congeladodospesos. Legadointacto. run.py aceita
 operator_policy eusa temporáriosprivados `.scheduled.POLICY`; método final
 `scheduled_MDM_POLICY`. Baselinepadrão nãoescalonado permanece opcionaldefault.
- Políticas original/static/size/generic/joint/adaptive mantêm15operadores,
 reordenam apósmelhoria,consomem permutaçãoRNGoriginal. Static/size/adaptive
 evitam descritoresdesnecessários;generic não calcula cobertura. Jointusa
 tempoGENERIC,isolandopredição deprobabilidade. Adaptativo sucesso/tempo com
 priorfixo1sucesso/.001s; simplescontrole, nãoalgoritmopublicadoreproduzido.
- PILOTO70iterações (NÃOtempoigual) COMPLETO sessão95789 terminou:
 `scheduled_development/protocol.16776a180cf7.json`,18runs,3geometriasDEV,
 seed42. Adapterchecksanterioresforamdiferentesruns,semusar seus tempos.
 Auditoria `audit_scheduled.py`PASS18; custos2650/5186/20466emtodaspolíticas.
 Jointmenortempoobservado nas3, mas diferençaspequenas frentecontroles. Tabela
 completa emresearch_notes/operator_scheduling_protocol.md. Original5.516/7.503/
3.196s;joint4.101/5.101/2.393;generic4.478/5.327/2.631;adaptive4.120/5.542/2.481.
- LIMITAÇÃOIMPORTANTE:wrappercronometra ecalcula custosantes/depois paraTODAS
 políticas; podeinflar comparação original. Próximobuilddeveevitarwrapper para
 políticassemadaptação e incluir controle MDMcorrigidolimpo. Nãoalegar speedup
 final. Piloto é3geometrias1semente,semtesteestatístico e semtempoigual.
- Nadaativo desta retomada. Builds31132/64542/7269 echeck80450 encerrados.
 Últimolog `/tmp/pcctp_scheduled_pilot.log`. Include/Utils.hpreservado.

Próximo concreto: instrumentarprogressodentroda busca local (melhor aceito),
 adicionar prazo comumcooperativoregistrandoultrapassagens e construircontroles
 legacybare/original/static/size/generic/joint/adaptive. Não usar checkpoint
 esparsocomo melhoratdeadline. Removeroverheadextra deoperaçãonaspolíticasnão
 adaptativas. Sóentão ampliar sementes/geometriasDEV. Congelartudopré-reservado.
 Preservar pilotoatual semsobrescrever. Sinalpreliminar,end-to-end nãoéainda
 contribuição pronta; novidadegenérica limitada; artigoesqueleto,CTPnãoimplementar.


## Retomada às 6h51 — coleta final e inferência nativa verificada

- Coleta093590aaa738 TERMINOU:20geometriasTRAIN planejadas,semfalhasdegeração.
 Sessão26945 encerrada. Total no banco23perfis/22geometrias(2DEV+20TRAIN),
36334chamadas;21perfisTRAIN/32808chamadas porvarianteextra degr229.
10perfis checkpoint-only e12capados em2000chamadas. Nãoalegar70iterações
 para interrompidos. Todasasrotas/custos antes/depois eproveniênciaauditados.
- LOGO completo todosTRAIN terminou: `operator_prediction.eb53cc02a54c.json`
 compara3modelos. Acrescentadaablação genérica que retira quota,prêmio,cobertura
 mantendo tamanho/seleção/obrigatórios/custo/variação/economia. Resultado4modelos
 em `operator_prediction.830f5d4ca7e1.json` (sessão29524 encerrada).
 Brier operador_only0.117835,size0.116365,generic0.114996,full0.113794.
 Logloss0.377717/0.372386/0.368022/0.364934. LogtempoMAE
1.272162/0.468763/0.456017/0.470808. FullvencegenericemBrier14/20geometrias;
 sinalpequeno,semtesteconfirmatório,semvalor de política demonstrado.
- Implementados `export_operator_models.py` e `operator_inference.cpp`:
 ajustes finais SOMENTE20TRAIN,mesmos pesosporgeometria e penalidadesfixas.
 Modelos congelados em `operator_models/models.52b8e9ae6089.json`.
 Verificação nativa vsPython102vetores(inclclipping),15ações,4modelos:
6120comparaçõesPASS,erro máximo5.33e-15. Custo viaPython/ctypes20–29us,
 EXCLUIextração de estado da rota. FonteexataPY/CPP arquivada porSHA.
 Sessões15470/13089encerradas; log `/tmp/pcctp_operator_models_ablation.log`.
- ABIATUAL0=operator_only,1=operator_size,2=operator_generic,3=operator_state.
 Modeloantigo3cabeças bfc2363c359b tinhaABIíndice2=full comCPPantigoarquivado;
 NÃOusá-lo comCPPnovo. Usar52b8e9ae6089 eíndicesexplícitos parafuturopiloto.
- Ainda NÃO implementados extração nativa dosdescritores nem política noMDM.
 O modelo prediz ocorrência demelhoria, não magnitude. Umrankprob/tempo é
 hipóteseagoraautorizada para piloto, não conclusão. Compararfullvsgeneric com
 MESMAcabeçadetempo genérica paraisolar cobertura/prêmio;includeinferência.
- Nota `research_notes/operator_scheduling_protocol.md` contém tabela,limitações
 einstruções donovopiloto. Semprocessosativosdesta retomada.

Próximo concreto: implementareverificar extraçãonativa contraPython nasrotas
 exportadas; criarbuildMDMescalonadoseparado comtodos15operadores (semomissão),
 reordenando somenteapós mudança de estado. Controlesordemoriginal,prioridade
 estáticaporoperador,tamanho,generic,full eadaptaçãosimples; usarinstânciasDEV,
 mesmosorçamentostotais/sementes. Nenhum reservado aberto. Se nãohaverganho
 preservarrefutação. Novidadegenérica limitada pelosantecedentesjáregistrados.
 Objetivofinalnãoatingido; artigo continuaesqueleto,CTPforadoescopo.


## Retomada às 6h41 — controle de tamanho e antecedentes diretos

- Coleta ampliada093590aaa738 ainda ATIVA sessão26945, log
 `/tmp/pcctp_move_profile_expanded.log`; não duplicar. Naúltimaleitura avança
 alémdord100; planos20TRAIN, timeout90s/caso. Há checkpoints interrompidos,
 preservados explicitamente. Nenhum sinal de falha do coletor.
- `operator_prediction.py --all-profiles --partial` agora permite snapshot
 explicitamentePARCIAL. Sem--partial, avaliaçãoampliada exigecompletionmanifest.
 Corrigida saída versionada parahashdeFONTES+CODIGO, evitando sobrescrever
 análises comcontrolesdistintos. Fonteexata salva emoperator_prediction_sources/
 porSHAantesdeavaliar. Fontesdosdoissnapshots anteriores recuperadas, hashes
 verificados. Default4geometriasrerodado comnovo controle/códigopreservado.
- Dois snapshots comparáveis,13geometriasTRAIN, fontesiguais, em
 `operator_prediction.046ea1c038fb.partial_snapshot.json`(2modelos) e
 `operator_prediction.84a3044c9ba7.partial_snapshot.json`(3modelos).
 Controleadicional operator_size usaSÓlog(nvisitas),log(tamanhorota), poroperador.
 Resultado: logtempoMAE operator_only1.21728,size0.45310,estado0.45626.
 Portanto ganhoaparente frente amédia simples é explicadoportamanho;estado
 nãoadiciona melhora dotempo. Brier0.10991/0.10841/0.10783;logloss
0.35927/0.35286/0.35083. Pequenosganhos preditivosnãoequivalemaganho de política.
 Sessões81380/82833/72567 encerradas. Nãoimplementar política antesdeexaminar
 coletaampliadacompleta e essescontrolesfortes.
- `operator_opportunity.py`: diagnóstico dosprefixosfalhos antesdeumoperador
 observado comsucesso, mantendoamesmarotadeentrada. Snapshot17geometrias
 (15TRAIN+2DEV) tevefraçãomédiaporgeometria0.60492dotempo nativo emtaisfalhas.
 Isto sugere espaço para testarordenação,mas NÃO ébounddespeedup ouganhode
 trajetória: outraordem mudaestadosfuturos. Excluiconstrução/overhead,censura
 mantida. Fontes/códigos congelados emoperator_opportunity.c4e78bbe1ddf.json.
- Literatura primáriaverificada:`research_notes/operator_scheduling_prior_art.md`.
 Dudel2026(JMS,doi10.5282/jums/v11i1pp1-26) jáestuda seleçãoaprendida sensível
 atempoemALNS; Zou,Hao,Wu2024 arXiv2403.14405 járeordena RVNDcomRL.
 Nãoalegarnovidadegenérica dessascombinações. Páginaeditorial/páginasDudel
 têmmetadadodepaginaçãoambíguoa conferirantesdoreferencialfinal.

Próximo: concluir eavaliartodos20TRAIN usandooperator_only,operator_size,
 operator_state,coma mesmaavaliaçãoporgeometria. Coletorraiz aindaativo.
 Casohajasinal contextual, testar política completa semomitir operadores,
 contraordemoriginal,estáticaportempo/ganho eadaptaçãosimples. Nenhumteste
 reservadoaberto; nenhumCTP; semvantagem/novidade demonstrada,artigo não pronto.


## Retomada às 6h31 — previsão observacional e ampliação TRAIN

- Coleta inicial deoperadores terminou (sessão78060):6perfis,11526chamadas
 auditadas antes/depois. gr229 é checkpoint-only no limite90s,custo61427;
 não misturar combaseline70iterações61278. Outros5perfis completaram e mantiveram
 custo final. Arquivos em `experiments/original_mdm/move_profiles/`.
- Implementado `operator_prediction.py`: auditoria de todas as rotas,custos,
 proveniência e papéis; previsão LOGO SOMENTE4geometriasTRAIN/8000observações/
1524melhorias. DEV auditado mas não treinado. Logistic por operador vs operador
 com10descritores deestado e interações. Penalidadesfixas e pesosporgeometria,
 padronização sótreino. Ridge prediz logtempo; piso1microssegundo explícito.
- Resultado `operator_prediction_development.json`: estado pioraBrier
0.10585→0.11334 e logloss0.35150→0.38015. Erro logtempo0.66310→0.56080,
 melhora2/4geometrias. Semvalor depolítica/speedclaim. ExtraçãoPython mediana
65–206us/estadoúnico,contada em auditoriadefeatures; nãoignorarna inferência.
 Sessões67328/11802 encerradas. Log `/tmp/pcctp_operator_prediction.log`.
- Perfil ampliado ATIVO sessão26945:
 `python3 experiments/original_mdm/profile_moves.py --train-limit 20 --dev-limit 0
 --include-incomplete`,log `/tmp/pcctp_move_profile_expanded.log`.
 Protocolo093590aaa738 (20TRAIN somente); timeout90s/caso, prefixo2000chamadas.
 Reutiliza mesmainstância/semente/build jávalidada. Primeiraa280 gera registros
 (raw.moves.jsonl cresceu); ainda não finalizou. Não duplicar processo.
- Coletor agora registra falhas RuntimeError emmanifest econtinua,sem excluir
 casos do plano. Não tenta salvar como completos os casos semcheckpoint.
 Fontes do protocolo ehashbinário congelados. Perfisgrandes podemtomar~30min.
- Avaliador tem `--all-profiles`, saída versionada porhash dos arquivos. Agrupa
 TODASas variantes damesmageometria fora dotreino emLOGO e ponderapesospor
 geometria,porque ampliação podeacrescentar variante diferente deuma jáperfilada.
 Não rodar avaliação final daampliação enquanto coleta ainda está ativa; snapshot
 parcial deve ser explicitamente descrito. Default continua quatrogeometrias.
- Nota `research_notes/operator_scheduling_protocol.md` cominterpretação, limites
 observacionais, censura e próximoscontroles. Nãohámodelo nativo/nova política
 implementada; melhoriaoffline nãoequivale acontribuição científica.

Próximo: acompanharcoleta semduplicar, auditar novosperfis, e rodar
 `OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3
 experiments/original_mdm/operator_prediction.py --all-profiles` ao terminar.
 Antes dereordenarMDM, verificarseprevisão temsinal contra médiaoperador e
 sóentãocomparar políticas completas comtodosoperadores,tempo total econtroles.
 Preservar falhas/resultadosnegativos; teste reservado intacto; nãoimplementarCTP.
 Artigo continuaesqueleto e objetivo finalnãoatingido.


## Retomada às 6h21 — blocos refutados; perfil de operadores em andamento

- Novos controles em `dp_trajectory.py`: bloco contíguo8visitados+4ausentes,
 seleção por proxy de economia e controle8visitadosdispersos+4ausentes. MesmaDP,
 cap18,0.2s/consulta,5s online/seisgeometrias. Protocolo780bd67687d6;18trajetórias.
 Sessão57189 terminou; auditoria passou em42trajetórias acumuladas.
- Nenhum superou checkpoint MDM anterior ao endpoint. Bloco uniforme empatou2/
 perdeu4,bloco por economia3/3,disperso4/2. Blocos fizeram22738/21638consultas
 sem limites;disperso764com143limites. Resolver mais consultas baratas não deu
 melhor qualidade. Preservar a refutação e não promover redução de q como ganho.
- Nova hipótese observacional: aprendizado para ordenar os operadores existentes
 doMDM, onde há mais transições com melhoria que no rankingDP. Protocolo em
 `research_notes/operator_scheduling_protocol.md`; ainda NÃO há política aprendida.
- Build isolado `build_move_profile/` via `build.py --initialize-decoder-index
 --profile-moves`, sem alterar fontes legadas, adiciona wrapper de tempo/rotas
 sobre `obter_vizinho`. Modelo de decisões e RNG intactos. Outputs só primeira
2000chamadas/instância; censura explícita, custo nativo exclui exportação.
 Profiling tem overhead: não comparar seu walltime ao baseline.
- NOVO processo ATIVO sessão78060: `python3 experiments/original_mdm/profile_moves.py`,
 log `/tmp/pcctp_move_profile.log`. DuasgeometriasDEV +quatroTRAIN,90s/caso,
 uma variante e semente42. Não duplicar. Build sessão52333 já terminou.
- Passaram att48(1526movimentos,semcensura),berlin52(2000,censurado),bier127(2000,
 censurado), com rotas/custos antes/depois recalculados e checkpoints comuns
 iguais aos baselines. Maisgeometriasdetreino ainda em execução.
- Protocolo4e1c42c7aa48 em `experiments/original_mdm/move_profiles/`, com sidecar
 `provenance.4e1c42c7aa48.json` que registra hashes dasfontes ebinário capturados
 durante coleta sem modificações posteriores. Raw ebuild ignorados noGit.
- run.py ganhou `profile_moves=False` mantendo padrão de algoritmointacto; perfis
 usam diretórios temporáriosprivados `.profile` ebináriosFPMaxisolados.

Próximo: terminar e auditar perfis, resumir ganho/tempo por operador PORgeometria.
Se houver sinal,predizerprobabilidade/ganho/tempo em TRAIN comcontrole médiaspor
operador;dadosobservacionais não são contrafactuais. Não saltar de acurácia para
alegação de melhoria; uma política será medida em trajetórias completas, com
orçamento total econtroles. NãoimplementarCTP. Artigo ainda não pronto.


## Retomada às 6h11 — busca DP completa e novos rótulos

- Coleta corrigida das58outrasvariantes terminou sem falhas de geração;
  timeout/checkpoint não significa conclusão das70iterações. Auditoria passou:
 86baselines,57completos,29checkpoint-only,440snapshots,104estados/2226consultas.
  Nada de validação/teste reservado foi rotulado. Sessões46473/72399 encerradas.
- `dp_trajectory.py`:24trajetórias, seis geometrias, quatro controles,5s online,
  limite18componentes e0.2s/consulta. Protocolo66e93b8c4d05 em
  `dp_trajectory_development/`. Sessão59026 terminou; auditoria independente
  `audit_dp_trajectories.py` passou em24. Não há processos ativos desta retomada.
- Nenhum controle superou checkpoint MDM anterior ao endpoint. Empates/piores:
  optional_uniform3/3, anchor_uniform4/2, anchor_alns4/2, anchor_random5/1.
  Em berlin52, âncoras próximas5244 versus aleatórias5186=MDM; em eil51 todos276
  versusMDM273. Não confundir vantagem contra controle restrito com ganho sobreMDM.
- Implementados `dp_labels.py` e `evaluate_dp_labels.py`. Nova coleta apenas em
 20geometrias TRAIN, primeira variante/fase1,12grupos opcionais x2estratégias de
 âncora(próxima/aleatória),480consultas. Fontes congeladas, trava, publicação atômica.
 123consultas excederamcap18e permanecem com motivo;24melhoraram. Auditoria de todas
 rotas,objetivos,ganhos,footprints eproveniência passou em `dp_training_audit.json`.
 Sessão22735 terminou. Log `/tmp/pcctp_dp_training_labels.log`.
- LOGO exploratório APENASentregeometriasdetreino: ganho médio normalizado por
 geometria:uniform0.0000711,capped_uniform0.0001024,regra menorq0.0004085,
 capped_joint_saving0.0002386; quatro modelos(linear/RFF,estrutural/conjunto)
 todos0.0003639. Portanto nenhum ML bateu o controlemenorq. Alvos escassos:
 apenas7das20geometrias têm algum ganho na amostra. Não é confirmação nem
 avaliação de trajetória aprendida. Outputs `dp_ranking_training_results.json`.
- Protocolo e limitações `research_notes/dp_learning_protocol.md`. Rótulos antigos
 MIP/route_context não foram reutilizados. Número de componentes ecap são dados
 tantoaMLquantoacontroles; não alegar valor de ML por evitar consultasimpossíveis.

Próximo passo: não insistir em ajustar o ranker sob rótulos tão esparsos. Investigar
se a decisão com sinal é liberar âncoras/conjuntos maiores de forma contígua para
controlar q, ou usar ML em construção/continuação do MDM. Antes de nova coleta,
formular hipótese e controle simples correspondente. Preservar todos negativos.
Uma avaliação adicional possível é estadofase10 de TRAIN, mantendo desenvolvimento
como diagnóstico e teste intocado. Não confundir efeito de DP/clique com aprendizado.
Artigo ainda esqueleto; sem vantagem/novidade suficiente, sem conclusão da tarefa.


## Retomada às 6h00 — nova formulação do reparo e DP

- Piloto anneal terminou. Continua sem vantagem de ML demonstrada. O diagnóstico
  de âncoras também não chegou ao ótimo do MDM em berlin52/eil51. Não insistir
  em aumentar parâmetros sem uma hipótese nova e controles.
- Corrigido `free_y` dos nove arquivos em `routing_release_development/` para
  contar apenas opcionais livres; os quatro obrigatórios de roteamento ficam
  selecionados. Metadado corrigido documentado; tempos/rotas/custos intactos.
- Implementados `path_dp.cpp` e `dp_repair.py`: DP exata sobre componentes de
  caminhos fixos, orientações e máscaras de seleção viável. Backend experimental
  usa Held–Karp clássico; não alegar que a DP ou contração são inéditas.
- Primeiro diagnóstico com `route_context`:108/108casos excederam18componentes,
  registrado em `dp_development/`. O excesso vinha da liberação de muitas
  arestas incumbentes, que destruía os caminhos a contrair.
- Novo operador `fragment_clique` em `pilot.free_edge_set`: remove arestas
  incidentes na liberação J, mantém os caminhos restantes e libera a fronteira
  completa, exceto as arestas internas fixas. Se J só tem ausentes, abre um
  intervalo do incumbente. Tem q<=2|J| componentes ativos. Argumento e recorrência
  estão em `research_notes/fragment_clique_dp.md`, ainda para revisão/novidade.
- Importante: a primeira versão libertava por engano a aresta interna de
  caminhos de tamanho2. Corrigido; outputs anteriores segregados em
  `fragment_clique_pre_boundary_fix/` e não usados na avaliação atual.
- `verify_dp.py` passou:30casos contra SCF completo,18vizinhanças por enumeração,
  timeout e retorno viável. O teste revela que checar o relógio pelos bits da
  máscara não funciona quando a raiz está no bit zero; isso já foi corrigido.
- Em108vizinhanças pareadas com `fragment_clique`, DP certificou106 e MIP104;
  duas excederam o limite de18componentes. Resultados em
  `fragment_clique_development/` e `fragment_clique_vectorized_development/`.
  Vetorização da seleção baixou mediana observada de38ms para17ms, contra36ms
  do MIP no lote pareado. Ainda é diagnóstico sob carga, sem claim final de speedup.
- Lotes DP e verificações terminaram; última checagem pós-otimização em
  `/tmp/pcctp_dp_final_verification.log` (sessão7962). Verificar conclusão.
  O caminho `dp_build/` é ignorado no Git; `library()` compila a dylib local
  com g++-15 na primeira utilização. Compilação não deve ser contada como
  inferência, mas construção/carregamento do modelo online deve ser medido.
- Coleta das outras variantes corrigidas (sessão46473) ainda avança, log
  `/tmp/pcctp_decoder_initialized_other_variants.log`. Primeiras variantes já
  encerradas; rótulos temporizados podem ser coletados incrementalmente, com trava.

Próxima hipótese: separar seleção opcional U de âncoras obrigatórias R e escolher
vizinhanças pelo ganho esperado e pelo tamanho combinatório q, usando a mesma
fronteira completa e o mesmo limite DP nos modelos e controles. Antes de ML,
testar algoritmo completo com essa reparação e controles sem aprendizado.
O ranking antigo treinado em MIP/route_context não é automaticamente validado
para o novo operador; coletar novos alvos ou declarar explicitamente transferência.
Controlar casos com limite de tamanho e NÃO omiti-los.

## Retomada às 5h50 — resultado negativo e novas verificações

- O piloto greedy de6geometrias/dezpolíticas terminou: nenhuma política venceu
  o MDM nos checkpoints de tempo registrados. O ranking conjunto não superou
  os controles pareados. Preservar esse resultado; não promover o ganho contra
  o estado inicial como superioridade frente à continuação do MDM.
- `audit_trajectories.py` passou em105trajetórias disponíveis: custo final e
  viabilidade recalculados, melhor custo monotônico, tempos ordenados, modelos
  sem a geometria avaliada no treino. Protocolos podem estar incompletos no
  instante da auditoria. Rodar novamente após cada lote.
- Foi acrescentada aceitação `anneal`, de forma igual para todos os métodos.
  Temperatura inicial1%do custo inicial, final0.01%, decaimento exponencial por
  tempo. A construção pode propor pioras; retorna-se sempre a melhor solução
  verificada, separada da solução corrente. Greedy permanece disponível.
  Não houve ajuste em validação/teste nem escolha de temperatura por resultado.
- Sessão63247: pilotoanneal seisgeometrias, log
  `/tmp/pcctp_hybrid_anneal_six_geometries_round2.log`, ainda em execução.
  A primeira tentativa falhou na serialização de um booleano NumPy antes de
  publicar qualquer resultado; corrigida e reiniciada. Não é falha do solver.
- Protocolos de futuras execuções agora congelam também os nomes e hashes dos
  casos de avaliação antes das trajetórias. Não refazer o lote só por esse
  acréscimo de rastreabilidade. O código em execução já foi carregado antes dele.

### Diagnóstico de liberdade de roteamento

Seleção condicional completa não implica bom espaço de roteamento: as arestas
entre obrigatórios podem continuar presas, apesar de toda seleção viável ter
alguma rota na vizinhança. Por isso começou um diagnóstico específico:

- `routing_release_probe.py`, sessão79330, log
  `/tmp/pcctp_routing_release_probe.log`. Três geometrias de desenvolvimento,
  primeiro estado do MDM,12consultas de0.5s por configuração.
- Compara8opcionais,12opcionais e8opcionais+4obrigatórios usados apenas como
  âncoras de roteamento. A visita dos obrigatórios continua fixada emum; só
  suas arestas incidentes são liberadas. Os dois últimos casos têm o mesmo
  número total de vértices liberados nas arestas e o mesmo orçamento adicional.
- Isso é diagnóstico de mecanismo, não vantagem de ML ou teste confirmatório.
  ATENÇÃO: o processo foi carregado quando `free_y` ainda contava todos os
  vértices liberados. Nos JSONs desse primeiro lote, derivar o número real de
  visitas livres de `optional_released`. Corrigir esse metadado após o lote
  (sem alterar soluções/tempos). A classe foi corrigida para os próximos lotes.

Próxima decisão: se liberar apenas opcionais não altera o roteamento necessário,
estudar vizinhanças que separem seleção e reordenação de obrigatórios, mantendo
controles de orçamento e capacidade iguais. Se ML continuar sem benefício,
investigar outra tarefa de aprendizado em vez de multiplicar configurações até
encontrar uma comparação favorável. A teoria de seleção ainda requer revisão
e confronto com literatura; não é novidade estabelecida.

## Retomada às 5h40 — busca iterativa e controles equivalentes

- `hybrid_lns.py` implementa trajetórias completas de desenvolvimento com
  custo de seleção/descritores/reparo incluído, controle adaptativo simples
  sem ML, ranking genérico e conjunto. A variante híbrida faz MIP a cada8consultas,
  até0.2s, com a solução construtiva como ponto inicial viável. Há híbridos
  aleatório, regra simples, adaptativo e aprendido: não atribuir ao ML o efeito
  de acrescentar MIP. ALNS aqui é um controle simples, não réplica de um método
  publicado. Aceitação atual é apenas de melhora, sem simulated annealing.
- `FragmentRepair` agora aceita `start_route` validada contra todas as fixações,
  conserva a melhor solução inicial em timeout e usa-a como MIP start nativo.
  `verify_native.py`:60casos comparados com SCF completo, todos passaram com
  soluções construtivas iniciais;32cortes foram adicionados nesses testes.
- Pilotos `/tmp/pcctp_hybrid_lns_pilot.log` (sessão33352) e
  `/tmp/pcctp_hybrid_matched_controls.log` (sessão84230) terminaram. ATENÇÃO:
  limit3 nesse piloto significou3instâncias da MESMA geometria att48, não3
  geometrias independentes. Nada de inferência estatística a partir disso.
- No primeiro caso att48, métodos construtivos sem MIP mantiveram2692;
  híbridos aleatório/regra/genérico/conjunto chegaram2650 em5s extras. Isso
  mostra utilidade do reparo exato nesse exemplo, não vantagem específica do ML.
  O MDM também atingiu2650 no próprio histórico. Demais variantes não deram
  uma vantagem geral; resultados completos nos JSONs.
- NOVA sessão37277 em execução:
  `hybrid_lns.py --budget 5 --limit 6 --one-per-geometry`, log
  `/tmp/pcctp_hybrid_six_geometries.log`. Agora são até6geometrias diferentes,
  primeira variante disponível de cada, dez políticas. Ordem dos métodos é
  aleatorizada deterministicamente por instância/semente. Protocolo salva hashes
  do código e dados; modelos treinados só em geometrias de treino e congelados
  antes das trajetórias. Não iniciar cópia desse piloto em paralelo.
- Sessão10330: coleta incremental com tempos corrigidos, log
  `/tmp/pcctp_timed_fixed_oracle2.log`, diretório `timed_oracle_decoder_initialized_states/`.
  Verificar término antes de reexecutar. Sessão46473 das outras variantes
  corrigidas continua; primeiras variantes corrigidas terminaram com um timeout
  sem primeira solução (rd400), explicitamente mantido no manifesto.

Próximos passos: auditar trajetórias e curvas por tempo; ampliar geometrias e
sementes; investigar diversificação/aceitação controlada, pois melhorias apenas
monótonas podem ficar presas. Usar os mesmos mecanismos e orçamento nos controles
para isolar ML. A coleta concorrente serve a desenvolvimento, não a uma alegação
final de aceleração. Continua faltando comparação final controlada e novidade
bibliográfica estabelecida.

## Retomada às 5h30 — lema e reparação construtiva

Foi implementada e verificada uma construção de rota para qualquer atribuição
condicional viável de visitas, dentro do operador `route_context`:

- `selection_lift.py`, `verify_selection_lift.py` e prova para revisão em
  `selection_lifting_lemma.md`.
- Orçamento de k arestas extras (em vez de2k) é suficiente; uma família
  construtiva de exemplos exige k no pior caso. Não afirmar novidade sem
  revisão bibliográfica. Foram verificadas1954 construções e o limite inferior
  por enumeração de ciclos para k=1..6.
- `best_lift` enumera a seleção condicional e reconecta de forma heurística as
  três melhores atribuições por substituto aditivo. O substituto não é bound
  do custo, e não há certificado de reparação ótima para essa heurística.
- `probe_lifting.py`, sessão85390, terminou: em20estados/480vizinhanças corrigidas
  com tempos, o MIP limitado melhorou64consultas, a construção melhorou69;
  igualou/superou52das64melhoras do MIP e foi estritamente melhor em14consultas.
  Nenhuma construção ficou abaixo de um ótimo local certificado pelo MIP
  (checagem obrigatória no script). Ainda é diagnóstico, não algoritmo completo.
- Triagem temporal:8das64melhoras locais já eram igualadas pelo MDM em checkpoint
  anterior ao tempo otimista concedido à consulta. As outras não são vitórias
  temporais provadas; checkpoints são esparsos e os tempos históricos.
- `lifting_prior_art.md` registra antecedentes relevantes: HDR(2023) e DRHG,
  AAAI2025. Contração de caminhos e aprendizado com contração já têm antecedentes.

Próximo passo: testar a construção como controle barato dentro de busca
iterativa, e como solução inicial viável para o MIP; incluir todo o custo de
seleção. Comparar qualidade por tempo com MDM continuado e ALNS antes de atribuir
ganho a ML. Não abandonar nem omitir o fato de a regra simples ser competitiva.

As coletas corrigidas das primeiras e demais variantes ainda avançam
(`/tmp/pcctp_decoder_initialized_states.log`,
`/tmp/pcctp_decoder_initialized_other_variants.log`). Coletor temporizado1
terminou com20estados. Pode retomar incrementalmente na mesma saída, cuja trava
evita duplicações. Teste e validação continuam reservados por geometria.

## Retomada às 5h20 — processos e auditoria temporal

- Sessão87177: primeiras variantes corrigidas, log
  `/tmp/pcctp_decoder_initialized_states.log`, ainda em execução.
- Sessão46473: outras58variantes corrigidas, log
  `/tmp/pcctp_decoder_initialized_other_variants.log`, ainda em execução.
  Casos diferentes dos da primeira coleta; não iniciar duplicatas.
- Sessão83335 do coletor corrigido terminou com13estados/312rótulos.
  Auditoria específica passou em8baselines corrigidos e312reparos, sem mistura
  de métodos. Usar `audit_dataset.py --state-dir .../decoder_initialized_states`
  para a referência corrigida; a auditoria padrão ainda trata o legado original.
- Novo coletor **com tempos por consulta** em execução: log
  `/tmp/pcctp_timed_fixed_oracle1.log`, saída
  `experiments/ml_lns/timed_oracle_decoder_initialized_states/`.
  Não misturar com o coletor anterior: rótulos schema1 não possuem todos os tempos.
  `strong_oracle.py` agora usa trava de arquivo por diretório para impedir
  dois coletores simultâneos na mesma saída.
- `time_dominance_audit.py` implementa uma triagem conservadora: concede a uma
  consulta MIP apenas seu próprio tempo, omitindo otimisticamente custo de
  selecionar entre candidatos/modelo; se o MDM já atingiu custo igual ou melhor
  em checkpoint anterior, a melhora local não demonstra vantagem temporal.
  Ausência de dominância registrada não prova vantagem: checkpoints são esparsos,
  tempos vêm de execuções separadas e a máquina tem carga concorrente.
  Executar com `--oracle-dir experiments/ml_lns/timed_oracle_decoder_initialized_states`.
  Ainda faltam comparações controladas do algoritmo completo e ALNS.

As alterações desta retomada preparam uma comparação justa de tempo e impedem
que melhoras sobre um estado inicial sejam confundidas com aceleração científica.

## Retomada às 5h11 — atualização prioritária

Esta seção substitui o estado antigo dos processos abaixo.

- As coletas originais terminaram: 29 primeiras variantes, um timeout sem
  primeira solução; 58 outras variantes, três SIGSEGV. As sementes 43/44 nas
  nove primeiras variantes de desenvolvimento terminaram com duas SIGSEGV.
- Auditoria `experiments/ml_lns/dataset_audit.json`: 99 baselines, 509 snapshots,
  151 estados rotulados e 3384 reparos verificados; 29 geometrias; nenhum rótulo
  reservado de validação/teste. Arquivos têm hashes. Falhas fora do corpus
  continuam nos manifestos, não foram convertidas em soluções válidas.
- `strong_oracle.py` rodada 5 (sessão 87985) terminou. Coletou as variantes
  disponíveis de seed42. Sementes43/44 ainda precisam de coleta específica.
- `nonlinear_ranker.py` e `evaluate_nonlinear.py` implementam controles lineares
  e random Fourier features (largura64, alpha1, seed2026), escala apenas no treino,
  ponderação igual por geometria informativa. Resultado ainda exploratório,
  sem vantagem assegurada frente à regra `joint_savings`. Não usar teste final.
- Exportações agora são publicadas atomicamente somente após validação. Novas
  execuções exportam tempo por snapshot e tempo total do processo. NÃO imputar
  tempos aos snapshots antigos que não têm esse campo.

### Falha do legado localizada e correção isolada

O decodificador em `mineracao.cpp` declara `int cc_atual;` sem inicializar e o
usa ao criar componente de padrão com uma única aresta. Relatórios locais de
crash apontaram `decode`; teste mínimo `verify_decoder.py` passa com inicialização
em zero. A versão corrigida é separada, nunca uma substituição silenciosa:

```sh
python3 experiments/original_mdm/build.py --initialize-decoder-index
python3 experiments/original_mdm/verify_decoder.py
```

Diretório `build_decoder_initialized/`, método JSON
`original_MDM_legacy_decoder_index_initialized`. Fontes fornecidos intactos.
Os três casos de seed42 que antes falhavam concluíram e passaram na validação:
eil51_VT9_T11_W31_50 custo248, hk48_VT9_T10_W29_75 custo7803,
rd100_VT20_T60_W20_50 custo6626. Controle brazil preservou20466.
Resultados em `experiments/original_mdm/diagnostics/decoder_initialized_reproduction.json`.
Não inferir daí equivalência universal com execuções de comportamento indefinido.

ASAN/UBSAN com clang concluiu dois casos sem reportar erro e atingiu timeout no
terceiro; isso não invalida o defeito de inicialização nem certifica ausência
de outros bugs. GCC não dispõe de libasan local. LLDB foi encerrado depois da
identificação por inspeção e teste mínimo; não há depurador útil a manter.

### Processo atualmente em execução

`strong_states.py --roles development train --seeds 42 --initialize-decoder-index`,
log `/tmp/pcctp_decoder_initialized_states.log`. Coleta as primeiras variantes
com a correção EXPLICITAMENTE identificada; arquivos em `decoder_initialized_states/`.
Não misturar com `original_states/` nem reaproveitar seus rótulos como se fossem
do mesmo algoritmo. Split por geometria foi copiado, preservando reservas.

O coletor aceita banco separado e escolhe saída separada automaticamente:

```sh
python3 experiments/ml_lns/strong_oracle.py --state-dir experiments/ml_lns/decoder_initialized_states
```

Antes de comparar algoritmos novos, estabilizar essa referência e verificar
consistência da correção nas outras instâncias. O teste final continua intocado.
Próximo estudo necessário: qualidade por tempo desde o primeiro estado do MDM,
comparado à continuação do próprio MDM, ALNS sem aprendizado e regra simples.

Atualizado em 18 de outubro de 2026. Não é um relatório de contribuição concluída.

## Decisões vigentes

- Branch `codex/pcctp-learning-submission`; commits `f0cc6a6` e `ac048e4`
  registram memória, esqueleto editorial e pesquisa acumulada.
- Usuário pediu explicitamente não aguardar aprovação e prosseguir até paper
  completo/revisado. Autorizou retomadas a cada 10 minutos durante sua ausência.
- Heartbeat ativo `avan-ar-pesquisa-e-artigo-pcctp`, a cada 10 minutos neste chat.
  Não criar duplicata. Desativar apenas ao cumprir os critérios ou pedido do usuário.
- CTP continua fora do escopo atual. `include/Utils.h` é modificação preexistente
  do usuário, não incluir em commits da pesquisa.

## Processos disparados nesta etapa

Os IDs abaixo são sessões de exec da conversa; se indisponíveis na retomada,
inspecione os arquivos de saída e logs antes de iniciar outra execução.

- 25250: `strong_states.py --roles development train --seeds 42`, retomada dos
  casos inicialmente interrompidos. Log `/tmp/pcctp_original_states_round2.log`.
  Arquivos completos ou checkpoints validados em `experiments/ml_lns/original_states/`.
  Manifesto final `last_generation.json`.
- 1342: `strong_states.py --roles development train --seeds 42 --other-variants`.
  Log `/tmp/pcctp_original_other_variants.log`; duas execuções paralelas.
  Coleta TODAS as outras variantes das geometrias permitidas, não apenas as
  que pareçam favoráveis. Manifesto `other_variant_generation.json` ao terminar.
  Não colide com a primeira seleção: instâncias possuem nomes distintos.
- 63842: `strong_oracle.py`, log `/tmp/pcctp_strong_oracle_round3.log`.
  Só enxerga os arquivos disponíveis no início. Reexecutar incrementalmente
  quando terminarem novas instâncias, sem iniciar dois coletores para os mesmos
  estados simultaneamente. Arquivos `.partial` não são resultados completos.
- 40974: diagnóstico de raios com backend highspy, warm start e cortes externos;
  log `/tmp/pcctp_radius_cuts.log`. Provavelmente já concluído na retomada.
- Rodadas anteriores de raio com scipy, highspy frio/quente e busca iterativa
  de 3 geometrias são curtas e já produziram seus JSONs. Não repetir sem motivo.

## Verificações realizadas

- `verify_signature.py`: 60 casos aleatórios, 960 atribuições, dez descritores.
- `verify_fpmax.py`: padrões maximais de seis limiares coincidiram com enumeração.
- `checkpoint_verification.json`: mesma semente preservou rota e custo original;
  timeout artificial reteve solução viável validada independentemente.
- `verify_fragment.py`: 37 casos, incluindo 15 por enumeração e 20 contra SCF
  completo; dois casos liberando apenas vértices ausentes.
- `verify_native.py`: reproduz 60 comparações da formulação com cortes externos
  e highspy, com/sem solução inicial; resultado salvo passou, 35 cortes adicionados.
- `fixed_path_repair.md`: argumento de equivalência inteira e limite de até 4k
  componentes para liberação incidente em k opcionais mais 2k arestas. Revisar
  formalmente antes de entrar no artigo; não alegar novidade da contração.

## O que os resultados dizem até aqui

- Piloto anterior sem MDM mostrou vantagem do contexto de rota, mas não sustenta
  vantagem científica frente à heurística original.
- MDM original foi recuperado. A primeira coleta teve nove timeouts em 29 casos;
  agora checkpoints impedem a perda total das observações. Não excluir geometrias.
- Em `strong_development_results.json`, o ranking linear com descritores conjuntos
  supera o genérico em uma avaliação exploratória, mas a regra simples
  `joint_savings` permanece competitiva ou superior. Ainda não há valor adicional
  de ML comprovado nem significância confirmatória.
- Busca iterativa adicional de 10 s em att48, berlin52 e brazil58 após 70
  iterações do MDM: nenhuma das quatro políticas melhorou o custo.
- Raios 8/16/24 em gr137, pr144, u159 também não melhoraram nos diagnósticos.
  Nos artigos fornecidos vários desses valores coincidem com ótimos reportados;
  igualdade numérica deve ser distinguida de certificado computacional próprio.
  Nessas instâncias é relevante reduzir tempo até o ótimo, não melhorar o custo.
- Contração não assegura aceleração: tempos preliminares são mistos e rodados
  sob carga concorrente. Comparações finais devem usar execução controlada.
- Gurobi instalado tem licença expirada; continuar com HiGHS disponível sem
  aguardar usuário, comprar licença ou alterar configuração de licença.

## Próximos passos concretos

1. Verificar logs/manifestos e ausência de processos duplicados. Aproveitar as
   variantes adicionais já concluídas para coleta incremental de rótulos.
2. Auditar os resultados e medir também tempo até alvo e integral de qualidade.
   Comparar desde inicialização de mesmo orçamento, não atribuir a ML uma
   melhora que o próprio MDM alcançaria mais rapidamente continuando sua busca.
3. Ampliar rótulos em estados e variantes difíceis de treino/desenvolvimento.
   Testar modelos não lineares com os MESMOS controles de capacidade e regras
   simples. Evitar acrescentar uma rede apenas para dar aparência de novidade.
4. Integrar política completa ao fluxo GRASP/MDM; controles ALNS sem aprendizado
   e MDM sob tempo igual ainda faltam. O script `iterative_study.py` é um piloto,
   não a avaliação final. Seu modelo é congelado com hashes antes da trajetória.
5. Se ML não superar os controles, registrar refutação e investigar alternativa
   específica, mantendo validação/teste reservados e verificando literatura.
6. `submission/main.tex` é apenas esqueleto. Ainda faltam contribuição sustentada,
   protocolo final, experimento confirmatório, bibliografia auditada e texto.
   Não preencher vinte páginas com alegações não demonstradas.

## Comandos úteis

```sh
python3 experiments/ml_lns/strong_oracle.py
python3 experiments/ml_lns/evaluate_strong.py
python3 experiments/ml_lns/verify_signature.py
python3 experiments/ml_lns/verify_fragment.py
python3 experiments/ml_lns/verify_native.py
python3 experiments/ml_lns/radius_probe.py --backend highspy --warm-start --connectivity cuts
```

Antes de escrever o manuscrito: aplicar scientific-writing (o SKILL.md instalado
é um ponteiro textual para outro arquivo; localizar a fonte real), usar fontes
primárias verificadas e copiar/aplicar template editorial de forma rastreável.
Ao gerar PDF, seguir skill de PDF, compilar e inspecionar todas as páginas.
