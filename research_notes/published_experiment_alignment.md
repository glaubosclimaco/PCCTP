# Auditoria solicitada pelo usuário: tempos e protocolo dos artigos

## Evidência conferida em19/10

RAIRO2021, PDF publicado, seção4.3: maxIt=70; calibração37instâncias,
10execuções por configuração; comparação nas107restantes com10sementes
pareadas. Fonte https://www.numdam.org/item/10.1051/ro/2021057.pdf,
linhas extraídas618-625. Critério heurístico é iterações, não300segundos.

MDM, manuscrito fornecido artigos/mdm/main_review.tex: linhas397hardware
Intel i7-6700/3.40GHz/16GB;408mesmo maxIt da heurística-base;411dez
execuções no ajuste;415comparação107instâncias;837custosdezexecuções
porheurística. Linha811explicita70iterações. Fonte principal.cpp linha95
numIt=70. MDM20%suporte; legadocap10,supportabsolute2; alteraçõesde
ocupaçãoelite e decodificador seguem explicitamente declaradas.

Tabela MDM desse manuscrito: a280_R56_T168_W56_50, linha724:
MDMmean2001.40, tempo médio946.62s; rd400_R80_T240_W80_25, linha747:
MDMmean10568.30,tempo médio4641.00s. São tempos observados naquele
hardware, não limites a impor literalmente em máquina diferente.

IJISE2023, manuscrito local sample(single col).tex linha633: solver
Gurobi6.5.2, limite3600s, opções/threads informadas no texto. Esse limite
exato não é o critério de parada das heurísticas MDM.

## Estado atual e correção de interpretação

Os lotes60/300s,9hardTRAIN e3sementes NÃOreproduzem o experimento dos
artigos. São estudos exploratórios de orçamento igual entre métodos nesta
máquina. Ganho versus MDMinterrompido cedo não implica ganho versus
MDM70iterações completo, principalmente se mineração não ativar.

Jáhouve perfis70iterações:7completos e2censurados no lote fixed_work.
a280seed42recovered70iterations retornou em1721.512420s nesta máquina;
éum caso/uma semente, não tempo médio comparável aos946.62spublicados.
Não tratar os censurados como70iterações completas nem os perfiscomo
replicaçãoestatísticacom10sementes.

A decisão sobre artigo precisa de uma referênciaMDM70iterações, sementes
pareadas, mineraçãoobservada e mesma máquina. Próximaetapa de alinhamento:
primeiro garantir reprodução completa dos9casos, incluindo os2censurados;
depois10sementes previamente fixadas porcaso e candidato com orçamento
totalcomparável, debitando inicialização. Registrarprotocolo/hashesantesda
coleta. A duração de uma iteração dos métodos diferentes não é comparável;
70iterações do híbrido não igualam esforço a70iterações doMDM.

Não mudar o lote300s em andamento nem apagar resultados. Reclassificá-lo
como estudo anytime/complementar, não reprodução definitiva ou evidência
suficiente para abandonar/publicar oproblema. Não usar dadosreservados
para ajustar métodos nessa correção. Teste12continua fechado.
