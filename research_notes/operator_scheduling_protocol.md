# Próxima hipótese: escalonamento dos operadores do MDM

Estado: diagnóstico observacional; nenhum benefício de ML estabelecido.

As tentativas anteriores de ranking de vizinhanças exatas tiveram rótulos
esparsos e não superaram controles simples. A busca com blocos contíguos
também não resolveu a limitação: no piloto de seis geometrias, permitiu milhares
de reparos baratos, mas nenhum método superou os checkpoints comparáveis do
MDM. Uniforme contíguo empatou2e perdeu4; proxy de economia empatou3e perdeu3;
disperso empatou4e perdeu2. São apenas diagnósticos, sem testes reservados.

Hipótese diferente: prever o ganho e o custo de operadores do próprio MDM
durante a descida, onde existem mais transições com melhoria, para ordenar
operadores pelo retorno previsto por tempo. Não atribuir originalidade genérica
à seleção de operadores aprendida; verificar a literatura se houver sinal.

`build.py --initialize-decoder-index --profile-moves` cria uma cópia isolada em
`build_move_profile/`. O wrapper observa `obter_vizinho`: salva a rota de entrada,
cronometra somente a execução do operador e salva rota de saída/custos. Não
altera escolhas nem usa o RNG. A exportação e suas alocações têm overhead;
tempo total da versão instrumentada não serve para comparar velocidade.

`profile_moves.py` coleta as primeiras2000chamadas por instância, duas geometrias
de desenvolvimento e quatro de treino com baseline completo. Essa amostra é
censurada e concentra a fase inicial: não representa necessariamente a execução
inteira. Cada rota e custo de entrada/saída é recalculado independentemente;
checkpoints comuns precisam manter o custo do baseline corrigido.

Até agora att48, berlin52 e bier127 passaram e preservaram seus custos finais.
Nos dois primeiros, os operadores8,9,11,12 concentraram grande parte do tempo;
também produziram melhorias. Portanto não há justificativa para simplesmente
suprimi-los. Custo alto isolado não significa operador dispensável.

Os rótulos observacionais dizem o que ocorreu no operador executado, não o
ganho contrafactual dos outros operadores naquela rota. Antes de aprender uma
política, comparar previsão fora da geometria de treino com médias por operador,
regras de tamanho/estado e adaptação online. Os dados de desenvolvimento não
entram no treino. Uma primeira política poderá reordenar todos os operadores,
mantendo uma varredura completa após a última melhoria, em vez de omitir
operadores por previsão. Ganho de tempo e qualidade precisa aparecer em busca
completa sob orçamento total, com inferência medida e exploração controlada.

Não usar precisão de classificação como substituto de melhoria do algoritmo.
Não usar movimentos da mesma trajetória como réplicas independentes. Separar
geometrias, registrar mudança de distribuição após reordenação e manter os
resultados negativos. CTP continua fora do escopo.

## Primeiro teste preditivo, às 6h31

Seis perfis concluídos, 11526chamadas auditadas: duas geometrias de desenvolvimento
e quatro de treino. gr229 foi interrompido no limite de90s; a amostra de2000
chamadas é válida, mas seu custo final disponível61427não equivale às70iterações
do baseline61278. Os demais cinco perfis completaram e preservaram custo final.

`operator_prediction.py`: LOGO exclusivamente em quatro geometrias TRAIN,
8000chamadas e1524melhorias. Compara regressão logística por operador com
interceptos e interações dos mesmos operadores com dez descritores do estado.
Prediz também logtempo com regressão ridge. Padronização apenas no treino,
pesos iguais por geometria, penalização fixa sem ajuste por resultados.

O modelo com estado piorou Brier médio0.10585→0.11334 e logloss0.35150→0.38015.
Erro absoluto médio do logtempo caiu0.66310→0.56080, mas melhorou em apenas duas
das quatro geometrias. Não é demonstração robusta de benefício computacional.
A extração Python por estado único teve medianas65–206microssegundos nesse
diagnóstico; ela precisa entrar no custo da política, se implementada.

Foi iniciada ampliação para até20geometriasdetreino, incluindo baselines com
checkpoint incompleto. Usa apenas uma variante por geometria nesta coleta;
perfis anteriores de mesma instância/semente/build são reutilizados. Variantes
anteriores da mesma geometria podem existir: LOGO deve reter todas as variantes
da geometria fora do treino simultaneamente e ponderar por geometria.
Falhas sem primeiro checkpoint são registradas e não apagam geometrias do plano.
Protocolo093590aaa738. Desenvolvimento/validação/teste não entram nessa ampliação.

## Coleta ampliada e ablação, às 6h51

Coleta terminou sem falhas de geração nas20geometriasTRAIN planejadas. Incluir
perfil anterior de outra variante de gr229 produz21perfis de treino,32808chamadas.
LOGO retém todas as variantes da geometria avaliada e pondera por geometria.
Isso é avaliação exploratória no treino, não avaliação reservada.

Resultados `operator_prediction.830f5d4ca7e1.json`:

| Modelo por operador | Brier | Logloss | MAE logtempo |
|---|---:|---:|---:|
| Somente intercepto |0.117835|0.377717|1.272162|
| Tamanho da instância/rota |0.116365|0.372386|0.468763|
| Estado genérico sem cobertura/prêmio |0.114996|0.368022|0.456017|
| Estado completo |0.113794|0.364934|0.470808|

Os descritores de cobertura/prêmio melhoram um pouco a previsão de melhoria
(14de20geometrias emBrier) e pioram o tempo frente ao estado genérico. Não há
inferência confirmatória, nem justificativa para alegar ganho de busca a partir
dessa tabela. O modelo genérico preserva tamanho, fração selecionada/obrigatória,
custos das arestas e proxy de economia; retira quota, prêmio e cobertura crítica.

Modelos finais de treino congelados em `operator_models/models.52b8e9ae6089.json`,
semDEV/validação/teste no treino. `operator_inference.cpp` foi comparado com
Python em102vetores incluindo extrapolação/clipping,15ações e4modelos:
6120comparações,erro máximo5.33e-15. Chamadas viaPython/ctypes~20–29us excluem
extração dos descritores da rota, que ainda precisa ser implementada/verificada
no caminho nativo e cobrada na política completa.

Mapa ABI atual:0=operator_only,1=operator_size,2=operator_generic,3=operator_state.
O modelo anterior `models.bfc2363c359b.json` tem apenas3modelos e foi verificado
com a versão anterior do C++ (arquivada porSHA); não usar seu índice2 com o ABI
atual. Usar explicitamente nomes/contextos e o modelo de4cabeças acima.

Próximo experimento: mesma previsão de tempo genérica para todos os modelos
com estado, variando a previsão de melhoria. Comparar ordem original, regra
estática por operador, tamanho, estado genérico, estado com cobertura/prêmio,
e adaptação simples. Manter todos os operadores elegíveis até completar a
varredura após a última melhoria; incluir todo overhead e preservar censura.
Uma política pode mudar a trajetória, portanto avaliar em execução completa.

## Integração e primeiro piloto, às 7h01

`operator_state.cpp` reproduziu os dez descritores Python em460estados de23
perfis, erro máximo6.66e-16. O adaptador do legado, que monta índices, custos,
prêmios e incidência de cobertura, foi verificado em20estados de cada uma de
att48/berlin52/brazil58; erro máximo2.22e-16. Esse segundo teste verifica a ligação
com os tipos e conjuntos do código fornecido, além da rotina portátil.

Build separado `build_scheduled/`, gerado com `--initialize-decoder-index
--scheduled`, incorpora os modelos congelados e mantém os15operadores. A
permutação aleatória original é sorteada antes da reordenação; a nova ordem
é aplicada no início da busca local e após melhoria. Não omite operadores.
Controles estático, tamanho e adaptativo evitam os descritores que não usam;
o genérico não calcula cobertura crítica. O modelo completo usa o mesmo
preditor de tempo genérico. A adaptação simples acumula sucessos/tempo com
prior fixo; não reproduz um algoritmo publicado específico.

Piloto `scheduled_development/protocol.16776a180cf7.json`:3geometriasDEV,
uma semente42,6políticas,70iterações. Todos18completaram e mantiveram os custos
2650/5186/20466. A auditoria independente recalculou rotas, custos, monotonicidade,
método e proveniência. Nenhuma geometria reservada foi usada.

| Política | att48(s) | berlin52(s) | brazil58(s) |
|---|---:|---:|---:|
| Ordem original instrumentada |5.516|7.503|3.196|
| Estática por operador |4.232|5.971|2.759|
| Tamanho |4.668|5.424|2.490|
| Estado genérico |4.478|5.327|2.631|
| Estado com cobertura/prêmio |4.101|5.101|2.393|
| Adaptação simples |4.120|5.542|2.481|

São observações preliminares sob número fixo de iterações, não evidência final
de vantagem sob orçamento de tempo. O wrapper mede todo operador e calcula
custos extras em todas as políticas; isso pode inflar a comparação com a ordem
original. É preciso comparar também com o MDM corrigido sem esse wrapper e
remover custos desnecessários das políticas que não adaptam online. Diferenças
pequenas, uma semente e apenas3geometrias não sustentam inferência estatística.

Próximo: coleta de progresso dentro da busca local e prazo comum, com limites
e overruns registrados, mais sementes/geometriasDEV e controle legado limpo.
Não substituir isso por testes de significância sobre18execuções correlacionadas.
