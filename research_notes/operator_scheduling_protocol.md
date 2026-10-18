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
