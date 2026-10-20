# Imitação do conjunto ótimo de visitas

Regra escrita antes de comparar custos nas geometrias de treino. Não é método
novo. É regressão logística, aplicação já usada em otimização combinatória:
prever o suporte de uma solução ótima e resolver o que resta.

O aprendizado anterior escolhia operador, padrão minerado ou aresta a reabrir.
Aqui a decisão é outra: quais vértices opcionais entram no ciclo.

## Rótulos

Somente as oito geometrias de desenvolvimento com ótimo provado. A rota do
elite de semente 42 entra como rótulo só se o custo dela for exatamente o ótimo
provado. Opcional visitado nessa rota é classe 1. Opcional de fora é classe 0.
pr76 fica de fora porque o ótimo não foi provado. Validação e teste permanecem
fechados.

Ótimos usados apenas para aceitar o rótulo: att48 2650, berlin52 5186,
brazil58 20466, eil51 273, eil76 332, gr48 3571, hk48 7593, st70 419.

## Modelo

Cinco descritores, sem seleção posterior:

1. prêmio dividido pelo maior prêmio dos visitáveis
2. fração das demandas cobertas pelo vértice
3. fração das demandas em que ele é o único cobertor mais próximo
4. distância ao obrigatório mais próximo, dividida pela mediana dos custos positivos
5. distância média aos obrigatórios, na mesma escala

Regressão logística com intercepto, penalidade L2 igual a 1 nos coeficientes
exceto o intercepto, pesos inversos à frequência da classe. Padronização só com
as linhas de treino. Nenhuma busca de hiperparâmetro.

## Decisão e controle

Obrigatórios entram sempre. Os opcionais são ordenados pela probabilidade.
Entra o menor prefixo que cumpre quota e cobertura. Depois entra todo opcional
restante com probabilidade pelo menos 0,5. Conjunto inviável conta como pior
que o MDM.

O controle sem aprendizado ordena por cobertura mais prêmio menos distância ao
obrigatório mais próximo e fica com o mesmo número de opcionais que o conjunto
aprendido. Se esse conjunto empatado em tamanho for inviável, conta como pior.
Não se acrescentam vértices depois para torná-lo viável.

O ciclo do conjunto é o modelo de fluxo resolvido em no máximo 15 segundos.
Se o solver devolver solução viável, o custo dela entra na comparação mesmo
sem prova. Falha ou ausência de solução conta como pior.

## Onde se mede

Somente as primeiras variantes de treino cuja execução perfilada completou 70
iterações: rd100 VT20_T20_W60_75, pr144 VT28_T29_W87_50, u159 VT31_T32_W96_50,
kroB150 VT30_T30_W90_25, gr137 VT26_T28_W83_25, kroC100 VT20_T20_W60_75,
pr136 VT26_T28_W82_25, bier127 VT24_T26_W77_50, lin105 VT21_T21_W63_25,
gr96 VT18_T20_W58_75, kroA150 VT30_T30_W90_50. Semente 42. O custo de referência
é o custo final desse arquivo. Instância que não tiver 70 iterações completas
sai da conta e não é substituída.

O resultado melhora só se o conjunto aprendido ficar estritamente abaixo do
MDM em pelo menos três dessas geometrias e, além disso, ficar estritamente
melhor que o controle de mesma cardinalidade em mais geometrias do que o
contrário. Caso contrário, não há contribuição de aplicação. Não se altera
limiar, descritores nem a lista depois de ver os custos.

## Resultado

Arquivo `experiments/ml_lns/visit_set_prediction.json`, código
`c66c87c0de17b1f85c8fbd2d391285bfeadea5389b53254876ffdad8212633fe`.
Rótulos: 63 opcionais positivos em 120. As 11 geometrias foram resolvidas.
A métrica é o custo da solução viável, e o problema é de minimização. O tempo
não entra na comparação. O limite de 15 segundos só corta o modelo inteiro.
A referência do MDM é o custo ao fim das 70 iterações, não um prazo pareado.

O conjunto aprendido ficou estritamente abaixo desse custo em zero geometrias.
Ficou igual em três (pr136, bier127, gr96) e estritamente pior em oito. Contra
o controle de mesma cardinalidade: 1 melhor, 3 piores, e o restante empate ou
controle inviável. A regra não foi satisfeita.

Dos oito piores, só rd100 e u159 têm o ciclo do conjunto previsto provado.
Neles o ótimo do conjunto fixado é pior que o MDM: 5621 contra 5421, e 29147
contra 28613. Nos outros seis, o número é o incumbente aos 15 segundos, e em
vários o limiar 0,5 tinha mantido todos os opcionais. Isso não é warm start
do modelo completo.

O limiar 0,5 manteve quase todos os opcionais em várias geometrias, e o
modelo de 15 segundos não alcançou o custo das 70 iterações do MDM. Onde o
conjunto de fato encolheu, o custo piorou: rd100 5621 contra 5421 e u159
29147 contra 28613. Em bier127 o aprendizado empatou com o MDM em 74745 e
o controle ficou em 74944. Isso não é melhora sobre o MDM.

Não há contribuição de aplicação. Limiar, descritores e lista não serão
alterados neste lote.
