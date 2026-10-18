# Realização mista da família de facetas

18 de outubro de 2026. Corolário da prova em partition_facet_theorem.md.
Não é uma nova família de politopos nem prova de novidade bibliográfica.

## Construção

Para q>=4, use raiz obrigatória r, q-2 blocos opcionais singleton e dois
blocos opcionais A e B com dois vértices cada. Todas as uniões de pares
de blocos são vizinhanças de cobertura, EXCETO A união B. Prêmio unitário
em cada vértice opcional, prêmio zero na raiz, quota Q=q-1.

O conjunto de seleções viáveis é exatamente o mesmo da construção com
TODAS as demandas por pares e quota não ativa. Portanto o politopo de rotas
é idêntico e a desigualdade x(Eentre)>=q define faceta por aquele teorema.

## Prova da equivalência de seleções

As demandas presentes formam o grafo de obrigações Kq menos a aresta AB.
Um conjunto de blocos vazios deve ser independente nesse grafo. Logo ele
tem tamanho no máximo dois; o único par vazio possível é {A,B}.

Se A e B estão vazios, todos os q-2 blocos singleton precisam estar ocupados,
pois há demanda para cada par envolvendo um deles. Essa seleção tem
exatamente q-2 vértices e prêmio q-2, violando Q=q-1. Nenhuma seleção com
menos de q-1 blocos é portanto conjuntamente viável.

Reciprocamente, ocupar pelo menos q-1 blocos garante todas as demandas
por pares, inclusive a omitida, e seleciona pelo menos q-1 vértices de
prêmio unitário. Assim a quota é satisfeita. Isso prova a equivalência.

## As duas fontes de obrigação são essenciais ao posto conjunto

Somente cobertura permite selecionar os q-2 singleton, logo seu posto é
q-2. O posto conjunto é q-1.

Somente prêmio também tem posto menor que q-1. Para q=4, ambos os blocos
duplos dão prêmio 4 em dois blocos, enquanto Q=3; um bloco tem prêmio
máximo 2, então o posto de prêmio é exatamente 2. Para q>=5, os dois
blocos duplos e q-5 singleton dão prêmio q-1 em q-3 blocos. São os blocos
de maior capacidade: q-4 blocos têm capacidade no máximo q-2. Portanto
o posto de prêmio é q-3 para q>=5.

Consequentemente o posto conjunto supera estritamente tanto o de cobertura
quanto o de prêmio em toda esta família infinita. Para q=4: 2,2,3; para
q=5: 3,2,4; para q>=5 em geral: q-2,q-3,q-1.

Isso não prova que o corte não seja dedutível de outras famílias conhecidas
em (x,y), nem que seja uma nova faceta: trata-se de uma realização mista
do MESMO politopo já descrito na construção de cobertura pura. O ganho
sobre dois números de posto separados não é ganho sobre os fechos completos
das formulações separadas. Preservar essa distinção no eventual artigo.

## Conferência exaustiva e revisão da prova-base

O verificador `verify_mixed_facet_realization.py` enumerou 1984 seleções
no total para q=4,...,8; os conjuntos viáveis coincidem exatamente e todos
os postos acima foram conferidos. Resultados em `mixed_facet_realization.json`.

Revisão interna da prova-base: conferidos a inversibilidade afim da eliminação
de y, os dois sentidos da inserção intrabloco, a existência das rotas com
sucessor externo prescrito, o lema para q=4, a unicidade de potenciais em
triângulos e a rota fora da face para a dimensão total. Não foi identificado
erro nessa revisão. Ela não substitui parecer independente.

## Comparação com lifting já publicado

Sundar–Rathinam, preprint arXiv:1508.01813, proposição 3.1, descreve lifting
sequencial ao liberar variáveis de visitas fixadas, atribuído a Padberg1975,
e sua relação com Fischetti et al.1995. Essa operação mantém o conjunto de
vértices, enquanto a prova aqui usa inserção de vértices num bloco e colagem
de quocientes. Essa diferença de apresentação NÃO demonstra que o resultado
seja novo ou não dedutível de lifting geral. A fonte trabalha com clusters
disjuntos obrigatórios, hipótese distinta das obrigações por pares usadas aqui.

Fonte primária: https://arxiv.org/pdf/1508.01813 (seção 3, proposição 3.1).
Ainda falta confronto completo com cardinality-constrained cycle polytopes e
lifting por duplicação. Não reivindicar contribuição publicável nesta etapa.
