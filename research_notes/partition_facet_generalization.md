# Generalização de facetas: diagnóstico construtivo

18 de outubro de 2026. Não envolve geometrias de validação/teste.

Considere raiz isolada no bloco P0, q blocos opcionais não vazios, e uma
demanda com vizinhança Pi união Pj para cada par de blocos opcionais.
Prêmios unitários, quota 1. Cobertura equivale a ocupar pelo menos q-1
blocos opcionais, portanto a desigualdade candidata é x(Eentre)>=q.

## Obstrução demonstrada para q=3

Cada união Pi união Pj é forçada pela demanda correspondente. Como a raiz
está fora, cada corte de fronteira dessa união tem peso pelo menos 2.
Somar os três cortes conta exatamente duas vezes cada aresta entre partes:
uma aresta raiz–Pi pertence a duas fronteiras; uma aresta Pi–Pj também
pertence a duas. Logo sua soma, dividida por dois, é x(Eentre)>=3.

Portanto, para quaisquer tamanhos dos três blocos, o corte é consequência
linear dos cortes de cobertura já conhecidos. Na face em igualdade, os três
cortes de cobertura são simultaneamente justos. Não se pode extrapolar uma
conjectura de faceta dos casos q>=4 para q=3.

Para tamanhos (2,2,2), o certificado exato dá dimensão 20 ao politopo e 17
à face: codimensão 3, logo não é faceta. As três equações de fronteira,
somadas às de grau, fornecem o limite superior 17, atingido pelos geradores.

## Certificados finitos adicionais

`experiments/pcctp_cuts/probe_partition_facets.py` gera rotas factíveis com
semente fixa. Nas rotas justas, omite exatamente um bloco opcional e percorre
os vértices escolhidos de cada bloco consecutivamente. Isso garante exatamente
q cruzamentos. Rotas gerais permitem todos os blocos e ordem arbitrária.

Independência afim é certificada por eliminação no primo 1000003: menor
inteiro não nulo módulo o primo é não nulo sobre os racionais. Limites
superiores vêm de equações de grau e da igualdade do corte. Assim, amostragem
é apenas uma forma de encontrar os geradores, não uma inferência estatística
de dimensão. Falha em alcançar o limite seria inconclusiva.

| Tamanhos dos blocos opcionais | Dimensão do politopo | Dimensão da face | Conclusão finita |
|---|---:|---:|---|
| 2,2,2 | 20 | 17 | Não faceta |
| 2,2,2,2 | 35 | 34 | Faceta |
| 2,2,2,2,2 | 54 | 53 | Faceta |
| 2,2,2,2,2,2 | 77 | 76 | Faceta |
| 3,2,2,2 | 44 | 43 | Faceta |
| 1,2,2,2 | 27 | 26 | Faceta |

Os resultados e todas as rotas geradoras estão em `partition_facet_probe.json`.
Execução reproduzível: `/usr/bin/python3 experiments/pcctp_cuts/probe_partition_facets.py`.

## Hipótese e próximo passo

Os resultados sugerem investigar q>=4, inclusive blocos desiguais e singleton.
Não constituem prova para todo q ou todo tamanho. O próximo passo é uma
prova por coeficientes de hiperplano: qualquer igualdade válida em todas as
rotas justas deve ser combinação dos graus e do corte, usando substituição
de vértices, inserção dentro de bloco e permutação de pelo menos três blocos
opcionais visitados. Essas operações precisam ser formalizadas antes de
afirmar teorema. O caso misto quota–cobertura exige condições próprias;
esta construção usa quota não ativa e não prova sinergia de prêmio.

Uma faceta genuína também pode ser uma faceta conhecida. A conferência
bibliográfica continua necessária, especialmente para partições de group
Steiner e lifting em GTSP. Não confundir estes certificados com evidência de
aceleração do solver ou de contribuição suficiente para submissão.

## Atualização posterior

Uma prova geral para q>=4 está escrita em `partition_facet_theorem.md`,
ainda sujeita a revisão da prova e da novidade. Ela substitui a hipótese
aberta desta nota; os certificados finitos permanecem como conferências.
