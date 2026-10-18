# Extensão do argumento a postos intermediários

18 de outubro de 2026. Prova proposta; novidade e revisão independente pendentes.

## Enunciado

Considere grafo completo, uma raiz obrigatória singleton e q blocos opcionais
não vazios de tamanhos arbitrários. Para inteiros 3<=r<=q-1, imponha que
a seleção ocupe pelo menos r blocos. Isso pode ser realizado somente por
cobertura no PCCTP: criar uma demanda para cada união de q-r+1 blocos,
com prêmio unitário e quota 1 não ativa. Seja m o número de arestas.

O politopo de rotas tem dimensão m-1 e a desigualdade

    x(Eentre)>=r+1

define uma faceta de dimensão m-2. O teorema de cobertura por pares é
o caso r=q-1. Aqui r denota posto, não a raiz, denotada por 0 na prova.

## Equivalência da cobertura

Uma seleção com menos de r blocos deixa pelo menos q-r+1 blocos vazios,
logo deixa alguma demanda inteira sem cobertura. Uma seleção com pelo
menos r blocos deixa no máximo q-r blocos vazios, portanto intercepta
toda união de q-r+1 blocos. Não há restrição adicional sobre quantos
vértices de cada bloco são selecionados.

## Prova da propriedade de faceta

Elimine y pelos graus como em partition_facet_theorem.md. A raiz dá o
limite superior de dimensão m-1. Rotas com exatamente r blocos opcionais,
cada um percorrido num único segmento, são factíveis e justas. Fixe uma
igualdade a·x=beta em todas essas rotas.

As duas inserções opostas de vértices u,v no mesmo bloco e antes de um
sucessor externo w dão a_uv+a_vw=a_uw e a_uv+a_uw=a_vw. As rotas de base
existem porque podemos escolher r blocos contendo os blocos prescritos e
ordenar representantes livremente. Portanto coeficientes internos são zero
e os externos dependem somente dos blocos. Escreva-os c_ij.

Para cada subconjunto S de r blocos opcionais, o completo em {0} união S
tem r+1>=4 vértices. Todas suas ordens Hamiltonianas são rotas justas,
logo pelo lema já provado c_ij=t_i^S+t_j^S nesse completo.

Considere dois conjuntos S,S' que diferem por troca de apenas um bloco.
Seus completos se intersectam em {0} união (S interseção S'), com r>=3
vértices. A unicidade dos potenciais num triângulo força igualdade de
suas representações na interseção.

Para qualquer bloco i, o grafo cujos vértices são os subconjuntos S de
tamanho r contendo i, adjacentes por uma troca que preserva i, é conexo:
troque sucessivamente cada elemento de S\\S' por um elemento de S'\\S.
Assim t_i^S é independente de S. O mesmo vale para t_0, pois o grafo
de todos os subconjuntos de tamanho r é conexo. Isso define potenciais
globais t_0,...,t_q. Toda aresta de blocos aparece em algum desses
completos, logo c_ij=t_i+t_j globalmente.

O custo de qualquer rota justa com suporte S é

    beta=2t_0+2 sum_{i em S} t_i.

Dados i!=j, escolha r-1 outros blocos: isso é possível por r<=q-1.
Comparar os suportes obtidos ao acrescentar i ou j mostra t_i=t_j=t.
Portanto

    a=2t*1_Eentre+(t_0-t)*1_delta(0),
    beta=2t(r+1)+2(t_0-t).

As únicas igualdades nessas rotas são, assim, as duas indicadas, independentes.
Sua dimensão afim é m-2. Por validade do corte, a face inteira justa
tem a mesma dimensão: contém essas rotas e satisfaz as mesmas duas equações.
Uma rota com r+1 blocos opcionais é factível, cruza r+2 vezes e está fora
da face. A dimensão do politopo é portanto m-1, provando a afirmação.

## Por que excluir posto 2

Para r=2 a cobertura exige todas as uniões de q-1 blocos. Denote por S_i
a união de todos os blocos opcionais exceto Bi. Cada S_i é forçado,
e a raiz está fora: x(delta(S_i))>=2. Somando em i, uma aresta da raiz
é contada q-1 vezes, e uma aresta entre blocos opcionais é contada duas.
Usar x(delta(0))=2 fornece

    sum_i x(delta(S_i))=2*x(Eentre)+2(q-3)>=2q,

isto é, x(Eentre)>=3. Logo o corte de posto 2 é combinação dos cortes
de cobertura e do grau da raiz para qualquer q>=3. Não alegar que
isso, sozinho, prova não facetness em toda degeneração possível.

## Conferências e limites

Certificados finitos adicionais em general_rank_facet_checks.json: q=5,r=3;
q=6,r=3; q=6,r=4; e q=5,r=3 com tamanhos(1,2,1,3,1). Todos atingem
dimensões m-1 e m-2 por menores modulares. Não substituem a prova acima.

A fonte Stephan2009 confirma que cardinalidade é antecedente conhecido.
O argumento aqui ainda pode ser caso de lifting geral por blocos; sua
novidade permanece por verificar. A realização exige potencialmente muitas
demandas e não garante uma configuração Euclidiana plana. Nada aqui demonstra
ganho no benchmark, separação eficiente para partições arbitrárias ou aceleração.
