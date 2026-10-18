# Faceta para a construção de cobertura por pares: prova geral

18 de outubro de 2026. Resultado matemático proposto nesta investigação;
novidade bibliográfica não estabelecida. Não implica vantagem computacional.

## Enunciado e domínio exato

Seja G um grafo completo não orientado com V={r} união B1 união ... união Bq,
q>=4, blocos opcionais disjuntos e não vazios. A única visita obrigatória é r.
Para cada par i<j existe uma demanda cuja vizinhança é Bi união Bj.
Não há outras restrições de seleção: por exemplo, prêmio unitário nos opcionais
e quota 1, automaticamente satisfeita pela cobertura. Seja P o casco convexo
dos vetores de incidência de ciclos simples factíveis, em coordenadas (x,y),
com y_r omitido por estar fixo em 1. Ponha n=|V| e m=n(n-1)/2.

Então P tem dimensão m-1 e

    x(E_entre({r},B1,...,Bq)) >= q

define uma faceta de dimensão m-2, para quaisquer tamanhos positivos dos blocos.

Esta afirmação não se estende automaticamente a quota ativa, outros obrigatórios,
ou instâncias com demandas adicionais. O caso q=3 tem a decomposição em cortes
de cobertura apresentada em partition_facet_generalization.md.

## Redução afim às coordenadas de arestas

Todo ciclo satisfaz x(delta(r))=2 e y_v=x(delta(v))/2 para v!=r.
Logo a projeção em x é bijetiva por uma transformação afim inversível sobre P;
dimensões e facetas relativas são preservadas. Trabalhamos em R^m, onde
a única equação inicialmente necessária ao limite de dimensão é
x(delta(r))=2. Isso fornece dim(P)<=m-1.

Cobrir todos os pares equivale a deixar no máximo um bloco opcional vazio.
Toda rota visita portanto pelo menos q blocos contando a raiz e tem pelo
menos q arestas entre partes. Existem rotas justas: omitir um bloco, escolher
um vértice em cada bloco restante e percorrê-los com r em qualquer ordem.
Inserir outros vértices consecutivamente dentro de seu próprio bloco mantém
exatamente q cruzamentos e preserva factibilidade.

Seja F a face justa. Para determinar seu casco afim, tome qualquer igualdade
a·x=beta satisfeita por todas as rotas de F.

## Passo 1: coeficientes internos nulos e coeficientes externos uniformes

Escolha u,v distintos num mesmo bloco Bi e w em outro bloco ou igual a r.
Pode-se escolher uma rota justa contendo u, não v, e com sucessor w de u:
selecionar q-1 blocos opcionais incluindo Bi e o bloco de w quando aplicável,
omitir outro bloco e ordenar livremente os representantes com a raiz.

Substituir a aresta uw pelo caminho u,v,w mantém a rota em F. Comparando
as igualdades das duas rotas, obtemos

    a_uv + a_vw = a_uw.                         (A)

Analogamente, comece com v em vez de u e insira u antes de w:

    a_uv + a_uw = a_vw.                         (B)

Somar A e B dá 2a_uv=0. Portanto a_uv=0 e a_uw=a_vw.
Como w é qualquer vértice externo, o coeficiente de uma aresta entre dois
blocos depende apenas dos blocos, não das extremidades. Denote-o c_ij,
incluindo o índice 0 para o bloco singleton {r}. Blocos singleton não
exigem a operação de inserção; a conclusão entre blocos permanece válida.

## Lema elementar de custos constantes em ciclos Hamiltonianos

Num grafo completo com s>=4 vértices, se todos os ciclos Hamiltonianos têm
o mesmo custo para pesos c, existem potenciais t_v tais que c_uv=t_u+t_v.

Prova: trocas de duas arestas em ciclos Hamiltonianos dão, para quaisquer
quatro vértices distintos, igualdade entre as três somas de emparelhamentos
c_ij+c_kl, c_ik+c_jl e c_il+c_jk. As trocas existem colocando os quatro
vértices na ordem apropriada num ciclo e invertendo o segmento intermediário.
Fixe âncoras 1,2,3 e defina

    t1=(c12+c13-c23)/2, t2=c12-t1, t3=c13-t1,
    tv=c1v-t1 para v fora das âncoras.

As igualdades dos emparelhamentos fornecem c2v=t2+tv e c3v=t3+tv;
para v,w fora das âncoras, c12+cvw=c1v+c2w fornece cvw=tv+tw.
Assim a representação vale em todas as arestas. Esses potenciais são únicos
em qualquer subgrafo completo com pelo menos três vértices, pela solução
das três equações de um triângulo. Fim da prova do lema.

## Passo 2: compatibilidade das representações locais

Fixe um bloco h a omitir. Escolher um representante de cada outro bloco
opcional e a raiz fornece todos os ciclos Hamiltonianos do grafo quociente
com q vértices, todos com custo beta. Pelo lema, seus coeficientes c_ij
admitem potenciais t_i^(h).

Os grafos quocientes associados a duas omissões distintas se intersectam
num completo com q-1>=3 vértices. Como os coeficientes nas arestas comuns
são os mesmos e potenciais num triângulo são únicos, suas representações
coincidem em toda a interseção. Portanto há potenciais globais t0,...,tq
com c_ij=t_i+t_j. Cada aresta pertence a algum quociente; cada par de
representações contendo um vértice coincide nele.

## Passo 3: comparar omissões

O custo de uma rota justa com omissão h é

    beta = 2(t0 + sum_{i=1}^q ti - th).

Comparar h e k mostra th=tk para todos os blocos opcionais. Escreva esse
valor t. Então a_ij=2t entre blocos opcionais e a_rv=t0+t em arestas da raiz;
os coeficientes internos são zero. Consequentemente

    a = 2t * 1_Eentre + (t0-t) * 1_delta(r),
    beta = 2t*q + 2(t0-t).

Logo toda igualdade em F é combinação da igualdade da partição e da
equação de grau da raiz. Elas são independentes: existem arestas entre
blocos opcionais, onde a primeira tem coeficiente 1 e a segunda 0.
Portanto dim(F)=m-2.

Por fim, uma rota com um representante de TODOS os q blocos opcionais
tem q+1 cruzamentos; é factível e está fora de F. Assim dim(P)>=m-1.
Junto ao limite superior, dim(P)=m-1, concluindo a prova de faceta.

## Auditoria de escopo e interpretação

- A prova usa q>=4 exatamente para aplicar o lema e colar representações
  locais em interseções com triângulos. A falha em q=3 é consistente com
  a decomposição e a codimensão 3 já certificadas.
- A quota não é ativa. Não reivindicar que este teorema demonstra interação
  essencial de prêmio e cobertura.
- Com todos os blocos singleton, o corte é um limite de cardinalidade de
  visitas. A passagem para blocos maiores tem natureza de lifting por
  duplicação; é especialmente importante confrontá-la com resultados antigos
  de lifting em GTSP e de cardinalidade antes de afirmar novidade.
- A prova caracteriza uma família infinita de instâncias, não certifica que
  o corte seja faceta em uma instância arbitrária do benchmark PCCTP.
- O resultado precisa de revisão independente da prova e da literatura;
  não é, isoladamente, critério de encerramento do artigo.
