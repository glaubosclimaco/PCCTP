# Auditoria da origem do corte de posto de partição

18 de outubro de 2026. Análise teórica, sem abrir validação ou teste.

## Derivação por variáveis de ocupação

Fixe uma partição P0,...,Pk, raiz em P0, e posto conjunto rho >= 1.
Para uma rota inteira, defina zi=1 se ela visita algum vértice de Pi, i>0.
As seguintes relações são válidas no espaço estendido:

    x(delta(P0)) >= 2,
    x(delta(Pi)) >= 2 zi, i=1,...,k,
    sum_i zi >= rho.

A primeira relação usa rho>=1: a rota obrigatoriamente sai do bloco da raiz.
Somar as relações de fronteira com multiplicador 1/2 e a relação de posto
com multiplicador 1 produz

    x(E_entre(P)) >= 1 + rho.

Cada aresta entre partes aparece em exatamente duas fronteiras. Portanto,
o corte candidato é consequência linear de conectividade por ocupação e de
um posto de seleção em blocos. Não é preciso invocar paridade adicional nem
uma comb para essa derivação. A definição de rho continua exigindo o oráculo
conjunto de cobertura e prêmio.

Isso NÃO prova redundância na formulação original em (x,y): zi é indicador
de uma disjunção sobre vários vértices, não a soma de yi no bloco. Tampouco
prova equivalência entre a projeção completa do sistema estendido e o único
corte agregado. Limites zi<=1 podem gerar condições adicionais.

## Conferência racional da testemunha existente

Na testemunha de nove vértices, x(delta(P0))=2 e cada uma das quatro
fronteiras opcionais tem peso 5/4. Consequentemente zi<=5/8 e sum zi<=5/2,
incompatível com rho=3. A soma das cinco fronteiras dividida por dois é
7/2, enquanto o corte exige 4. A deficiência no sistema estendido é 1/2.

O casco completo das seleções em y, verificado no experimento anterior,
não resolve essa deficiência: a projeção em y perde a informação sobre
qual bloco foi ocupado em cada seleção da combinação convexa. Essa é uma
explicação do fortalecimento, não uma demonstração de novidade bibliográfica.

## Fontes primárias consultadas nesta retomada

- Baïou, Barahona e Mahjoub (2000), *Separation of partition inequalities*,
  Mathematics of Operations Research 25(2), 243–254. A introdução define
  desigualdades multicorte com RHS afim no número de partes. O algoritmo
  anunciado separa essa classe com parâmetros fixos; não foi demonstrado
  que se aplica quando o RHS é o posto combinatório rho(P).
  https://www.lamsade.dauphine.fr/~mahjoub/Papers/BBM00.pdf
- Maziero, Usberti e Cavellucci, *Branch-and-cut algorithms for the covering
  salesman problem*, preprint arXiv:2104.01173, seção 3, equações 7–10.
  Apresenta conectividade generalizada e interseção de coberturas. Essa seção
  não apresenta a fórmula com posto conjunto sobre blocos. A ausência nessa
  seção não estabelece ausência na literatura. A cautela sobre as hipóteses
  da CI impressa, registrada em prior_art.md, continua válida.
  https://arxiv.org/html/2104.01173

## Consequência para a decisão científica

Não apresentar a operação de contrair blocos e contar arestas como inovação.
Uma eventual contribuição precisa demonstrar algo além dessa derivação:
caracterização geral de facetas, levantamento não trivial, separação com
garantias para uma subfamília relevante, ou benefício computacional robusto.
As facetas finitas e o ganho de LP já medidos não satisfazem sozinhos esses
critérios. Aceleração inteira e vantagem de ML continuam não demonstradas.
Não mudar limiares, k, seletores nem reabrir as linhas de aprendizado rejeitadas.
