# Recombinação na união de arestas: resultado exploratório

Lote completo de nove geometrias hardTRAIN, auditado independentemente: 79
rotas verificadas. Duas melhorias, ganho médio de 0.1548048% sobre a melhor
rota de entrada. O critério pré-fixado (pelo menos três melhorias e média de
0.2%) FALHOU. Não abrir teste nem apresentar superioridade confirmatória.

Todos os nove modelos restritos retornaram Optimal dentro do orçamento;
esses certificados pertencem apenas à união das arestas de entrada. A falta
de novas arestas limita a vizinhança, não demonstra otimalidade PCCTP.

gil262: 1998→1981, mas 1981 já consta no cabeçalho da instância: NÃO é
novo melhor conhecido. si175: 7006→6968, ainda acima de 6871 no cabeçalho
e de uma rota anteriormente auditada. Não há recorde sustentado neste lote.

O custo de gerar o conjunto de rotas foi excluído: não comparar com algoritmo
completo sob tempo igual. Próxima hipótese possível, ainda não executada:
ampliar explicitamente o grafo de recombinação com arestas ausentes, avaliando
a troca entre variedade e dificuldade do modelo. Requer protocolo próprio
e controle contra reparação inteira existente; não retunar validação/teste.
