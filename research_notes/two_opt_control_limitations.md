# Limitação de implementação no controle 2-opt atual

19/10/2026, lote com fallback ainda em execução. O controle support_2opt usa
Python e recalcula o custo completo de cada reversão, enquanto LKH é compilado.
Portanto, comparar os dois não isola apenas a escolha de vizinhança: velocidade
do kernel e custos de interface também mudam. A comparação com MDM compilado
continua sendo o controle principal de desempenho; não apresentar vantagem
contra esse 2-opt como evidência de superioridade intrínseca do LKH.

Antes de uma avaliação final, preparar um controle2opt compilado, com custo
delta O(1) em matriz simétrica e orçamento de tempo cobrado. Manter construção
do pool e regras de aceitação iguais; validar independentemente custo e ausência
de movimento2opt melhorante quando terminar sem timeout. A preparação desse
controle não altera o lote atual, que deve ser concluído e preservado como piloto.
Não reabrir validação/teste para ajustar esse componente.
