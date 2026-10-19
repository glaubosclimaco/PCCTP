# Busca de suporte em incumbentes fortes — protocolo antes da coleta

19/10/2026. Nove hardTRAIN, mesma rota de entrada native+LKH do lote
fixed_work (sete70iterações, dois censurados, explicitamente identificados).
Três braços: LKH somente ordem, mudanças de suporte com2opt C++ verificado,
mesmas mudanças comLKH. Todos recebem30s adicionais; o limite é teto, não
obrigação de consumir tempo após convergência. Mesma enumeração de remoção,
inserção e troca1por1 já implementada, oito candidatos pré-ordenados.
Sem alegar nova vizinhança. Rotacionar ordem entre instâncias.

Este teste distingue uma limitação dos estados iniciais muito curtos de uma
limitação da busca; não é avaliação end-to-end nem repete validação. Gate:
supportLKH melhora>=3/9partidas e média>=0.2%; supera controleC++>=5/9.
Registrar clocks/CPU, preservar resultado bruto e rejeitar melhorias retornadas
após30sCLOCK_MONOTONIC. Auditar rotas e resultados por matriz crua. Caso
negativo, não aumentar novamente orçamento da mesma busca sem outra hipótese.

Literatura consultada: Reactive GRASP, DOI10.1051/ro/2021057, página Numdam
confirmaGRASPreativo+VND. Vidaletal arXiv1401.3794v2, resumo confirma busca
com seleção implícita via caminho mínimo com recursos em representação
exaustiva. Não foi implementada essa seleção aqui nem alegada equivalência
aPCCTP: cobertura pode exigir estados adicionais. Leitura atual dos resumos.
