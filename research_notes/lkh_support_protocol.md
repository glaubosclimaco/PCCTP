# Piloto diagnóstico de troca de suporte e ordenação

18/10/2026, antes das execuções. Somente9hardTRAIN,seed42, mesmos incumbentes
MDM10s já salvos. Cinco segundos ADICIONAIS por braço, debitando geração,
validação e ordenação; não é comparação de algoritmo completo com MDM15s.
Três braços: suporte fixo+LKH; suporte variável+2opt; suporte variável+LKH.
Vizinhança: remover opcional, adicionar opcional, trocar opcional visitado por
não visitado, mantendo factibilidade original. Ordenar até8candidatos de menor
custo bruto por rodada, LKH até0.2s/candidato (mesmas opções exceto prazo),
controle2opt primeira melhora; aceitar só melhoria estrita antesprazo.
Sem perturbação aleatória ou sintonia usando validação. Se não há melhoria,
encerrar rodada. Começar todosbraços com mesmo incumbente, incluindo ordenação
inicial em cada um. Este piloto testa complementaridade de suporte/ordenação,
não novidade de ILS, não recordes nem eficácia endtoend. Gate diagnóstico:
suporte+LKH supera suporte fixo+LKH em>=3/9 e nenhuma inviabilidade; relatar
controle2opt independentemente. Não libera teste reservado em hipótese alguma.
