# Fixação certificada por custos reduzidos

18/10/2026, antes de medir. Mesmas11TRAIN completas e separador GSEC congelado.
Não reivindicar método inédito. Certificar as fixações por um dual racional:
arredondar multiplicadores de linhas à escala10^6, corrigir sinais segundo
limites finitos, recomputar custos reduzidos e bound usando só inteiros.
Dados desta formulação devem ter coeficientes/limites finitos inteiros;
rejeitar se isso não ocorrer ou se houver risco de overflow.

Para min c*x e limites l<=x<=u, multiplicadores pi dão
B=pi*b+sum_j min(r_j*l_j,r_j*u_j), r=c-A^T*pi.
Escolher b_i=limite inferior da linha se pi_i>=0, superior se negativo.
O sinal de pi deve ser compatível com limite finito. O bound vale exatamente
para qualquer pi assim escolhido. Em binária livre, forçar valor oposto ao
que minimiza r_j*x_j aumenta o bound em abs(r_j). Se B+abs(r_j)>U,
fixar no minimizador preserva todas soluções com custo<=incumbent U e,
portanto, pelo menos um ótimo. Não chamar fixação válida para TODAS soluções.

Verificação antes do piloto: exemplo enumerado; cada solução de custo<=U
deve respeitar toda fixação. Confirmar vetor SCF MDM nas fixações emitidas.
Gate de tamanho, antes dos resultados: fixar>=30% das binárias inicialmente
livres em>=6/11 geometrias, com preparação completa<=2s. Se passar, comparar
com/sem fixação sob mesmo10sTOTAL, incluindo preparação e GSEC em ambos,
para isolar efeito de fixação. Se falhar, não aumentar escala/prazo/rodadas
para salvar o mecanismo. Teste reservado fechado.

## Resultado

Gate falhou:4/11 contra6 exigidas.4097fixações com certificados exatos,
auditadas por implementação independente em Pythoninteiro,sem resolverLP.
Não executar piloto inteiro deste mecanismo nem ajustar o limiar.
