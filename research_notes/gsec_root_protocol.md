# Diagnóstico da relaxação com GSEC clássicos

Congelado antes da medição,18/10/2026. Mesmas11TRAIN completas. Relaxação SCF
sem integridade; acrescentar x(delta(S))>=2y_v, raiz fora de S, via mincut
raiz–v em capacidades x. Converter capacidades por piso a inteiros de escala
10^6 somente para gerar candidatos; conferir violação na solução original
em ponto flutuante antes de acrescentar. Isso não altera validade dos cortes.

No máximo3rodadas,2s totais por geometria incluindo leitura/modelo/separação.
Não usar objetivo primal de LP interrompido como bound certificado; registrar
apenas ótimo de LP. Não abrir teste reservado. Gate para piloto inteiro:
reduzir pelo menos10% da folga MDM–LP inicial em pelo menos6/11geometrias,
com LP final ótimo dentro do prazo. Se passar, comparar algoritmo com/sem
separação usando10s totais,incluindo cortes,mesmo HiGHS e startMDM. Isso não
é nova desigualdade nem garante aceleração. Não ajustar rodada/escala/prazo
após resultados. Verificar separador por enumeração de cortes em exemplo
pequeno e cortes emitidos contra rota MDM e startSCF completo.
