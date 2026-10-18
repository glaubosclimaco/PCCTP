# Diagnóstico de seleção e roteamento por cortes de subciclo

18/10/2026, antes de executar. O piloto de ordem fixa mostrou11/11 suportes
MDM com ordem ótima;0melhorias. Não aumentar o prazo daquele mecanismo.

Mesmas11geometrias TRAIN completas, semente42,10s totais por caso. Agora
permitir TODAS as visitas e arestas, preservando obrigatórios, quota e cobertura.
Formulação binária x,y com graus, linking e cortes condicionais
x(delta(S))>=2*y_v para cada v em componente S sem a raiz.
Nunca impor x(delta(S))>=2 incondicionalmente a um componente opcional.
Separação externa de soluções inteiras retornadas por HiGHS. Cutoff igual
ao custo MDM. Sem MIP start. Guardar incumbent MDM se não houver ciclo
melhor validado antes do prazo. Não confundir ausência de melhoria com prova.

Critério para investir numa comparação completa: pelo menos2/11 melhorias
estritas dentro do prazo. Resultado de pós-processamento não demonstra
vantagem com tempo total justo; se passar, comparar MDM com tempo extra
equivalente, solver direto e combinação, com sementes pareadas. Sem mudar
limiar, instâncias ou prazo após resultados. Teste reservado fechado.
