# Piloto diagnóstico de otimização da ordem da rota

Congelado antes do piloto em 18/10/2026. Não é comparação final de algoritmo.

Onze geometrias TRAIN com execução MDM completa de70 iterações, semente42,
em move_profiles_monotonic_elite_bounded. Fixar exatamente o conjunto de
visitas de cada rota. Resolver o TSP não dirigido por grau2 e eliminação
iterativa de subciclos, com HiGHS/SciPy. Arestas binárias, cortes de cada
componente desconexa. Um limite superior igual ao custo MDM é válido;
não há MIP start nessa interface. Preservar a rota MDM se não houver
ciclo inteiro melhor retornado e validado dentro do prazo.

Orçamento10s por geometria incluindo leitura, montagem e separação; tempo
de importação comum registrado separadamente. Não usar melhoria que termine
após prazo; registrar sobretempo e status. Custos e viabilidade são conferidos
no PCCTP original, com o mesmo conjunto de visitas.

Hipótese: a ordem GENIUS/VND ainda deixa melhorias em pelo menos2 das11
rotas. Se passar, implementar comparação por tempo TOTAL com MDM recebendo
o mesmo tempo adicional e sementes pareadas, antes de qualquer alegação de
superioridade. Caso contrário, abandonar este mecanismo isolado sem aumentar
prazo ou trocar instâncias neste lote. Instâncias incompletas são listadas
como limitação, não contabilizadas como sucesso. Teste reservado fechado.

Esta é uma aplicação conhecida de eliminação de subciclos, sem reivindicação
de novidade matemática. Não é a antiga previsão de suporte: aqui o conjunto
é exatamente o incumbent MDM e todas as arestas entre ele podem mudar.

## Resultado

0/11 melhorias, 11/11 ordens provadas ótimas no suporte fixado; a regra
para continuar este mecanismo isolado falhou. Não é prova de ótimo global
do PCCTP. Resultado em route_polishing_results.json. Referência do método
clássico: Dantzig–Fulkerson–Johnson1954, DOI10.1287/opre.2.4.393.
