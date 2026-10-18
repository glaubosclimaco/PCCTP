# Hipótese de seleção conjunta de opcionais e âncoras

Protocolo de desenvolvimento, 18/10/2026. Não constitui configuração final
congelada nem evidência de novidade. Teste e validação reservados intactos.

O piloto de seleção opcional com reparação anterior não mostrou vantagem de
ML. A DP em fronteira completa torna possível testar outra decisão: quais
obrigatórios liberar para mudar o roteamento, além dos opcionais cuja seleção
pode mudar. Quatro obrigatórios próximos dos opcionais não são necessariamente
uma boa escolha: no piloto inicial berlin52, ficaram em5244, enquanto quatro
aleatórios chegaram a5186. São resultados de uma semente, para formular hipótese,
sem inferência estatística nem vantagem sobre MDM.

`dp_trajectory.py` compara quatro controles com a mesma DP, limite18componentes,
0.2s por consulta e5s adicionais totais:12opcionais;8opcionais+4obrigatórios
próximos;8+4aleatórios;8+4próximos com adaptação dos operadores opcionais.
Toda construção e seleção online entram no relógio. A biblioteca nativa é
compilada antes do orçamento, mas a construção do modelo fica dentro dele.
O piloto usa uma variante por geometria de desenvolvimento. Checkpoints MDM
são históricos e esparsos; não sustentam alegação final de aceleração.

Nova coleta `dp_labels.py`: apenas geometrias de treino, uma variante por
geometria, estado da primeira iteração. Para12conjuntos opcionais, coleta as
duas escolhas de âncoras (próximas/aleatórias), em mesma vizinhança clique e DP.
Alvos: ganho normalizado observado e tempo do reparo; consultas limitadas por
tamanho permanecem, com ganho zero e motivo explícito. Certificados distinguem
ótimos locais de reparos interrompidos. Hashes congelam fontes e baselines.

Descritores: oito opcionais básicos; número de componentes; número de âncoras;
proxy de economia de remoção das âncoras (somente proxy: continuam obrigatórias);
distância das âncoras aos opcionais; dez assinaturas condicionais de cobertura
e prêmio. Medir custo de extração. Não transferir silenciosamente os rótulos
antigos de MIP/route_context nem seus modelos.

Antes de alegar valor de ML, comparar com regra de economia, escolha aleatória,
adaptação simples e regra que conhece o mesmo número de componentes/limite.
Avaliar ranking por geometria independente, depois trajetórias com orçamento
total incluindo inferência. Se não houver sinal adicional, preservar a refutação
e investigar outro papel para ML, sem ajustar continuamente até obter significância.

As escolhas próximas/aleatórias não são intervenções idênticas em footprint;
o número de componentes deve entrar tanto no modelo quanto nos controles.
Novidade de contração, DP e destroy/repair exige confronto com antecedentes
HDR/DRHG e outros trabalhos; não apresentar mecanismos clássicos como novos.

## Resultado da primeira coleta

480consultas em20geometriasdetreino,24melhorias e123limites de tamanho. Auditoria
passou. LOGO exploratório: ganho médio0.0003639 dosquatro modelos versus0.0004085
da regra de menor número de componentes. Nenhum ganho adicional de ML demonstrado.
O piloto completo teve24trajetórias em6geometrias e nenhuma superou os checkpoints
MDM comparáveis; todos os resultados negativos estão preservados.
