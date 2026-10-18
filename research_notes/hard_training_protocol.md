# Mapeamento dos casos TRAIN incompletos

18/10/2026. Incluir TODOS os9perfis de treino incompletos,não selecionar por
resultado. Esses perfis têm instrumentação observacional de movimentos; seus
tempos NÃO são baseline limpo nem representam70iterações concluídas.

Extrair BS/AS/AT da coluna MDM da tabela explicitamente rotulada
tab:mdm-grasp em artigos/mdm/main_review.tex. Guardar linha e não comparar
tempos publicados com a máquina atual. Valor publicado é referência local
do artigo, não prova de melhor conhecido em toda a literatura atual.

Executar baseline original com as duas correções isoladas,sem profiling,
budget10s,seed42,coorteinteira9casos. Usar wrapper temporal já verificado:
inicialização cobrada,rotas validadas,melhorias tardias excluídas,overrun e
timeouts registrados. Não confundir ausência de incumbente em10s com inviabilidade.
Objetivo desta fase é diagnosticar margem,construção e custo de busca; não
selecionar limiar estatístico ou declarar contribuição. Teste reservado fechado.
