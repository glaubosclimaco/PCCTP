# Plano de execução e critérios de evidência

Atualizado em 18 de outubro de 2026. Branch: `codex/pcctp-learning-submission`.

O usuário autorizou continuidade autônoma até uma contribuição científica e um
artigo completo, revisado e pronto para submissão. Não aguardar aprovação para
etapas dentro desse escopo. ML já integra a linha de vizinhanças; se ela falhar,
investigar alternativas com ML sem apagar resultados negativos.

1. Recuperar e verificar o MDM original, com instrumentação observacional isolada.
   Manter a reconstrução anterior identificada como diagnóstico. Registrar casos
   que excedam tempo e corrigir o protocolo de coleta sem omiti-los.
2. Verificar reparação inteira e descritores condicionais por enumeração
   independente em instâncias pequenas. Medir custo de construir os descritores.
3. Investigar em desenvolvimento/treino o valor adicional dos descritores
   conjuntos de cobertura e prêmio: regras simples, ranking genérico, mesmos
   modelos com descritores brutos e conjuntos. Agregar por geometria.
4. Implementar a política iterativa completa. Comparar com MDM original e busca
   adaptativa sem aprendizado, sob orçamento total de tempo que inclua seleção,
   inferência e reparação. Usar sementes pareadas e estados iniciais idênticos.
5. Congelar configurações e hipóteses antes da avaliação reservada. Separar
   geometrias inteiras. Validação serve à seleção limitada; teste não orienta
   ajustes. Se uma hipótese falhar, documentar e definir outra investigação;
   não reutilizar o teste como confirmação independente.
6. Avaliar tamanho de efeito, incerteza por geometria, múltiplas comparações,
   curvas de tempo/qualidade, falhas e censura. Não tratar vizinhanças, variantes
   correlacionadas ou snapshots como instâncias independentes.
7. Confrontar a contribuição com fontes primárias da literatura atual. Redigir
   em `submission/` apenas afirmações sustentadas. Alvo de 20–30 páginas com
   formulação, análises, diagramas e gráficos necessários.
8. Auditar dados → tabelas/figuras → texto; contas e formulação; existência e
   metadados das referências; compilar e inspecionar o PDF completo. Só então
   classificar o manuscrito como pronto para submissão. Não prometer aceitação.

## Limitações já encontradas

- O benefício dos descritores brutos não apareceu no piloto anterior.
- Melhorias em uma reparação de um estado inicial não demonstram superioridade
  de um algoritmo completo nem valor adicional frente ao MDM após 70 iterações.
- Algumas execuções do MDM excedem 180 segundos antes da exportação final.
  Essas geometrias não podem ser silenciosamente excluídas da conclusão.
- As facetas verificadas em duas instâncias pequenas não constituem um teorema
  geral, e a novidade das desigualdades ainda não foi estabelecida.

CTP permanece fora da implementação atual por decisão explícita do usuário.
