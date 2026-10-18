# Comparação inteira após diagnóstico de relaxação

18/10/2026, fixado antes dos custos inteiros com cortes. O gate da raiz
passou9/11; cortes clássicos reduziram>=10% da folga, sem retuning.

Mesmas11TRAIN,10s totais por braço,highspy1.13.1,SCF,seed42,1thread,
startMDM em ambos. Controle sem cortes, tratamento com exatamente o separador
pré-fixado (máximo3rodadas/2s). Preparação,cortes,reconstrução e captura de
incumbentes debitadas do orçamento10s. Alternar ordem dos braços por geometria.
Não usar bounds da vizinhança como globais; aqui os modelos são completos.

Critério para estudo maior: pelo menos2melhorias de custo MDM no tratamento
e mais melhorias que controle; OU pelo menos2geometrias adicionais com prova
de ótimo no tratamento, sem perda de qualidade primal. Se falhar, não aumentar
rodadas/prazo nem selecionar subconjunto favorável. Tempos/provas isolados são
diagnósticos, não estatística confirmatória. Sem testes reservados.

## Resultado auditado

11pares:plain3provas e1melhoria;cuts4provas e0melhorias. Gate não passou.
GanhoLP não se converteu em vantagem suficiente no inteiro. AuditoriaPASS
emgsec_integer_audit.json. Não retunar rodadas/prazo neste lote.
