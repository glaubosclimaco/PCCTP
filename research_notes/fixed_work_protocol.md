# Diagnóstico de convergência com trabalho fixo do MDM

19/10/2026, antes das execuções. Pilotos15s têm censura variável e ganhosfrágeis;
comparação com valores publicados usualmente obtidos em orçamentosmaiores
exige conhecer a convergência da baseline. SemnovoML ou tuningreservado.
Todos9hardTRAIN,seed42,MDM70iterações com correçõeslegadas jádocumentadas,
semprofiling e semdeadlinealgorítmico(PCCTP_BUDGET0). Mesmo binário temporizado
combudgetinativo mantém70iterações originais. Limite externo1800s por caso
é apenas salvaguarda; se expirar,preservarcheckpoint e marcar CENSURADO,
jamais chamar de70iterações completas. Ordem geométrica fixa do protocolo.
Registrar CPUfilhos/parent eCLOCK_MONOTONIC,manifesto/hashbinário,checkpoints.

Apósretorno,rota de melhorqualidade(validada) recebe LKH até2sADICIONAIS,
mesmasopçõesMAX_TRIALS100. Relatar ganho de pósprocessamento com esse custo
adicional; nãoalegarcomparaçãoigualtempo ou trabalhosidênticos entreMDM eLKH.
Objetivo: medirfolga após baseline maisconvergida e verificar se ganhos curtos
persistem;nenhumgate parateste. Preservartodoscasos ecompararvalorespublicados
somentecomo referência descritiva,semclaim de tempo crosshardware ou ótimo.
