# Diagnóstico de persistência do ganho no treino

18/10/2026, antes de executar orçamento maior. A validação10s teve ganho médio
0.89319%,3vitórias/4empates/0derrotas; FALHOU o gate de4vitórias. Não abrir teste,
não alterar esse gate nem usar validação para calibrar tempo/seleção de casos.

Hipótese diagnóstica: ganho de ordenar visitas por LKH persiste quando o MDM
recebe mais tempo. Repetir TODOS os mesmos9casos difíceis TRAIN, sementes42/43/44,
MDM30s contra MDM29s+LKH até1s restante. Mesmas opções LKH, mesma ordem pareada
reprodutível, inicialização cobrada, aceitação até30s e tempos reais registrados.
Nenhuma geometria reservada neste diagnóstico. É extensão exploratória temporal,
não tentativa de fazer passar o gate reservado nem autorização para abrir teste.

Relatar todos27pares, ganho médio agregado por geometria, vitórias/empates/perdas,
factibilidade, tempos e comparação descritiva com10s. Não escolher posteriormente
horizonte vencedor como resultado confirmatório. Depois avaliar se há mecanismo
e efeito suficientes para estudo mais abrangente e controles adicionais.
