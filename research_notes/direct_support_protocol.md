# Viabilidade de construção rápida + busca de suporte sem MDM

19/10/2026, antes das execuções. Motivação: consumo de10s na construçãoMDM
limita saldo da busca. Diagnóstico distinto: construçãofallbackcomum verificada
seguida de support_lkh, orçamentoTOTAL15s incluindo leitura/construção/busca.
Todos9hardTRAIN, determinístico, nenhuma instância reservada. Mesmos candidatos,
top8,ordenaçãoLKH0.2s,aceitaçãoestrita; sem tuning. Registrartrajetórias globais.

Comparar descritivamente custos com médiasMDM15s jácoletadas, explicitando
execuções em momentosdiferentes e sem inferência pareada de velocidade.
Gate apenas para justificar novo experimento justo: custo menor que médiaMDM
em>=6/9 e ganhomédio>=0.5%,0falhas. Não liberar teste nem alegar confirmação.
Se falhar, preservar e não ajustar candidatos/tempo comvalidação.
