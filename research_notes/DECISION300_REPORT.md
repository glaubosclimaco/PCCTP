> **LOTE INTERROMPIDO: SUSPENSÃO DO COMPUTADOR.** Os dados abaixo são históricos e não sustentam comparação temporal controlada. Veja `ruin_total300_train/environment_interruption.json` e `awake300_repetition_protocol.md`.

# Avaliação experimental PCCTP: orçamento total de 300 segundos

Estado: **INCOMPLETE**. 43/81 execuções auditadas; 4/9 geometrias completas.

Nenhuma conclusão de publicação, significância ou generalização é automática.
Dados de treino selecionados por dificuldade; teste reservado permanece fechado.

| Método contra MDM300 | Ganho médio (%) | Vitórias | Perdas | Critério |
|---|---:|---:|---:|---|
| hybrid10 | 1.165714 | 4 | 0 | pendente |
| hybrid60 | 1.209160 | 4 | 0 | pendente |

## Custos e consistência nas geometrias completas

| Instância | Média MDM | Média H10 | Média H60 | H10: sementes melhores/iguais/piores | H60: sementes melhores/iguais/piores |
|---|---:|---:|---:|---|---|
| g_a280_VT56_T168_W56_50 | 2029.333 | 1998.000 | 2001.333 | 3/0/0 | 3/0/0 |
| g_gil262_VT51_T158_W53_75 | 2016.333 | 1982.667 | 1983.667 | 3/0/0 | 3/0/0 |
| g_gr229_VT45_T138_W46_50 | 98218.000 | 97060.000 | 96722.333 | 2/0/1 | 3/0/0 |
| g_kroA200_VT40_T120_W40_25 | 23569.000 | 23505.333 | 23495.000 | 2/0/1 | 3/0/0 |

As contagens por semente são descritivas dentro da mesma instância;
não multiplicam o número de geometrias independentes.
Geometrias com falha de factibilidade: 0.

O desfecho primário é a média das reduções relativas das médias por
geometria. Curvas em outros tempos são descritivas e não substituem300s.
As contagens acima referem-se apenas a geometrias completas; resultados
parciais não fundamentam aprovação ou rejeição do método.

| Método | Execuções | Mineração observada | Mineração pontual comprovada | Limite superior desconhecido |
|---|---:|---:|---:|---:|
| hybrid10 | 15 | 0 | 0 | 0 |
| hybrid60 | 14 | 0 | 0 | 0 |
| mdm | 14 | 1 | 1 | 0 |

Mineração observada em export final pode incluir sobretempo cooperativo.
A prova de ativação pontual usa export/checkpoint anterior ao orçamento nativo.
Ativação não demonstra eficácia da mineração. Tempo de CPU do processo pai
não deve ser apresentado como CPU total: LKH e mineradores são subprocessos.
A comparação principal usa tempo de parede debitado e incumbentes pontuais.

Este orçamento não reproduz as70iterações dos artigos. A decisão sobre
publicação ou abandono exige confronto com a referência completa;
published70_calibration_protocol.md define a calibração inicial, que
também não substitui a avaliação final com dez sementes.

Se o filtro de treino passar, a próxima etapa é desenhar confirmação
prospectiva mais ampla, com parâmetros congelados e controles pertinentes.
Se ambos falharem, preservar o diagnóstico e reavaliar essa linha antes
de consumir o teste reservado ou continuar a redação como se houvesse ganho.
