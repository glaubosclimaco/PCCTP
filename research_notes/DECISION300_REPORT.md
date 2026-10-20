# Avaliação experimental PCCTP: orçamento total de 300 segundos

Estado: **INCOMPLETE**. 9/81 execuções auditadas; 1/9 geometrias completas.

Nenhuma conclusão de publicação, significância ou generalização é automática.
Dados de treino selecionados por dificuldade; teste reservado permanece fechado.

| Método contra MDM300 | Ganho médio (%) | Vitórias | Perdas | Critério |
|---|---:|---:|---:|---|
| hybrid10 | 1.544021 | 1 | 0 | pendente |
| hybrid60 | 1.379763 | 1 | 0 | pendente |

## Custos e consistência nas geometrias completas

| Instância | Média MDM | Média H10 | Média H60 | H10: sementes melhores/iguais/piores | H60: sementes melhores/iguais/piores |
|---|---:|---:|---:|---|---|
| g_a280_VT56_T168_W56_50 | 2029.333 | 1998.000 | 2001.333 | 3/0/0 | 3/0/0 |

As contagens por semente são descritivas dentro da mesma instância;
não multiplicam o número de geometrias independentes.
Geometrias com falha de factibilidade: 0.

O desfecho primário é a média das reduções relativas das médias por
geometria. Curvas em outros tempos são descritivas e não substituem300s.
As contagens acima referem-se apenas a geometrias completas; resultados
parciais não fundamentam aprovação ou rejeição do método.

| Método | Execuções | Mineração observada | Mineração pontual comprovada | Limite superior desconhecido |
|---|---:|---:|---:|---:|
| hybrid10 | 3 | 0 | 0 | 0 |
| hybrid60 | 3 | 0 | 0 | 0 |
| mdm | 3 | 0 | 0 | 0 |

Mineração observada em export final pode incluir sobretempo cooperativo.
A prova de ativação pontual usa export/checkpoint anterior ao orçamento nativo.
Ativação não demonstra eficácia da mineração. Tempo de CPU do processo pai
não deve ser apresentado como CPU total: LKH e mineradores são subprocessos.
A comparação principal usa tempo de parede debitado e incumbentes pontuais.

Se o filtro de treino passar, a próxima etapa é desenhar confirmação
prospectiva mais ampla, com parâmetros congelados e controles pertinentes.
Se ambos falharem, preservar o diagnóstico e reavaliar essa linha antes
de consumir o teste reservado ou continuar a redação como se houvesse ganho.
