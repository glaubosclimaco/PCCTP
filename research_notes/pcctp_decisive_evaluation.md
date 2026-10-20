# Decisão sobre investir em um artigo PCCTP

O usuário reafirmou em19/10 que deseja experimentos reais para avaliar se vale
escrever o artigo. Manter foco no PCCTP; a triagem de outros problemas deixa de
ser prioridade enquanto este estudo roda. Não trocar por redação especulativa.

## Estudo em execução

81 execuções novas: nove instâncias de treino difíceis, três sementes e três
métodos,300segundos totais cada. Controle MDM recuperado com correções declaradas;
híbridos com10ou60segundos de MDM antes da busca. Mesmo fallback, inicialização
debitada, seleção determinística da ordem, apenas um otimizador de cada vez.
Não usar rotas offline como partidas nem permitir incumbentes após o prazo.
Todos os parâmetros e o critério primário já foram fixados antes da coleta.
Não alterá-los com base em resultados parciais.

As regras formais estão em ruin_long_horizon_protocol.md. Avaliar separadamente
os dois híbridos contra o controle: média por geometria>=0.2%, pelo menos6/9
vitórias, no máximo1/9perda e nenhuma regressão de factibilidade. Esse é um
filtro exploratório de investimento, não uma garantia de publicação.

O relatório deve conter custos de todas as sementes, médias por geometria,
perdas/empates, distribuições de efeito, curvas descritivas nos tempos1,10,30,
60,120,300segundos, fallback/timeouts e duração real dos processos. O desfecho
primário continua300segundos; não escolher outro ponto porque parece melhor.
Incluir quantas execuções realmente fizeram mineração e quais ativações têm
tempo comprovado. Contagem cumulativa sem timestamp não prova ativação pontual.

## Interpretação e decisão posterior

Se um híbrido passar, há motivo para continuar com uma avaliação prospectiva
mais ampla e controles fortes. Ainda será necessário investigar diversidade de
instâncias/quota/cobertura, ativação da mineração, custo total e comparação
compatível com o algoritmo publicado. Reservar dados independentes antes de
qualquer confirmação. As sete geometrias de validação já foram utilizadas e
não se tornam independentes novamente por trocar o orçamento ou as sementes.
As doze de teste continuam fechadas sob o protocolo atual.

Se ambos falharem, não expandir repetidamente essa mesma configuração esperando
um resultado favorável. Preservar o relatório e decidir entre uma hipótese
estrutural distinta com justificativa concreta ou migrar de problema. Não
interpretar um gate falho como prova de que nenhuma pesquisa PCCTP é possível.

Os ganhos curtos já observados e o MIP negativo serão relatados integralmente.
Nenhum novo recorde ou superioridade ao estado da arte será declarado com base
em valores de cabeçalho sem verificar origem, definição e comparabilidade.
O rascunho de artigo não substitui essa avaliação; a decisão de concluir a
submissão depende de evidência suficiente, e não do número de páginas escrito.
