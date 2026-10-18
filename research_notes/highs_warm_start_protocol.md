# HiGHS nativo: início MDM e captura de incumbentes

18/10/2026, antes da comparação. Usar highspy1.13.1 em ambos os braços,
formulação SCF idêntica de Repair,11TRAIN completas, semente42,1thread,
10s totais por chamada incluindo leitura, montagem e configuração. Importação
comum informada separadamente. Braços: solver sem start e solver com start MDM
completo, incluindo fluxos. Nos dois, incumbent MDM é fallback externo.
Não comparar a versão SciPy com a versão highspy para atribuir efeito ao start.

Verificar o vetor inicial contra limites, integralidade e todas as linhas.
Callbacks de solução melhor validam ciclos no PCCTP original e guardam somente
incumbentes cuja validação termine antes do prazo. Captura de incumbente deve
funcionar mesmo quando o solver retorna após10s. Registrar atraso, nós e limites;
limites da SCF completa são globais, mas só declarar ótimo se certificado.

Pré-requisito: instância sintética com enumeração completa e teste explícito
de entrega de callback, em ambos os braços. Critério para comparação maior:
start melhora MDM em pelo menos2/11 geometrias. Se falhar, não aumentar prazo
ou trocar instâncias neste lote. Se passar, comparar por tempo TOTAL contra
MDM com tempo extra, solver sem start e híbrido em sementes pareadas. Não
atribuir novidade ao MIP start. Validação e teste reservados fechados.

Documentação primária do start: https://github.com/ERGO-Code/HiGHS/blob/master/docs/src/guide/further.md
Consulta da API instalada confirmou setSolution e callbacks de incumbentes.
Isso não prova suporte de cortes dinâmicos: embora o enum cbMipDefineLazyConstraints
exista, a discussão oficial3079 (20/06/2026) informa ausência de lazy/user-cut
callbacks utilizáveis naquele momento. Não implementar branch-and-cut usando
o nome do enum como evidência de funcionalidade.
https://github.com/ERGO-Code/HiGHS/discussions/3079

## Resultado auditado

0/11 melhorias warm e0/11 cold.3provas warm e2cold. Gate de ganho primal
falhou; não ampliar prazo para este mecanismo. Auditoria de rotas/bounds
PASS em highs_warm_start_audit.json. Infraestrutura validada não equivale
a contribuição de desempenho. Próxima hipótese deve ser distinta.
