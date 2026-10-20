# Instruções para qualquer assistente que retome este projeto

Leia primeiro `diario.md` e `research_notes/CONTINUATION.md`. Depois confira
`CLAUDE.md`, `memory/projects/pcctp.md` e `research_notes/EXECUTION_PLAN.md`.
O diário explica decisões e evidências; CONTINUATION registra o estado
operacional. Confirme esse estado nos arquivos e processos antes de agir.

Durante todo o projeto, atualize `diario.md` após cada avanço material e antes
de interromper ou entregar o trabalho a outra ferramenta. Registre data/hora
e fuso, o que fez, por quê, hipótese, arquivos/comandos relevantes, evidências,
limitações, falhas e próximo passo concreto. Distingua planejamento, teste
sintético, experimento real parcial e resultado real completo. Preserve as
entradas anteriores; corrija erros com uma nota explícita. Não invente uma
cronologia para preencher lacunas. Mantenha o diário legível para uma pessoa
ou IA sem acesso ao chat e atualize seu resumo de retomada quando necessário.

Trabalhe na branch `codex/pcctp-learning-submission`. Preserve a alteração
preexistente do usuário em `include/Utils.h`. Não implemente CTP. Não abra as
12 geometrias de teste reservadas durante desenvolvimento ou seleção.
Não inicie processos duplicados; verifique workers, locks e filas antes de
lançar otimização. Não altere fontes congeladas de experimentos em execução.

O usuário autorizou investigação autônoma e aceita combinações de métodos da
literatura para melhorar resultados. ML é opcional. Preserve resultados
negativos, cobrança justa do tempo total e separação por geometria. O objetivo
é contribuição sustentada e artigo completo, auditado e compilado; um rascunho
ou efeito parcial não permite declarar conclusão ou garantir publicação.

O usuário pediu, em 20/10/2026, acompanhamento do experimento em curso e um
texto de situação pelo menos a cada 30 minutos. Grave no diário o que mudar
a leitura científica ou o plano operacional. Faça commit nessa branch quando
houver avanço relevante de método, resultado auditado, falha que mude o
plano, ou decisão sobre o que entra no artigo. A mensagem deve dizer o que
aconteceu e qual foi o avanço. Não faça commit de conferência sem mudança,
de `include/Utils.h`, de `*.queue.lock`, nem de arquivo que um processo
ativo ainda esteja escrevendo. A mesma regra está em
`.cursor/rules/commitar-avanco-relevante.mdc`.
