# Cache de serializacao: resultado completo e decisao editorial

Registro em 22/10/2026 11:29 Fortaleza (UTC-3). Experimento real completo.
Fonte: experiments/lkh/cached30_train_audit.json, status PASS, batch_complete
true; protocolo SHA24247ef583a90b6b76e9bb03166d2a7cead1c01f323a836f0e05b80cda43137a.
O desenho prospectivo cached30_conditional_design.md permanece historico,
sem reescrever seu status anterior a coleta.

Foram auditados todos os27pares/54fases e501verificacoes de rotas. Cada
braco recebeu30s totais desde setup; construcao lazy do cache foi cobrada.
As27rotas iniciais sao reutilizadas de TRAIN; sua geracao nao foi cobrada.
Portanto esta comparacao mede apenas a fase condicional de reparacao.

O estimando predefinido compara medias de custo por geometria (tres estados)
e depois faz media equiponderada das nove reducoes percentuais. Resultado:
400/10607 pontos percentuais =0,03771094560196097%; uma geometria positiva,
zero negativas e oito empates. rd400: media uncached10607 e cached10571,
reducao3600/10607%=0,33939851041764874%. As outras oito medias empatam.
O gate prospectivo exige media>=0,2%, pelo menos seis positivas e no maximo
uma negativa, com coleta completa/valida: FAIL. Nao selecionar apenas rd400.

Decisao: preservar como controle negativo de beneficio de qualidade. Nao
integrar cache na coleta paired70 em andamento nem expandir sementes para
resgatar este gate. Equivalencia dos bytes e velocidade de serializacao nao
estabelecem melhoria de solucao. Este resultado tampouco demonstra ausencia
de beneficio em todo outro regime: limita-se a estes estados/orcamento.

Para o artigo: relatar desenho, contabilizacao, denominador completo, efeito
e gate negativo em uma analise de componentes. Nao apresentar como speedup
end-to-end, confirmacao independente, superioridade ao MDM70 ou teste final.
Nao ha inferencia populacional/p-valor; as tres sementes sao repeticoes dentro
de geometria. Os dados completos e negativos continuam acessiveis no auditor.
