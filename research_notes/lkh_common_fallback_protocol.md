# Quatro braços com incumbente externo comum

19/10/2026, antes das execuções. O estudo anterior falhou por ausência de
incumbentes dentro do prazo; preservar integralmente os resultados. Novo estudo
exploratório em TODOS9hardTRAIN, sementes42/43/44, mesmos4braços e15s totais.

Antes do MDM, TODOSconstroem cicloNN sobre todosvisitáveis e removem opcionais
factíveis com capsuave0.25s, exatamente feasible_fallback.py verificado.
Leitura,construção e validação cobradas desde relógio externo. Registrar rota,
custo,tempo de aceitação e origem do incumbente. Reter melhor entre fallback
aceito até15s e incumbente nativo aceito até seu prazo(15sMDM,10soutros).
Pós-processar esse melhor incumbente usando somente saldo externoaté15s.
Sem passar fallback para rotina internaMDM; é uma solução de segurança externa.
Retorno do processo pode exceder prazo; registrar e não aceitar melhoria tardia.

Mesmos top8,opçõesLKH,2opt e ordem reprodutível do estudo anterior. Não mudar
configurações após observar resultados. Gate exploratório: ganho médio por
geometria>=0.5%vsMDM,>=6wins,<=1loss,0falhas de factibilidade,>=3winsvsfixo e
>=3winsvs2opt. Relatar resultados individuais,origens e casos que só usamfallback.
Não comparar como se fossem emparelhados com batch anterior sob carga distinta.
Não libera teste reservado; nenhuma alegação de recorde ou conclusão do artigo.
