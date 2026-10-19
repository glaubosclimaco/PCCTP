# Hipótese seguinte preparada: recombinação inteira em união de arestas

19/10/2026, antes de qualquer coleta PCCTP desta hipótese. Não iniciar enquanto
MDM70 estiver ativo. Objetivo: superar a limitação de fixar visitas/arestas locais,
permitindo recombinar rotas viáveis com diferentes suportes. Sem alegar novidade
na ideia clássica de recombinação; literatura específica ainda deve ser conferida.

Todos9hardTRAIN, após concluir o lote70 (incluindo censurados). Para cada caso,
unir arestas de TODAS as rotas native+LKH de fixed_work e das seis rotas salvas
em support_common_fallback: modosMDM/fixed_lkh,seeds42/43/44. Validar cada rota,
deduplicar sem filtrar por qualidade e iniciar da melhor disponível. Todas visitas
opcionais livres; somente arestas fora da união ficam fixadas emzero. Isso mantém
incumbente viável. ModeloSCF existente com highspywarmstart,threads1,seed42,
60s totais por caso incluindo leitura/união/modelo. Preservar fallback e descartar
melhorias retornadas fora do prazo. Medir CPU e os dois relógios por duração.

É pilotoOFFLINE: custo de produzir todas as rotas é anterior e não está nesse60s.
Não comparar como algoritmo completo com MDM60s. Limites/status ótimo pertencem
SOMENTE ao modelo restrito à união, nunca assumir boundPCCTPglobal.
Gatepotencial: melhora>=3/9sobre melhorrotadeentrada e ganhomédio>=0.2%,
nenhuma perda de factibilidade. Sem tunar em validação e sem abrir teste.
Antes da coleta: verificar modelo em grafos pequenos por enumeração independente,
auditando arestas e rejeição de partida fora do pool.
