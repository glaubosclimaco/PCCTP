# Comparação exploratória por tempo total de 15 segundos

19/10/2026, antes de executar. Todos9hardTRAIN, sementes42/43/44; quatro braços:
MDM15s, MDM10s+LKHfixo, MDM10s+suporte2opt, MDM10s+suporteLKH.
Mesmas implementações/opções do diagnóstico; o pós-processamento recebe apenas
o saldo real até15s, incluindo inicialização, leitura e retorno nativo. Aceitar
melhorias apenas antes do prazo externo; preservar incumbente nativo válido
em caso de excesso. Registrar tempos reais e falhas. Ordem aleatória reproduzível
por instância/semente. Não usar incumbentes em cache; cada braço começa do zero.

Gate exploratório para justificar estudo adicional: ganho médio entregeometrias
>=0.5% frenteMDM15s,>=6/9vitórias e<=1derrota, nenhuma inviabilidade; vantagem
sobre fixoLKH em>=3/9 e sobre suporte2opt em>=3/9. Relatar todos resultados,
sem tratar27sementes como27instâncias. Não libera teste reservado; validação
anterior do método distinto9+1 falhou3/7. Não alegar recorde:si1756871 já estava
no campo valor_otimo da entrada; auditado independentemente com Pythonescalar.
