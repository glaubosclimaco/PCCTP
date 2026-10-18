# Avaliação de validação congelada para MDM9s+LKH

18/10/2026, antes de abrir resultados reservados. Executar SOMENTE se o gate
pré-fixado do piloto TRAIN passar após auditoria completa. Não usar resultados
parciais como autorização de execução. Modelo9/1,versão/opções LKH,budget10s,
seeds42/43/44 e baselineMDM permanecem congelados.

Incluir TODAS as7geometrias com role=validation no protocolo original
decoder_initialized_states/protocol.json e suas instâncias já escolhidas nesse
protocolo. Não filtrar por tamanho ou dificuldade e não abrir role=test.
Ordem de braços reprodutível por instância+seed,mesma implementação temporal
do piloto e inicialização debitada. Não aproveitar incumbentes do treino.

Sem tuning após validação: gate para teste final é ganho médio positivo por
geometria,melhor média em>=4/7geometrias,pior em<=1/7 e nenhuma perda de
factibilidade. Relatar magnitude/incerteza e todos casos; o gate não substitui
evidência científica. Comparação pareada agregada por geometria;não confundir
21sementes com21instâncias independentes. Nenhuma afirmação de novo recorde
antes de confrontar literatura e nenhuma declaração de conclusão do artigo.
