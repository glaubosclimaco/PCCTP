# Expansão de arestas: critério não atingido

Nove geometrias hardTRAIN,18 execuções; auditoria independente de36 rotas
PASS. A auditoria reconstrói ambos os conjuntos de candidatos usando a matriz
crua, verifica hashes, custos, factibilidade e regra de descarte tardio.

Nearest melhora2/9 entradas, ganho médio0.2015418%, vence random2/9. O
critério exigia três melhorias e cinco vitórias sobre random; FALHOU.
gr22996791→96707 iguala referênciaMDM; si1757006→6885 ainda acima do
alvo6871 já conhecido. Nenhum novo recorde sustentado. gil262 nearest teve
melhoria tardia descartada conforme protocolo; não substituir o custo aceito
pelo resultado bruto ao calcular desempenho.

Ampliação não é monotônica em desempenho sob tempo limitado: embora o ótimo
do grafo maior não possa piorar, o solver pode não recuperar uma melhoria do
pool puro no prazo. Portanto, não inferir ganho do algoritmo pelo conjunto
viável maior. Não abrir validação/teste com esses critérios falhos.

Próxima decisão científica: antes de outra coleta, comparar os gargalos de
qualidade e de tempo por geometria entre as hipóteses já executadas. Evitar
continuar aumentando o grafo por tentativa e erro sem controle de custo.
