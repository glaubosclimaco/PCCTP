# Seleção aprendida: triagem interna por geometria

19/10/2026, antes de treinar. OráculoreparaçãohardTRAIN passou6/9,0.50094%.
Noveestados,9candidatoscada: amostra pequena, não confirma valorendtoend.
Ridgefixo(alpha1), descritores existentes pilot.features + onehotrandom/coverage/
prize, alvo melhoriarelativapercentual do custo. Leave-one-geometry-out: escalar
ajustado somente8geometriasdetreino porfold,semseleçãohiperparâmetros. Predizer
únicocandidato melhor emgeometriaexcluída, sem olharseurótulo paraescolher.
Comparar ganhoobservado comexpectativauniformeTODOS e expectativauniforme
apenasaleatórios (controleforte),ambosumsubproblema. Não compararcomoracle9
como se mesmocusto. GateML:>=0.1pontopercentualdeganhomédio sobreALEATÓRIO,
>=6/9vitórias contraexpectativarandom;senão nãointegrar seletor.
Relatartodosfolds,empates,negativos. SomenteTRAIN,múltiplashipótesesexploratórias,
nãoestatísticaconfirmatória e não libera teste. Tempo deextração/treinooffline
nãoésubstituto de compararalgoritmocompleto cobrandoinferência/reparaçãoonline.
