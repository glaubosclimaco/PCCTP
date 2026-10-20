# Identificadores de dados: pendência antes de afirmar réplica completa

A página primária do artigo RAIRO2021 lista nas referências5 e6 os registros
Mendeley10.17632/8yjbfgcfvn.1 (dataset2020) e10.17632/8yjbfgcfvn.2
(resultados2021). O registro2019 antes inspecionado tem outro identificador,
10.17632/gpnmh9kdr9.1, e sua descrição diz resultados completos de formulações
e heurísticas. As páginas não demonstram igualdade dos arquivos entre depósitos.

Não citar o registro2019 como prova de que os nove arquivos locais reproduzem
todo o benchmark publicado. A comparação atual usa variantes locais identificadas
por nome e hash, com fonte recuperada fornecida pelo usuário. Antes de alegar
replicação dos107 casos, reconciliar lista de instâncias, versões, parâmetros,
arquivos e resultados do depósito efetivamente citado pelo artigo. Não abrir os
doze testes reservados para resolver essa tarefa de metadados.

Fontes primárias consultadas em20/10/2026:
- https://numdam.org/articles/10.1051/ro/2021057/
- https://www.numdam.org/item/10.1051/ro/2021057.pdf
- https://data.mendeley.com/datasets/gpnmh9kdr9/1

A consulta direta ao depósito8yjbfgcfvn/1 falhou nesta ferramenta; não se obteve
sua lista de arquivos. Trata-se de pendência de identidade, não de erro comprovado
nos dados locais. Nenhum dado reservado foi baixado ou inspecionado.

## Atualização 20/10 — depósito oficial acessível no navegador; 3/9 TRAIN iguais

A página v2 foi lida na web e no navegador: publicada em 9/3/2021, descrição de
144 instâncias, DOI10.17632/8yjbfgcfvn.2, licença CC BY4.0. O índice contém
arquivos individuais, além de pastas; o contador de147 arquivos da interface
não deve ser interpretado como147 instâncias. A referência ao conjunto de107
casos é o subconjunto de avaliação após37casos de ajuste, conforme o manuscrito
MDM fornecido;144 e107 não são contagens contraditórias.

A consulta HTTP direta retornou403 e v1 continuou inacessível pelo web tool;
isto não impediu ler os painéis públicos de informação dos arquivos TRAIN em v2.
Os SHA-256 públicos de a280, gil262 e gr229 coincidiram byte a byte com os
arquivos locais, também com seus hashes no perfil independente. Registro com
nome/UUID/URL/hash em benchmark_metadata/mendeley_train_identity.json;
auditor restrito a TRAIN em experiments/lkh/audit_benchmark_train_identity.py.
Estado PARTIAL_TRAIN:3verificados,6pendentes. Nenhum download ou pasta de
resultados foi aberto, nem conteúdo das12geometrias reservadas.

Isso resolve a identidade desses três arquivos com o depósito v2 efetivamente
citado em RAIRO2021. Não prova equivalência das versões1/2, igualdade do depósito
2019gpnmh9kdr9, nem réplica dos107casos. Conferir os seis TRAIN restantes antes
de ampliar a afirmação. O PDF17p atual não foi alterado nesta retomada.

## Atualização seguinte — identidade9/9 TRAIN comv2 concluída

Os seis painéis restantes foram localizados pelo nome no índice público e seus
SHA-256 lidos no campo File information: kroA200, pr226, rd400, si175, ts225 e
tsp225. Todos coincidem com os bytes locais e com o perfil independente. O auditor
agora informa COMPLETE_TRAIN,9verificados,0pendentes. Permanecem fora desse
escopo os outros135arquivos de instância, a composição exata do subconjunto107,
os resultados históricos e a equivalência de hashes entre versões1/2.

A página pública de comparação entre versões foi também acessada:
https://data.mendeley.com/datasets/compare/8yjbfgcfvn
Exibe título, descrição144 e DOI de v1 e v2; a lista é virtualizada e não expõe
checksums nessa tela. Não inferir igualdade de todos os arquivos com base em
nomes/tamanhos. A informação de dia de publicação de v1 aparece16/5 no link da
página v2 e17/5 no comparador; conservar o ano2020 sem arbitrar essa diferença.
Nenhum arquivo adicional ou conteúdo reservado foi aberto. Não houve download.
