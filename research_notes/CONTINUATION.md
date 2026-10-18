# Estado operacional para retomada autônoma

## Retomada às 5h11 — atualização prioritária

Esta seção substitui o estado antigo dos processos abaixo.

- As coletas originais terminaram: 29 primeiras variantes, um timeout sem
  primeira solução; 58 outras variantes, três SIGSEGV. As sementes 43/44 nas
  nove primeiras variantes de desenvolvimento terminaram com duas SIGSEGV.
- Auditoria `experiments/ml_lns/dataset_audit.json`: 99 baselines, 509 snapshots,
  151 estados rotulados e 3384 reparos verificados; 29 geometrias; nenhum rótulo
  reservado de validação/teste. Arquivos têm hashes. Falhas fora do corpus
  continuam nos manifestos, não foram convertidas em soluções válidas.
- `strong_oracle.py` rodada 5 (sessão 87985) terminou. Coletou as variantes
  disponíveis de seed42. Sementes43/44 ainda precisam de coleta específica.
- `nonlinear_ranker.py` e `evaluate_nonlinear.py` implementam controles lineares
  e random Fourier features (largura64, alpha1, seed2026), escala apenas no treino,
  ponderação igual por geometria informativa. Resultado ainda exploratório,
  sem vantagem assegurada frente à regra `joint_savings`. Não usar teste final.
- Exportações agora são publicadas atomicamente somente após validação. Novas
  execuções exportam tempo por snapshot e tempo total do processo. NÃO imputar
  tempos aos snapshots antigos que não têm esse campo.

### Falha do legado localizada e correção isolada

O decodificador em `mineracao.cpp` declara `int cc_atual;` sem inicializar e o
usa ao criar componente de padrão com uma única aresta. Relatórios locais de
crash apontaram `decode`; teste mínimo `verify_decoder.py` passa com inicialização
em zero. A versão corrigida é separada, nunca uma substituição silenciosa:

```sh
python3 experiments/original_mdm/build.py --initialize-decoder-index
python3 experiments/original_mdm/verify_decoder.py
```

Diretório `build_decoder_initialized/`, método JSON
`original_MDM_legacy_decoder_index_initialized`. Fontes fornecidos intactos.
Os três casos de seed42 que antes falhavam concluíram e passaram na validação:
eil51_VT9_T11_W31_50 custo248, hk48_VT9_T10_W29_75 custo7803,
rd100_VT20_T60_W20_50 custo6626. Controle brazil preservou20466.
Resultados em `experiments/original_mdm/diagnostics/decoder_initialized_reproduction.json`.
Não inferir daí equivalência universal com execuções de comportamento indefinido.

ASAN/UBSAN com clang concluiu dois casos sem reportar erro e atingiu timeout no
terceiro; isso não invalida o defeito de inicialização nem certifica ausência
de outros bugs. GCC não dispõe de libasan local. LLDB foi encerrado depois da
identificação por inspeção e teste mínimo; não há depurador útil a manter.

### Processo atualmente em execução

`strong_states.py --roles development train --seeds 42 --initialize-decoder-index`,
log `/tmp/pcctp_decoder_initialized_states.log`. Coleta as primeiras variantes
com a correção EXPLICITAMENTE identificada; arquivos em `decoder_initialized_states/`.
Não misturar com `original_states/` nem reaproveitar seus rótulos como se fossem
do mesmo algoritmo. Split por geometria foi copiado, preservando reservas.

O coletor aceita banco separado e escolhe saída separada automaticamente:

```sh
python3 experiments/ml_lns/strong_oracle.py --state-dir experiments/ml_lns/decoder_initialized_states
```

Antes de comparar algoritmos novos, estabilizar essa referência e verificar
consistência da correção nas outras instâncias. O teste final continua intocado.
Próximo estudo necessário: qualidade por tempo desde o primeiro estado do MDM,
comparado à continuação do próprio MDM, ALNS sem aprendizado e regra simples.

Atualizado em 18 de outubro de 2026. Não é um relatório de contribuição concluída.

## Decisões vigentes

- Branch `codex/pcctp-learning-submission`; commits `f0cc6a6` e `ac048e4`
  registram memória, esqueleto editorial e pesquisa acumulada.
- Usuário pediu explicitamente não aguardar aprovação e prosseguir até paper
  completo/revisado. Autorizou retomadas a cada 10 minutos durante sua ausência.
- Heartbeat ativo `avan-ar-pesquisa-e-artigo-pcctp`, a cada 10 minutos neste chat.
  Não criar duplicata. Desativar apenas ao cumprir os critérios ou pedido do usuário.
- CTP continua fora do escopo atual. `include/Utils.h` é modificação preexistente
  do usuário, não incluir em commits da pesquisa.

## Processos disparados nesta etapa

Os IDs abaixo são sessões de exec da conversa; se indisponíveis na retomada,
inspecione os arquivos de saída e logs antes de iniciar outra execução.

- 25250: `strong_states.py --roles development train --seeds 42`, retomada dos
  casos inicialmente interrompidos. Log `/tmp/pcctp_original_states_round2.log`.
  Arquivos completos ou checkpoints validados em `experiments/ml_lns/original_states/`.
  Manifesto final `last_generation.json`.
- 1342: `strong_states.py --roles development train --seeds 42 --other-variants`.
  Log `/tmp/pcctp_original_other_variants.log`; duas execuções paralelas.
  Coleta TODAS as outras variantes das geometrias permitidas, não apenas as
  que pareçam favoráveis. Manifesto `other_variant_generation.json` ao terminar.
  Não colide com a primeira seleção: instâncias possuem nomes distintos.
- 63842: `strong_oracle.py`, log `/tmp/pcctp_strong_oracle_round3.log`.
  Só enxerga os arquivos disponíveis no início. Reexecutar incrementalmente
  quando terminarem novas instâncias, sem iniciar dois coletores para os mesmos
  estados simultaneamente. Arquivos `.partial` não são resultados completos.
- 40974: diagnóstico de raios com backend highspy, warm start e cortes externos;
  log `/tmp/pcctp_radius_cuts.log`. Provavelmente já concluído na retomada.
- Rodadas anteriores de raio com scipy, highspy frio/quente e busca iterativa
  de 3 geometrias são curtas e já produziram seus JSONs. Não repetir sem motivo.

## Verificações realizadas

- `verify_signature.py`: 60 casos aleatórios, 960 atribuições, dez descritores.
- `verify_fpmax.py`: padrões maximais de seis limiares coincidiram com enumeração.
- `checkpoint_verification.json`: mesma semente preservou rota e custo original;
  timeout artificial reteve solução viável validada independentemente.
- `verify_fragment.py`: 37 casos, incluindo 15 por enumeração e 20 contra SCF
  completo; dois casos liberando apenas vértices ausentes.
- `verify_native.py`: reproduz 60 comparações da formulação com cortes externos
  e highspy, com/sem solução inicial; resultado salvo passou, 35 cortes adicionados.
- `fixed_path_repair.md`: argumento de equivalência inteira e limite de até 4k
  componentes para liberação incidente em k opcionais mais 2k arestas. Revisar
  formalmente antes de entrar no artigo; não alegar novidade da contração.

## O que os resultados dizem até aqui

- Piloto anterior sem MDM mostrou vantagem do contexto de rota, mas não sustenta
  vantagem científica frente à heurística original.
- MDM original foi recuperado. A primeira coleta teve nove timeouts em 29 casos;
  agora checkpoints impedem a perda total das observações. Não excluir geometrias.
- Em `strong_development_results.json`, o ranking linear com descritores conjuntos
  supera o genérico em uma avaliação exploratória, mas a regra simples
  `joint_savings` permanece competitiva ou superior. Ainda não há valor adicional
  de ML comprovado nem significância confirmatória.
- Busca iterativa adicional de 10 s em att48, berlin52 e brazil58 após 70
  iterações do MDM: nenhuma das quatro políticas melhorou o custo.
- Raios 8/16/24 em gr137, pr144, u159 também não melhoraram nos diagnósticos.
  Nos artigos fornecidos vários desses valores coincidem com ótimos reportados;
  igualdade numérica deve ser distinguida de certificado computacional próprio.
  Nessas instâncias é relevante reduzir tempo até o ótimo, não melhorar o custo.
- Contração não assegura aceleração: tempos preliminares são mistos e rodados
  sob carga concorrente. Comparações finais devem usar execução controlada.
- Gurobi instalado tem licença expirada; continuar com HiGHS disponível sem
  aguardar usuário, comprar licença ou alterar configuração de licença.

## Próximos passos concretos

1. Verificar logs/manifestos e ausência de processos duplicados. Aproveitar as
   variantes adicionais já concluídas para coleta incremental de rótulos.
2. Auditar os resultados e medir também tempo até alvo e integral de qualidade.
   Comparar desde inicialização de mesmo orçamento, não atribuir a ML uma
   melhora que o próprio MDM alcançaria mais rapidamente continuando sua busca.
3. Ampliar rótulos em estados e variantes difíceis de treino/desenvolvimento.
   Testar modelos não lineares com os MESMOS controles de capacidade e regras
   simples. Evitar acrescentar uma rede apenas para dar aparência de novidade.
4. Integrar política completa ao fluxo GRASP/MDM; controles ALNS sem aprendizado
   e MDM sob tempo igual ainda faltam. O script `iterative_study.py` é um piloto,
   não a avaliação final. Seu modelo é congelado com hashes antes da trajetória.
5. Se ML não superar os controles, registrar refutação e investigar alternativa
   específica, mantendo validação/teste reservados e verificando literatura.
6. `submission/main.tex` é apenas esqueleto. Ainda faltam contribuição sustentada,
   protocolo final, experimento confirmatório, bibliografia auditada e texto.
   Não preencher vinte páginas com alegações não demonstradas.

## Comandos úteis

```sh
python3 experiments/ml_lns/strong_oracle.py
python3 experiments/ml_lns/evaluate_strong.py
python3 experiments/ml_lns/verify_signature.py
python3 experiments/ml_lns/verify_fragment.py
python3 experiments/ml_lns/verify_native.py
python3 experiments/ml_lns/radius_probe.py --backend highspy --warm-start --connectivity cuts
```

Antes de escrever o manuscrito: aplicar scientific-writing (o SKILL.md instalado
é um ponteiro textual para outro arquivo; localizar a fonte real), usar fontes
primárias verificadas e copiar/aplicar template editorial de forma rastreável.
Ao gerar PDF, seguir skill de PDF, compilar e inspecionar todas as páginas.
