# MDM original fornecido pelo usuário

Compilação isolada e instrumentada dos fontes em `mdm-grasp-genius-vnd/`.
Não confundir com `experiments/legacy_grasp/`, uma reconstrução anterior usada
apenas para diagnóstico. O algoritmo original faz 70 iterações e ativa FPMax
quando se estabiliza o conjunto elite. Suporte absoluto 2, estabilidade 5%.

```sh
python3 experiments/original_mdm/build.py
python3 experiments/original_mdm/verify_fpmax.py
python3 experiments/original_mdm/run.py data/g_brazil58_VT11_T12_W35_75 --seed 42 --output /tmp/mdm-brazil.json
```

`build/manifest.json` registra hashes dos fontes originais, compilador, comando
e alterações observacionais. As mudanças são exportação de soluções/checkpoints,
contador de mineração, argumento de saída, expansão do buffer de caminho e troca
de `malloc.h` por `cstdlib` no FPMax para compilar no macOS. Fontes originais
permanecem intactos. Requer compilador GNU C++; foi usado `g++-15`.

Cada execução usa diretório privado em `runs/` para impedir colisões dos arquivos
temporários do minerador. Todas as rotas exportadas passam por verificação
independente de cobertura, prêmio, obrigatórios, ciclo simples e custo.

## Limites e verificações

**Correção isolada do decodificador:** `build.py --initialize-decoder-index`
cria `build_decoder_initialized/` e inicializa `cc_atual=0`, usado anteriormente
sem inicialização em componentes de uma única aresta. `verify_decoder.py`
reproduz esse caso mínimo. O método tem identificador próprio nos JSONs; não
misturar seus resultados com os do legado intacto. Três casos que falhavam no
legado concluíram e foram validados após a correção. Outros defeitos ainda
podem existir. Ver registros em `diagnostics/`.

- `fpmax_verification.json`: maximalidade, suporte e completude dos padrões em
  um banco pequeno para seis limiares, por enumeração independente. Não certifica
  o algoritmo de seleção dos dez padrões no wrapper.
- `checkpoint_verification.json`: a instrumentação preservou rota e custo para
  uma execução de mesma semente; um timeout artificial reteve rota viável.
- Os checkpoints são exportados nas iterações 1, 5, 10, 20, 35 e 70. Em timeout,
  `completed=false`, `timed_out=true` e `iterations` é a iteração do último
  checkpoint disponível. Não apresentar seu custo como melhor custo exato no
  instante do timeout: pode ter havido melhora após o checkpoint.
- A referência de custo em artigo anterior deve distinguir ótimo provado de
  melhor valor conhecido; igualdade numérica sozinha não é certificado.
- A licença local do Gurobi consultada expirou em fevereiro de 2026. Os novos
  reparos usam SciPy/HiGHS ou highspy, explicitamente identificados por resultado.

## Coleta por geometria

`experiments/ml_lns/strong_states.py` congela os papéis de geometria antes da
coleta. A seleção inicial foi a primeira variante lexicográfica. A opção
`--other-variants` coleta todas as demais variantes das mesmas geometrias sem
alterar os papéis. Não escolher variantes por seus resultados. Validação e teste
permanecem reservados durante o desenvolvimento.
