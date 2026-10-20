# Conferência estática do MDM recuperado e do artigo fornecido

A conferência verificou os16arquivos instrumentados contra os hashes do
manifesto de compilação. Não recompilou nem modificou o otimizador em execução.
O certificado é experiments/original_mdm/published_parameter_inventory.json.

| Parâmetro | Código recuperado | Evidência/limitação |
|---|---|---|
| maxIt |70|principal.cpp; critério também descrito nos artigos |
| alpha |0.6 fixo|não é a política alpha reativa do R-GRASP2021 |
| capacidade elite |10|grasp.cpp |
| máximo de padrões |10|grasp.cpp |
| suporte |2 absoluto|mineracao.cpp311/320;20% apenas para banco cheio10 |
| estabilidade |ceil(70*5/100)=4iterações|grasp.cpp; preservar no protocolo |
| política |original|retorna antes de features/reordenação aprendida |
| budget da referência70 |0|desabilita saída temporal; permanece numIt70 |

O manuscrito MDM descreve suporte20%, mas os trechos conferidos não fixam todos
os valores numéricos alpha/capacidade/estabilidade. Não inferir identidade do
executável publicado somente por usar o código recuperado. A correção ocupação
elite e a inicialização do decodificador permanecem diferenças explícitas.
Não passar a recalcular suporte20% sobre ocupação parcial silenciosamente:
isso mudaria o algoritmo/configuração de todos os lotes já congelados.

A política original retorna de schedule_order antes da extração de descritores
e ordenação. O wrapper dos movimentos também não calcula recompensa/tempo
adaptativo nessa política. Isso foi verificado no fonte compilado; não é uma
medição de overhead nulo. Hooks de progresso/checkpoint ainda são instrumentação.

Hardware local observado durante o estudo: AppleM2,arm64,16GiB,8núcleos
físicos/lógicos; compilador do manifesto GCC15.1.0. O manuscrito MDM informa
Inteli7-6700/3.40GHz/16GB e g++5.4.0. Não normalizar tempos por GHz nem
comparar diretamente os segundos como medida de aceleração. Snapshot com
momento UTC em experiments/original_mdm/research_machine_snapshot.json.
