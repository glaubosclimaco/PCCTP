# Limite de raiz com o corte de partição

Registrado antes de ler os limites deste lote. As linhas de aprendizado, de escolha de padrão e de reabertura não superam o MDM na qualidade. O corte de partição já tem prova e testemunha verificada. Falta saber se ele fortalece o limite linear nas instâncias reais.

O modelo é a relaxação linear da formulação com fluxo usada no reparo, sem variáveis inteiras. A partição sai das componentes em que a aresta fracionária vale pelo menos 0,5. A componente da raiz é o bloco da raiz. Se restarem mais de 8 blocos, permanecem os 8 de menor índice mínimo e os outros entram no bloco da raiz. O posto é o mínimo exato de blocos fora da raiz numa seleção que cumpre obrigatórios, cobertura e quota, por enumeração de no máximo 256 subconjuntos. Se esse posto é pelo menos 1 e a soma fracionária das arestas entre blocos fica abaixo do posto mais 1, o corte entra e a relaxação é reotimizada. No máximo 5 rodadas. Não há outro limiar.

Cada corte emitido tem de ser satisfeito pela rota final do MDM quando esse custo coincide com um ótimo já provado. Se uma rota ótima violar o corte, o separador está errado e o lote não serve.

O corte é computacionalmente útil neste conjunto somente se o limite final for estritamente maior em pelo menos 5 das 9 geometrias de desenvolvimento. Caso contrário, a família continua válida na testemunha, mas este separador não justifica um artigo computacional. Validação e teste ficam fechados. O limiar não será mudado depois do resultado.

## Prova inteira, regra registrada antes dessa segunda medição

O limite de raiz pode subir sem acelerar a prova. A segunda medição usa o mesmo separador, sem mudar limiar, blocos ou rodadas. Para cada geometria de desenvolvimento há duas provas inteiras de 30 s: o modelo puro e o modelo com os cortes acumulados na raiz. O tempo do segundo inclui a separação. Uma prova conta como mais rápida só quando termina e gasta no máximo 80% do tempo da outra. A família segue para um texto computacional somente se provar um conjunto estritamente maior, ou se, entre as que ambas provam, for mais rápida em mais geometrias do que o contrário.

## Resultado do limite

O separador não emitiu corte inválido para a rota viável do MDM. O limite subiu em 7 das 9 geometrias, então a regra de utilidade do limite passa. A folga relativa ao ótimo provado, antes e depois:

- att48: 6,24% para 5,26%, cinco cortes, ótimo 2650
- berlin52: 17,33% para 15,52%, cinco cortes, ótimo 5186
- brazil58: 36,75% para 25,19%, cinco cortes, ótimo 20466
- eil76: 3,15% sem corte, ótimo 332
- gr48: 6,56% para 4,22%, dois cortes, ótimo 3571
- hk48: 6,82% para 6,38%, um corte, ótimo 7593
- st70: 20,44% para 19,24%, quatro cortes, ótimo 419

eil51 não recebeu corte. pr76 subiu o limite, sem ótimo provado nesta etapa. O ganho grande está em brazil58. Nos outros, a subida é pequena. O teto de cinco rodadas não será aumentado depois desse resultado.

## Resultado da prova

As duas versões provaram as mesmas oito geometrias. pr76 ficou sem prova nos 30 s, nas duas. Nenhuma versão com cortes foi mais rápida pela regra dos 80%. A versão pura foi mais rápida em berlin52 (12,8 s contra 18,4 s), gr48 (1,1 s contra 3,8 s) e hk48 (3,0 s contra 4,3 s). A regra de continuidade da prova falha.

Decisão: não escrever artigo de aceleração com este separador. O efeito no limite de raiz fica registrado, em especial brazil58, e não será vendido como prova mais rápida. A novidade frente a cortes de GTSP e de covering salesman continua em aberto. Validação e teste não foram abertos.
