# Repetição integral após suspensão do computador

Em 20/10/2026 às 08h30 de Fortaleza, a inspeção identificou suspensão por
fechamento da tampa às 01h00min27s e repetidos ciclos de manutenção até
08h25min45s. O log do macOS está preservado em sleep_events_20261020.log.
O lote original foi interrompido com 44 registros completos preservados;
13 excederam 330 segundos de relógio. A auditoria de rotas e timestamps
permanece distinta da validade da comparação temporal. Esses dados não
fundamentam um teste de desempenho sob condições controladas.

Repetir todos os 81 ensaios em ruin_total300_awake_train, mantendo instâncias,
sementes 70–72, ordem dos métodos, binários, parâmetros e orçamento de 300 s.
Não selecionar repetições a partir do custo nem combinar as primeiras
geometrias do lote anterior com as novas. Preservar os dados anteriores.

Executar com caffeinate -i temporário, sem alterar preferências do sistema.
Isso reduz suspensão por inatividade e não impede suspensão por tampa fechada.
O computador precisa permanecer acordado. Um monitor amostra relógio de parede
e monotônico a cada segundo durante cada bloco de três métodos. Um intervalo
maior que 5 s, diferença entre relógios maior que 3 s ou evento de suspensão
do pmset no período invalida o bloco e interrompe a fila. Os limiares são
critérios operacionais conservadores definidos antes desta repetição;
não demonstram igualdade de carga, temperatura ou CPU.

Cada bloco grava hashes dos três resultados e seu diagnóstico ambiental.
Retomadas aceitam apenas blocos completos com diagnóstico PASS. Um bloco
parcial ou reprovado exige preservação e revisão explícita antes de repetir;
o programa não o transforma silenciosamente em observação válida.

Os observadores anteriores serão reutilizados com substituição explícita
dos nomes da pasta e relatórios, sem editar suas fontes congeladas. A
elegibilidade ambiental exige adicionalmente environment_audit.json PASS
com 27 blocos, além da auditoria de rotas, custos e deadlines.

O critério de treino permanece: ganho médio ≥0,2%, pelo menos seis vitórias
e no máximo uma perda entre nove geometrias, sem falhas ou ausências. As duas
variantes serão reportadas. Nenhum teste reservado será aberto. A referência
completa de 70 iterações ainda é necessária antes de julgar o artigo.
