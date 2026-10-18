# Piloto de tempo total: MDM e MDM+LKH

Congelado18/10/2026 após screen7/9,antes do piloto pareado. TODOS os9TRAIN
maiores,seeds42,43,44. ControleMDM10s;híbridoMDM9s seguido por LKH em no
máximo1s restante,atédeadlineTOTAL10s. Inicialização/leitura/exportações
cobradas. SeMDMretornar tarde,não executar LKH semtempo. IncumbenteMDM
emitido a tempo permanece fallback;overrun e hardtimeout reportados.

LKH3.0.6,mesmas opções do screen,seed42fixa no backend para isolar variação
doMDM,RUNS1,MAX_TRIALS100. Sem tuning9/1. Ordem dos braços randomizada de
forma reprodutível por instância+seed. Arquivos privados por braço/caso/seed.

Gate para avaliação maior: ganho relativo médio de custo por geometria>=0.5%,
melhor média em>=6/9geometrias e pior em<=1/9,sem perda de factibilidade.
Não tratar27runs como27geometrias independentes;análise pareada por9geometrias.
Relatar perdas/empates e intervalos exploratórios. Se falhar,não ajustarfração
de tempo nem selecionargeometrias deste lote. Teste reservado fechado.
Este piloto não autoriza alegar novo recorde ou vantagem estatística confirmada.
