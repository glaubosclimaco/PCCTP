# Persistência e confronto descritivo com MDM publicado

Custos mínimos de três sementes do híbrido. Valores publicados extraídos da tabela tab:mdm-grasp do manuscrito fornecido; não são necessariamente ótimos. Orçamentos/hardware diferem, portanto não há inferência de velocidade frente ao artigo.

| Instância | Melhor MDM publicado | Híbrido10s | Híbrido30s | Gap30s % |
|---|---:|---:|---:|---:|
| g_a280_VT56_T168_W56_50 | 2001 | 2003 | 2073 | 3.598 |
| g_gil262_VT51_T158_W53_75 | 1987 | 2010 | 2010 | 1.158 |
| g_gr229_VT45_T138_W46_50 | 96707 | 99311 | 98217 | 1.561 |
| g_kroA200_VT40_T120_W40_25 | 23495 | 23495 | 23495 | 0.000 |
| g_pr226_VT44_T136_W46_50 | 67286 | 67852 | 67852 | 0.841 |
| g_rd400_VT80_T240_W80_25 | 10529 | 10801 | 10798 | 2.555 |
| g_si175_VT35_T105_W35_50 | 6987 | 7031 | 7006 | 0.272 |
| g_ts225_VT45_T135_W45_50 | 108805 | 108989 | 108989 | 0.169 |
| g_tsp225_VT45_T45_W135_75 | 2447 | 2450 | 2447 | 0.000 |

Melhorias estritas sobre a tabela publicada: 0; igualdades: 2.
O pós-processamento isolado melhora a baseline de curto prazo em parte dos casos, mas esses dados não sustentam novos recordes ou uma contribuição principal suficiente. Próxima hipótese: alternar mudanças de conjunto visitado com ordenação, apenas em treino/desenvolvimento, com controles de busca sem LKH e sem mudança de suporte. O teste reservado continua fechado após falha do gate de validação.
