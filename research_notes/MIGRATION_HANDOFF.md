# Retomada em outra maquina

## Estado de entrega

Execucao local encerrada a pedido do usuario. Automacao recorrente removida.
Nao ha retomada automatica. Branch: codex/pcctp-learning-submission.
Leia diario.md, research_notes/CONTINUATION.md, CLAUDE.md,
memory/projects/pcctp.md e research_notes/EXECUTION_PLAN.md nesta ordem.

Coleta paired70_rebuilt_train: 2/90 pares completos, auditoria PASS,
236 verificacoes, zero geometrias completas. O terceiro par a280/82 foi
interrompido por solicitacao do usuario e arquivado integralmente em
experiments/lkh/paired70_rebuilt_interrupted/20261022T150221Z.
Nao reutilizar sua referencia ou candidata parcial como par completo.
Cinco pares de paired70_train pertencem ao backend anterior e ficam separados.

## Primeiro acesso

    git clone --branch codex/pcctp-learning-submission https://github.com/glaubosclimaco/PCCTP.git
    cd PCCTP
    python3 experiments/lkh/audit_paired70_rebuilt.py

O auditor e somente leitura dos resultados, com gravacao de relatorios;
nao inicia otimizacao. Consulte experiments/original_mdm/README.md e build.py
para compilacao. As fontes e executaveis de build_timed_elite_bounded foram
preservados. Backups dos quatro executaveis usados estao em
experiments/migration/runtime_macos_arm64/, com SHA256 no manifest.json.
Sao macOS arm64; nao executam diretamente em Linux/x86. LKH tem fonte/revisao
registrada em experiments/lkh/rebuilt_backend/provenance.json. Nao confundir
reconstrucao de codigo com equivalencia de binario ou de desempenho.

## Antes de iniciar otimizacao

O runner congelado contem caminhos absolutos do Mac original e referencias a
/private/tmp; nao executa automaticamente em outra maquina. NAO editar seus
hashes para forcar retomada. Criar runner/protocolo/diretorio novos para o novo
host, declarando CPU/OS/compilador/binarios, caminhos e dependencias.
Refazer verificacoes sinteticas, locks reais e calibracao no novo host.
Reiniciar todos90pares no novo host, mantendo os pares deste Mac separados;
medir MDM70 fresco e cobrar exatamente esse tempo de cada candidata no mesmo
host. Nao combinar tempos de equipamentos diferentes ou normalizar por GHz.
Preservar9TRAIN,seeds81-90,ordem/gates e todas falhas/censuras. O novo estudo
continua exploratorio, pois resultados TRAIN ja foram vistos.
Nao abrir12geometrias reservadas, implementarCTP ou duplicar otimizadores.
ML e opcional. Cache nao integra candidatos; gate condicional completo FAIL.

## Evidencias e artigo

awake300:81execucoes validas;H60+0,979441%,H10-1,279343%,TRAIN apenas.
Calibracao MDM70:9/9 completos,uma semente apenas,nao resultado final.
Cache30:27pares/54fases,501verificacoes,ganho0,037711%,gateFAIL.
Detalhes e limites nos relatorios e diario. Artigo submission/main.tex/PDF
parcial20paginas,abstract/conclusao/evidencias finais pendentes;nao pronto.
Manifesto historico do artigo espera LKH original perdido; nao reescrever
para esconder divergencia do novo binario.

include/Utils.h contem alteracao local preexistente do usuario: nao foi
incluida neste envio, conforme instrucao explicita. Nao e dependencia do
baseline isolado congelado. Locks *.queue.lock nao sao versionados.
