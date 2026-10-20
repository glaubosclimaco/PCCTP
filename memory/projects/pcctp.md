# Pesquisa PCCTP

## Decisões do usuário

Em 18 de outubro de 2026, o usuário escolheu continuar com aprendizado de
vizinhanças estruturadas por cobertura e prêmio, auxiliado por programação
inteira, e pediu atualizações dos resultados. A investigação atual permanece
no PCCTP; resultados preliminares não equivalem a contribuição inédita ou a
evidência de desempenho de publicação.

## Baseline legado recuperado

O usuário forneceu `mdm-grasp-genius-vnd/` em 18 de outubro de 2026, contendo o
MDM original e fontes do FPMax em `Debug/FPMax/`. Priorizar essa versão como
baseline, com compilação e instrumentação isoladas em `experiments/original_mdm/`.
A reconstrução anterior em `experiments/legacy_grasp/` é diagnóstica, não deve
ser apresentada como o método publicado nem misturada com estados do original.

Na retomada às 5h11 foi localizado um índice não inicializado (`cc_atual`) no
decodificador de padrões de uma aresta. Correção mínima isolada em
`build_decoder_initialized/`, método `original_MDM_legacy_decoder_index_initialized`.
O teste `verify_decoder.py` passou e três casos que falhavam concluíram com rotas
validadas. Preservar o legado intacto; não misturar resultados originais e
corrigidos. Estado e comandos atualizados em `research_notes/CONTINUATION.md`.

## Entregável final e critérios de encerramento

Em 18 de outubro de 2026, o usuário definiu como objetivo um artigo completo
em LaTeX, aproximadamente 20–30 páginas, pronto para submissão a periódico de
ponta, em pasta dedicada à submissão. Pode-se usar template Elsevier ou Springer
Nature. Pasta criada: `submission/`; `main.tex` é apenas um esqueleto inicial,
não um artigo concluído. O template editorial será aplicado durante a redação.
Incluir formulação matemática, tabelas, figuras, gráficos comparativos,
testes estatísticos, diagramas e casos quando forem necessários à argumentação.
Antes do encerramento, revisar integralmente o texto e sua consistência com os
dados; conferir contas, cálculos e experimentos; verificar existência, metadados
e adequação das citações e referências; compilar e inspecionar o PDF. Não inventar
resultados para preencher o artigo nem tratar publicação ou novidade como garantidas.

## Autonomia e isolamento Git

O usuário autorizou avançar autonomamente se estiver ausente ou se uma decisão
ficar sem resposta por cerca de 10 minutos. Escolher o próximo passo razoável
para alcançar a contribuição científica e o artigo revisado, documentando
suposições e decisões. Esse prazo não constitui aprovação para ações que
exijam autorização explícita por regras superiores. Não é necessário esperar
10 minutos para tarefas já autorizadas que podem continuar independentemente.
Reforçou depois: não esperar aprovação e seguir o plano até o final. Continuar
testando e aprimorando; se a linha atual não funcionar, investigar outras formas
de ML. A linha atual já usa ML para selecionar vizinhanças. Preservar controles,
resultados negativos e testes reservados; satisfação não é um limiar de p-valor.

Todas as alterações desta pesquisa devem permanecer na branch
`codex/pcctp-learning-submission`, criada em 18 de outubro de 2026 a partir de
`main`. A modificação preexistente em `include/Utils.h` pertence ao usuário;
preservá-la e não incluí-la em commits da pesquisa sem revisão específica.

## Extensão futura ao CTP

Instrução explícita do usuário, em 18 de outubro de 2026:

> Depois de termos resultados consolidados no PCCTP, aplicar também ao CTP
> (Covering Tour Problem) o que funcionar. Por enquanto, guardar a ideia e não
> implementar agora.

Essa é uma intenção futura condicionada à consolidação do PCCTP. Não iniciar
implementação, experimentos, geração de instâncias ou pesquisa específica de CTP
durante a etapa atual. Quando a condição estiver satisfeita, retomar a decisão
e avaliar a transferência dos mecanismos; não presumir que um efeito ligado à
quota de prêmio exista no CTP da mesma forma.

## Esclarecimento sobre contribuição científica

Em 18 de outubro de 2026, o usuário reforçou que aceita aplicação, adaptação
e combinação de métodos publicados para melhorar a resolução do PCCTP; não
exige técnica inédita. Priorizar evidência de desempenho e análise científica,
sem transformar novidade de faceta em requisito. ML é opcional. Preservar
resultados negativos e avaliações reservadas. Essa orientação atualiza a
prioridade operacional, sem dispensar baselines fortes ou rigor experimental.


## Prioridade reafirmada em19/10: avaliação experimental antes da decisão de artigo

Depois de discutir eventual troca de problema, o usuário pediu explicitamente
experimentos reais no PCCTP para decidir se vale escrever um artigo. Priorizar
fechamento/auditoria do estudo300s81runs e relatório honesto de ganhos/perdas,
curvas de qualidade e exposição efetiva à mineração. Não iniciar busca de outro
problema ou priorizar páginas do rascunho em detrimento desse diagnóstico.
Plano concreto em research_notes/pcctp_decisive_evaluation.md; respeitar gates
originais e teste reservado fechado. Não prometer publicação a partir de treino.


## Alinhamento com os artigos, após questionamento direto sobre tempo

Em19/10 o usuário pediu conferir se os experimentos eram os mesmos dosartigos.
A conferência confirmou70iterações/10sementes no protocolo heurístico, não
300segundos/3sementes. Os lotescurtos são complementares e não sustentam
superioridade sobreMDMcompleto. Antes da decisão de artigo, reproduzir referência
MDM70completa, sem confundir censura com término, e preparar comparação justa
comdezsementes. Calibração9TRAIN/seed80 está nafila após81runs300s;
éapenasdimensionamento/validação, não réplica final. Nãoabrirteste12/CTP.
