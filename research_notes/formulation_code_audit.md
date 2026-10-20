# Correspondência entre formulação, implementação e certificados

Auditoria de leitura em20/10/2026. Os hashes dos arquivos inspecionados estão
em formulation_code_audit.json. Não houve nova execução de solver nem alteração
do manuscrito/PDF ou dos experimentos ativos.

1. **Modelo de referência.** Em problem_formulation.tex, grau2y, raiz obrigatória,
   arestas binárias e cortes enraizados implicam um único ciclo simples de pelo
   menos três vértices. Visitas obrigatórias, cobertura por conjuntos e prêmio
   apenas de vértices visitados coincidem com pilot.validate e raw_route_audit.
   O argumento não exige desigualdade triangular. A hipótese de raiz obrigatória
   deve permanecer explícita; não generalizar o modelo para T vazio sem adaptação.

2. **Representação implementada.** pilot.Repair usa fluxo de uma mercadoria
   para conectividade, e não materializa todos os cortes enraizados do texto.
   O fluxo exige uma unidade por vértice selecionado diferente da raiz e é
   limitado por(n−1)x nas arestas. Essa distinção é de representação; o texto
   deve continuar chamando a formulação de referência, sem afirmar que foi
   exatamente a matriz submetida a todos os solvers.

3. **Contração.** fragment_repair forma componentes das arestas fixadas a1,
   exclui arestas fixadas a0, mantém grau residual por extremidade original
   e soma o custo fixo como constante. Componentes não triviais são obrigatórios;
   vértices opcionais isolados podem ser selecionados quando liberados. A cobertura
   e o prêmio são agregados sobre os vértices originais do componente. Fluxo ou
   separação de cortes asseguram conectividade. Uma volta inteira retida não tem
   capacidade de grau para incorporar outro ciclo. O argumento do texto corresponde
   a essas decisões, sem constituir evidência de aceleração ou novidade.

4. **Escopo dos limites.** Uma solução viável de reparação é um limite superior
   do PCCTP, após validação original. Um limite inferior ou prova de ótimo da
   vizinhança restrita não prova limite inferior/ótimo global do PCCTP. O modelo
   separado de inserção permite no máximo um novo vértice por lacuna retida;
   não deve ser descrito como solução exata de toda reinserção possível.

5. **Certificados históricos.** fragment_verification.json registra37 casos,
   sendo15 por enumeração,20 comparações com SCF e2 liberando apenas ausentes.
   native_cut_verification.json registra60 verificações e **32 cortes**.
   Uma nota operacional antiga dizia35; usar o artefato correto se o número
   for mencionado. Esses certificados não contêm hashes das fontes. O hash
   registrado agora identifica o arquivo lido, mas não demonstra retroativamente
   que o teste histórico usou o código atual. Antes de apresentar esses testes
   como certificação da versão final, reexecutá-los de forma serial, com fontes
   congeladas e hashes gravados. Não lançar outro solver enquanto a coleta está ativa.

O manuscrito continua parcial. Esta auditoria verifica correspondência e limites
de interpretação; não acrescenta resultado computacional, validação estatística
ou contribuição inédita.
