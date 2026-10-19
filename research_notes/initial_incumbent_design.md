# Inicialização e ausência de incumbentes em prazos curtos

Diagnóstico de código, 19/10/2026. No build temporizado do MDM, a primeira
chamada research_progress da construção ocorre depois de construir a solução
inicial (grasp.cpp, perto da linha3491). A rotina construtiva pode ultrapassar
o prazo antes de registrar qualquer solução. Assim, reservar10s ao MDM dentro
de um orçamento15s não assegura que o pós-processamento receba um incumbente.
Não se pode interpretar falta de incumbente como custo zero, excluir o caso ou
reaproveitar uma rota de uma execução anterior sem cobrar sua obtenção.

Uma próxima investigação possível é fornecer, a TODOS os braços, um incumbente
externo factível de obtenção rápida, debitado do orçamento total. Selecionar todos
os vértices visitáveis fornece cobertura se cada demanda possui um vizinho e
prêmio suficiente se a soma total atende à quota; ambos devem ser verificados.
Um ciclo desses vértices é factível com pelo menos três vértices e inclui todos
os obrigatórios. Isso não pressupõe desigualdade triangular nem fornece qualidade.
Melhorias por remoção factível podem ser interrompidas a qualquer momento,
retendo o último ciclo validado. A construção inteira, leitura e validação
precisam entrar no tempo total. Não executar esse novo estudo antes de auditar
por completo o lote15s atual; resultados ausentes devem permanecer preservados.

O fallback comum resolveria disponibilidade, mas não demonstraria contribuição
científica por si. Seria preciso comparar qualidade e tempo com esse mesmo
fallback em MDM, LKHfixo, suporte2opt e suporteLKH, e relatar a proporção de
soluções que ainda usam apenas o fallback. Nenhuma escolha de fallback pode usar
alvos ou soluções do conjunto de teste. Não abrir teste reservado nesta etapa.
