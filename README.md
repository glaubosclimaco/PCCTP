# PCCTP

A grasp heuristic for the Prize-Collecting Covering Tour Problem.

Créditos: Francisco Glaubod and Rogégio da Silva.

The definition of the PCCTP can be given as follows. Consider an undirected graph
G = (N, E), with the set of vertices N = V ∪ W and V = R ∪ T. Let D be a coverage
distance, ce be a cost associated to each e ∈ E and pi be a prise associated to each
vertex i ∈ V . Let T be a subset containing the vertices that must be visited and R a
subset of vertices that are optional, and therefore may or may not be visited. Finally,
let W be a subset containing the vertices that must be covered by another vertex from
V , i.e., a vertex i ∈ V covers a vertex j ∈ W if cij ≤ D. The goal of the PCCTP is to
find a simple minimum cost cycle that visits all vertices in T, covers all vertices of W,
and collects at least a minimum prise (PRISE).

---

## Como compilar

**Dependências:** Gurobi (solver MIP) e LEMON (grafos). No Debian/Ubuntu instale LEMON com:
```bash
sudo apt-get install liblemon-dev
```
O Gurobi deve estar instalado (ex.: em `/opt/gurobi*/linux64`). Se estiver noutro sítio, defina `GUROBI_HOME` antes de compilar.

Na raiz do projeto:
```bash
make clean
make
```
O executável gerado é `pcctp`.

---

## Como executar

Uso geral:
```bash
./pcctp <instancia> <seed> <aplica_cortes> <remove_E0> <gamma> <aplica_regras> <ordem> <perfix> [modo]
```

- **modo** (opcional): `0` = GRASP + solver MIP (padrão); `1` = só heurística (GRASP).

**Exemplos:**

| Modo        | Comando |
|------------|--------|
| Só solver (MIP) | `./pcctp data/instancia.txt 42 1 0 0.5 0 1 0.5` |
| Só heurística   | `./pcctp data/instancia.txt 42 1 1 0.5 0 1 0.5 1` |
| Híbrido (GRASP + MIP) | `./pcctp data/instancia.txt 42 1 1 0.5 0 1 0.5` |
| Regras de redução    | `./pcctp data/instancia.txt 42 1 1 0.5 1 1 0.5` |

Os resultados são gravados em `resultados/<metodo>/` com data e hora no nome do ficheiro. Nos modos **solver** e **hibrido**, o ficheiro contém todo o log do Gurobi e termina com a linha **MIP** (instância, ótimo, custo, tempo).

Mais detalhes em [BUILD.md](BUILD.md).
