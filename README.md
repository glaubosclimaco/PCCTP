# PCCTP
A grasp heuristic for the Prize-Collecting Covering Tour Problem

Créditos: Rogégio da Silva

The definition of the PCCTP can be given as follows. Consider an undirected graph
G = (N, E), with the set of vertices N = V ∪ W and V = R ∪ T. Let D be a coverage
distance, ce be a cost associated to each e ∈ E and pi be a prise associated to each
vertex i ∈ V . Let T be a subset containing the vertices that must be visited and R a
subset of vertices that are optional, and therefore may or may not be visited. Finally,
let W be a subset containing the vertices that must be covered by another vertex from
V , i.e., a vertex i ∈ V covers a vertex j ∈ W if cij ≤ D. The goal of the PCCTP is to
find a simple minimum cost cycle that visits all vertices in T, covers all vertices of W,
and collects at least a minimum prise (PRISE)

