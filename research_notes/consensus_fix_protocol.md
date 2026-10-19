# Pool-guided visit fix-and-optimize, before collection

Nine hardTRAIN, same original edge_pool starts and input routes. Free12 optional
visit variables (or all if fewer); fix others to incumbent. Consensus arm chooses
highest frequency*(1-frequency) over pool, ties by global vertex index. Random
arm samples same12, seed42. ALL route edges free in both; unlike earlier fragment
repair, mandatory-to-mandatory edges and all routing changes remain possible.
Standard SCF, warmstart,60s including setup per arm. Alternate arm order.
Gate: consensus improves>=3/9 entries and mean>=0.2%; beats random>=5/9 and
fullSCF>=3/9. No held-out data or increasing k after seeing results. Offline
input generation excluded. Fixings heuristic, never global reductions; bounds
local. Raw scalar audit must reconstruct releases, fixed visits and routes.
Two small neighborhoods checked against exhaustive enumeration; invalid
incumbent-incompatible fixing rejected. No novelty claim for fix-and-optimize.
