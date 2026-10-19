# Ruin/recreate with coverage-prize repair: before collection

Nine hardTRAIN, original edge_pool starts,60s extra/arm including read/setup.
Remove6 optional visits (cap preserves3vertices). Random arm uniform seed42;
related arm nearest to random removed seed by matrix cost, ties by index.
Repair temporary infeasibility using insertion cost divided by uncovered demand
count plus capped prize deficit normalized by mean visit prize; random score
factor uniform0.8..1.2. Stop once coverage/quota satisfied. LKH up to0.2s per
repaired support. Accept current within1% of best, reset to best every20attempts;
best never worsens. Raw feasibility/cost audit and native acceptance deadline.

Compared with earlier diversification: starts are strong, removals may violate
coverage/prize and are explicitly repaired, rather than only feasible1for1swaps.
No adaptive operator weighting, so do not call this ALNS implementation.
Background family: Ropke/Pisinger2006,DOI10.1287/trsc.1050.0135, primary
INFORMS/DTU metadata checked; no claim this reproduces their VRP algorithm.
Gate: related improves>=3/9,mean>=0.2%,beats random>=5 and fullSCF60>=3.
Random arm reported independently. No retuning after results, no test unlock.
Offline input-generation cost excluded; not an end-to-end claim or novelty proof.
