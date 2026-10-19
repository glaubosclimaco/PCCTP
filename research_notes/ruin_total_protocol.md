# Equal total60s: original MDM versus MDM + randomized coverage/prize ruin repair

Before collection, nine hardTRAIN x seeds42/43/44, two arms, randomized order
per geometry/seed. Same feasible fallback construction, charged to60s. MDM arm
continues original corrected elite-bounded MDM up to60s; hybrid arm MDM10s then
random-removal ruin/recreate for remaining time. No offline route pool allowed
as initial solution. Best valid native/fallback incumbent used, preserving exact
original-data feasibility. Native external safety timeout90s; overrun cannot
contribute late solutions, predeadline exports/fallback preserved explicitly.
If initialization overruns, keep accepted fallback and record lost repair time.

Fix removal6, repair/acceptance/reset rules from prior protocol; choose RANDOM
because related/random tied all nine previous costs, no evidence of related
advantage. No adaptive selection/ML claim. Gate: mean geometry paired gain>=0.2%,
at least6/9 geometry mean wins, at most1/9 losses, no feasibility regression.
Aggregate seed pairs within geometry before comparisons; no treating27 as
independent geometries. Train screen only; no test or reused validation access.
Record full wall and CPU metadata, original raw/checkpoint exports and source
hashes. No broad ranking or significance before audit and complete results.
