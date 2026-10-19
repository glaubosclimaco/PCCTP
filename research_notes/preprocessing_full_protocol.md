# Frozen preprocessing comparison,19October2026

Nine hardTRAIN, same starts as completed full_graph60,60seconds each including
reading, reducing and model setup. Standard HiGHS presolve remains enabled.
Apply only structural.py constraint reductions; no route or edge elimination.
Record removed coverage rows, newly mandatory vertices, preprocessing seconds
and source hash. Independently validate returned routes against ORIGINAL raw
instance, not only reduced model. External native-clock late discard unchanged.

Potential gate: at least three lower accepted costs than original full60 and
mean relative gain>=0.2%, no feasibility loss. Runtime/node comparisons remain
descriptive because these are sequential machine runs; no speedup claim from
row counts. If gate fails preserve negative; do not retune on held-out sets.
