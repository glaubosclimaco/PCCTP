# Conditional TRAIN screen: joint insertion versus greedy repair

Frozen before any benchmark outcome for joint insertion is collected. Run only
after the current validation batch is complete and its audit/summary are recorded.
If that candidate passes validation, prioritize its planned reserved evaluation;
this alternative is not a reason to tune it. Current partial results already make
its three-win gate unreachable; nevertheless finish all scheduled runs.

Question: does optimizing visits and insertion gaps jointly yield better repaired
routes than the current greedy rule, after the same destruction and LKH allowance?
This is a component screen, not a new end-to-end method, learned policy, novelty
claim or equal-total-time comparison including MDM generation.

Use only the nine hardTRAIN geometries, the three saved MDM10 pre_routes from the
audited corrected absolute batch (seeds42-44), and four deterministic removals
per start. Never use the improved final hybrid routes as starts. Remove uniformly
min(6,optional_count,route_length-3) optional vertices, preserving mandatory visits.
Use the same destroyed retained cycle for all three arms, including cases already
feasible after removal. Original incumbent is the common feasible fallback.

Arms: (a) existing randomized greedy repair; (b) joint insertion MIP with the three
cheapest retained gaps per absent vertex; (c) the same MIP with three uniformly
sampled distinct gaps per absent vertex. All absent vertices are candidates in
both MIP arms and both have the same variable count. No matching augmentation or
post-hoc candidate expansion. The restricted model may be infeasible even when
the original problem is feasible. Never classify original infeasibility from it.

Each arm receives a 1.0-second conditional allowance, with slot selection/model
construction/repair and polishing charged. MIP receives at most0.5seconds including
its model construction; greedy runs until the common deadline. Each timely feasible
repair receives at most0.2seconds fixed-support LKH, further limited by the remaining
allowance. These are maximum allowances, not a requirement to waste idle time;
report actual times and attempt counts. Keep the original incumbent if repair is
infeasible, late or worse; record that disposition separately. Timely feasible
repair improvements can be retained even if subsequent polishing returns late.
Use original scalar cost/feasibility checks and absolute CLOCK_MONOTONIC timestamps.
Data loading and common state/destruction generation are outside this conditional
screen and explicitly excluded. No parallel optimizers. Arm order is deterministically
shuffled per matched case.324scheduled arm executions, no optional stopping.

Primary screen: per geometry average all12costs per arm, compare nearest-gap MIP
with greedy by100*(mean_greedy-mean_MIP)/mean_greedy. Require equal-geometry mean
>=0.2% and at least3/9positive means. Specific slot-selection component additionally
requires at least3/9positive geometry means versus random slots. All324scheduled
records and route/source/deadline audits must be present. Report reconstruction
success, restricted infeasibility, fallback count, timing and negative outcomes.
No p-value or improvement claim from this exploratory TRAIN gate. If it passes,
investigate integration with a fresh equal-total-time TRAIN experiment; do not
open validation/test or assume repeated MIP calls pay for their overhead.

Preserve hashes of saved starts, source files, LKH and HiGHS version, all seeds and
removed vertex lists. Technical interruptions retain their artifacts/reasons;
completed unfavorable outcomes are never rerun selectively. No parameter changes
after collection; a failed gate is a preserved negative experiment.
