# Analysis safeguards before any reserved evaluation

All current hardTRAIN batches are exploratory. Repeated attempts, failed gates
and technical reruns must be disclosed; nominal training p-values would not
establish a confirmed contribution. Keep the unchanged gates for each batch.
Do not drop rd400, failed seeds, fallback starts or overrun records.

Primary geometry outcome for the current two-arm screen is the relative
reduction of seed-mean cost:100*(mean_MDM-mean_hybrid)/mean_MDM. Average this
over complete geometries with equal geometry weights. Retain all three seed
pairs in the artifact. Do not replace this by mean of per-seed percentages,
median or a weighted average after seeing an outlier. Secondary metrics may
be reported explicitly, never substituted for the primary gate.

No inferential analysis until all pairs pass route/deadline/source audits.
Future inferential unit is geometry, not a seed or quota variant. A publication
comparison should report paired effect sizes and uncertainty, feasibility and
cost-time tradeoffs, plus exact or permutation paired tests when assumptions
are defensible. If bootstrap intervals are used, resample whole geometries;
state the small-sample and benchmark-selection limitations. Preserve zero
and negative effects. Correct multiple planned method comparisons, and label
all post-hoc subgroups as exploratory. Seeds do not multiply independent cases.

The original12test geometries remain unopened. Validation7 was previously used
for another frozen LKH variant which failed; it is not a pristine independent
confirmation sample. No new parameter tuning there. Any new reserved evaluation
requires a frozen algorithm, source/binary hashes, primary metric, budget grid,
seed list, stopping rules, component controls and explicit handling of failures.
Do not unlock test solely because exploratory sensitivity reconstruction looks
positive. The current corrected end-to-end batch and component ablation must
first be complete, audited and assessed against their existing criteria.
