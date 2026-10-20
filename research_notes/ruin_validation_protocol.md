# Frozen validation screen for randomized removal/repair hybrid

Written before this candidate has any new validation output. Prerequisites:
corrected hardTRAIN gate passed; complete same-start component ablation passed.
Seven validation geometries from the original split, original one-instance
selection unchanged. Seeds45-49, both arms,60seconds total, same common fallback,
MDM60 control versus MDM10+random ruin/recreate remainder. Identical search
parameters, binaries and absolute event eligibility from the training batch.
Randomize arm order per instance/seed; record hashes before data collection.
No alternate removals, longer LKH call, adaptive cutoff or geometry-specific rule.

Primary screen outcome: relative reduction of five-seed mean cost, averaged
with equal geometry weight. Screen requires mean>=0.2%, at least3/7 positive
geometry means, at most1/7 negative, no feasibility regression or unaudited
missing runs. This is a new candidate-specific gate, not a retroactive change
to previous failed9+1LKH or related-vs-random gates. No inferential claim from
this screen. Preserve ties, losses and fallback cases; no selective restarts.

These validation geometries were used in the earlier different9+1LKH study;
therefore this is limited validation screening, not pristine confirmation.
Do not adjust the candidate after inspecting these outcomes. If the screen
fails, do not open the12test geometries for this configuration. If it passes,
prepare a fully frozen final test protocol before reading any test data.
Test outcomes remain unavailable to this runner; it filters roles=='validation'.

Publication scope remains problem-specific. The60second baseline rarely
activated mining in training, so report that exposure and later evaluate
longer horizons or published stopping criteria before broad claims. Novelty
cannot be asserted for GRASP+MIP, LKH or destruction/repair alone. Component
and generalization evidence are required. Algorithm changes following this
screen must be investigated on training, not tuned on these validation results.
