# Prospective final test protocol — conditional, not yet launched

Recorded 19 October 2026 while validation is still running. This document does
not authorize bypassing the existing validation gate. No test instance data or
test outcomes were read to prepare it. A failed validation closes this candidate;
do not run this protocol for that candidate or adjust its parameters on validation.

## Entry and configuration freeze

Require complete PASS audits for corrected training, the conditional fixed-support
ablation, and the existing seven-geometry validation screen. Validation must meet
its unchanged mean >=0.2%, >=3 positive geometry means, <=1 negative geometry mean,
and no feasibility regression/missing runs. Check source and binary hashes against
the validation protocol before any test instance is loaded. An unexplained hash
change blocks launch. Analysis and runner code may differ; optimizer dependencies
must not. Store this document's hash with the launch manifest.

Use exactly the twelve role=='test' geometries and their previously selected
instances in the original split; no choosing difficult/easy instances after access.
Seeds 100 through 109 inclusive, two arms, 60 seconds total each: 240 runs.
Control is the corrected, elite-bounded original MDM; hybrid is the unchanged
MDM10 + randomized removal/greedy repair/LKH remainder. Common charged fallback,
absolute CLOCK_MONOTONIC incumbent eligibility, identical frozen parameters and
binaries. Preserve the decoder and elite-bound corrections as disclosed differences
from the legacy code. Randomize paired arm order deterministically from instance
name and seed, as in validation. Execute sequentially, with no other optimizer.
Initialization, file construction, selection and optimization time are charged.
Cleanup time can exceed 60 seconds; later incumbents are ineligible. Report this
distinction and measured process duration.

## Outcomes and unit of analysis

For geometry g compute M_g and H_g as the arithmetic means of all ten eligible
control and hybrid costs. Primary outcome is the equally weighted mean of
G_g=100*(M_g-H_g)/M_g over all twelve geometries. Positive favors hybrid.
Do not substitute a mean of seed percentages or weight by route cost/instance size.
Report all twelve paired means, G_g, median G_g, range and wins/ties/losses
(tie tolerance 1e-9 percentage points). Preserve all seed results and feasible
fallbacks. Seeds quantify within-instance variability; they are not 120 independent
geometries. The existing single-instance-per-geometry selection limits generality
across coverage/quota regimes and must be stated.

Secondary summaries: paired seed differences, feasibility counts, elapsed process
time, fallback use, mining exposure, accepted improvement counts and calls/timeouts.
They do not replace the primary outcome. No post-hoc subgroup is confirmatory.
Do not compare header targets as certified best-known solutions without checking
their provenance and comparable problem definition.

## Uncertainty and decision

Report a descriptive 95% percentile interval for the mean G_g by resampling whole
geometries (10,000 replicates, NumPy default_rng seed 20261019). No seed-level
resampling presented as independent-instance evidence. With only twelve selected
benchmark geometries this interval describes sensitivity to benchmark composition;
population coverage requires an exchangeable sampling model that this fixed suite
does not establish. Show the individual effects even if the interval is favorable.

As an explicitly secondary analysis report a two-sided exact binomial sign test
on nonzero geometry effects, null probability 1/2, with ties excluded and their
count reported. If no nonzero effect exists, p=1. Its independence/equal-sign
assumptions and limited power must be explicit. A p-value is not a publication
criterion. There is one planned method comparison at this horizon; subsequent
methods/horizons need their own prospective design and multiplicity treatment.

Practical performance claim requires complete auditable data, no feasibility
regression, primary mean >=0.2%, median >=0 and more positive than negative
geometry means. Also report whether the descriptive interval includes zero.
If these conditions fail, retain the results, do not retune on this test, and do
not present the candidate as improving the reserved suite. Passing is evidence
for a scoped 60-second performance claim, not automatic sufficiency for a journal
submission or proof of superiority to all published methods. In particular,
longer-horizon/mining-active and published stopping-rule comparisons remain needed.

## Failure handling and stopping

Audit original scalar cost, quota, coverage, mandatory vertices and absolute event
timestamps for every retained route. Missing runs do not disappear from the
denominator. An algorithmic failure remains a failure; if a feasible common
fallback was produced on time it is the eligible result. If no feasible result
exists, report failure counts and do not compute a complete cost-only primary
claim by dropping that pair.

Document infrastructure interruptions and preserve partial artifacts. Restart
only genuinely missing/technically interrupted executions, never a completed
unfavorable run. Keep both interrupted and replacement records, reasons and
provenance. A changed optimizer invalidates this frozen comparison; preserve the
original batch and design another study without calling reused test data untouched.
Stop after the prespecified 240 valid scheduled executions, not when significance
or a favorable mean is reached. No optional seed expansion based on outcomes.
