# Prospective design: comparison against completed70-iteration MDM

Prepared before seed80calibration outputs or any main-study outputs. This design
does not mean that main collection has started. Implement and verify its runner,
deadline handling and independent auditor before launch. Complete the seed80
reference calibration first. The current300-second queue remains unchanged.

## Scope and schedule

Use the same nine hardTRAIN instances, ten fresh seeds81through90 inclusive.
This is a problem-specific TRAIN comparison against a paper-compatible stopping
criterion, not a replication of all107published evaluation instances. No reserved
test instance is read. Do not claim the original article used these numeric seeds;
the passages examined specify ten different seeds without listing them.

For each instance-seed pair, first execute the recovered corrected MDM until all
70iterations finish. Preserve native final cost and its measured total process
wall duration T, including startup and native input loading. Parent start is
recorded before launch setup, as in the current timing adapter. Independent
post-exit auditing is separately timed; it must not enlarge the search budget.
Use policy original, elite10, support2absolute, stability5%, alpha0.6, no learned
ordering, and the existing declared decoder/elite corrections. Budget0retains
the iteration stopping criterion. A safety timeout43200s means censoring, not
completion. Do not replace a censored reference by its last snapshot.

After a completed reference, execute both unchanged hybrid policies with a total
incumbent horizon T: MDM initialization of min(10,T) or min(60,T), then the
existing uniform removal6/greedy repair/LKH0.2s remainder, acceptance1%,reset20.
Common fallback and all candidate setup are charged. No offline reference route
or solution value seeds the candidate. Retain only candidate events with native
timestamp<=candidate_parent_start+T. Cleanup may overrun and is reported.

This is an explicitly baseline-defined, per-pair budget: it is conditioned on
the measured duration of that reference run. It is not an independently fixed
common-horizon trial. The reference necessarily precedes candidates within a
pair; disclose this chronological asymmetry and possible machine-state drift.
Randomize instance-seed pair order before collection and randomize the two
candidate orders deterministically. Preserve the separate300-second randomized
arm-order experiment as complementary evidence. Do not conceal this limitation
or call the paired design fully randomized. No concurrent optimizers.

270scheduled executions if all90references complete; no selective seed expansion
or early stopping after favorable outcomes. If a reference is censored/invalid,
its two candidate slots have unavailable T and remain explicitly unscheduled
because of the reference failure; the intended270 denominator and missing pairs
stay in the report. Do not make a complete cost-only primary claim by omitting
them. Infrastructure interruption and genuine solver/algorithm failure are kept
distinct, with original artifacts retained and reasons for any technical rerun.

## Outcomes and analysis

For each candidate, compute100*(mean_seed_MDM70-mean_seed_candidate)/
mean_seed_MDM70 per geometry over all ten seeds, then average equally over nine
geometries. Report both candidates, best and mean cost, dispersion, all seed
pairs, failures, actual search horizons, cleanup durations and native mining
activation. Do not compare raw seconds to the old paper's processor as a speedup.

TRAIN investment criterion per candidate: primary mean>=0.2%, >=6positive
geometry means, <=1negative geometry mean, all scheduled reference/candidate
pairs valid and no feasibility regression. This is an exploratory screen, not
confirmatory significance or automatic sufficiency for a top-journal paper.
Report the time allocation contrast descriptively. Seeds measure conditional
variability; they do not multiply the number of independent geometries. No
selecting a best seed, favorable variant or alternate endpoint for the primary
comparison. A passing candidate needs a separately frozen broader evaluation.

## Time-to-target limitation

The existing untimed70-iteration binary stores only selected iteration
checkpoints; it does not supply a complete timestamped improvement stream.
Do not compare a candidate's precise target times to coarse reference checkpoints
as if both were exact. Primary endpoints here are final cost and measured time.
A faithful time-to-target study requires prospective logging/targets/seed counts,
separate verification and explicit censoring; no synthetic interpolation.

The supplied MDM manuscript describes a100-run TTT experiment on kroA150 with
seeds1-100, but the inspected text gives target12918 while its caption gives
12962. Resolve this against the final published version before reusing that
target or reproducing its plot. Do not silently choose either value.
