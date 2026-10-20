# Prospective conditional comparison of cached matrix serialization

Status: prospective design with phase adapter and scalar audit verified on
synthetic stubs; no search observations collected and no runner queued.
Written after the serialization microbenchmark, before any cached optimizer run.
The ongoing 81-run collection and seed80 reference calibration take priority.

## Question and scope

Does exact reuse of integer distance strings improve the final feasible cost
obtained by the existing uniform-removal / greedy-repair / LKH phase when both
implementations receive the same conditional start and 30 seconds total?
The byte-equivalence and serialization timings do not answer this question.
This is a TRAIN component study, excluding generation of the shared MDM starts,
not a fresh end-to-end comparison with MDM or a new-method novelty claim.

Use all nine approved hardTRAIN variants. For each, use the three existing
pre-repair routes in ruin_total60_absolute seed42/43/44 records. Bind all 27
source records to the completed independent result audit and retain exact
initial costs/routes. Assign fresh repair RNG seeds91/92/93 to source seeds
42/43/44 respectively. Both implementations use identical starting routes and
RNG seeds. Do not open validation or reserved test inputs for this study.

## Arms and accounting

Two arms: unchanged screen.polish and isolated CachedPolisher. Both use the
unchanged random destruction size6, greedy repair, LKH3.0.6 with per-call
allowance at most0.2s, fixed LKH seed42, acceptance1% and reset every20 attempts.
The cached arm constructs one cache lazily within the first polish call;
construction consumes that call's allowance and the outer phase allowance.
Subsequent calls may reuse it only for the same unchanged input matrix.
No timing learned from a reference search enlarges either 30-second budget.

Record parent monotonic start before phase setup, all incumbent absolute event
stamps, raw phase elapsed time, cache build events and counters. The original
common start is the eligible fallback at time zero. Reconstruct the primary
solution using only feasible improvements observed by parent_start+30s, even
if cleanup or the last call returns later. Preserve late events and durations
separately. Cache construction, input loading, serialization, validation within
search and LKH invocation all consume the phase budget. Independent post-exit
auditing does not extend search. Exclude no adverse seed, timeout or loss.

54 scheduled phases:9geometries x3source states x2arms, nominal27minutes.
Shuffle the27pair order deterministically with seed20261020 and shuffle the
two arm orders per pair using the fixed instance/source-seed identifier.
Use one optimizer at a time, ClockGuard and pmset sleep checks per pair.
Hold the same OS study lock as the prospective paired70 runner throughout,
so neither study can overlap the other. Require both completed/audited awake81
and nine completed/audited seed80 references before launch. Preserve and stop
on environmental failure or a partial unreviewed pair; do not silently resume.

## Primary screen and interpretation

For each geometry, compare mean final costs over its three paired starts;
compute100*(mean_uncached-mean_cached)/mean_uncached, then average the nine
geometry effects equally. Retain all54slots and report both per-run costs and
attempt counts, normal-exit/timeout/other statuses, memory expense and cache
initialization time. More polish calls are a mechanism description, not the
primary outcome or evidence of better routes.

Exploratory investment rule: mean additional reduction>=0.2%, at least six
positive geometry means, at most one negative mean, no missing/invalid phases.
No seed expansion, best-seed choice, optional stopping or alternate gate from
favorable microbenchmark timing. Ties and nonpositive results are retained.
If the cost screen fails, do not label the cache a solution-quality improvement
merely because string assembly is faster.

A positive conditional screen would still require a separately frozen complete
search comparison with all cache costs charged. It does not change or replace
either planned paired70 candidate. Any later protocol amendment must precede
its main observations and disclose this additional TRAIN selection. The twelve
reserved test geometries remain untouched; reused training starts cannot be
called independent confirmation.

## Implementation prerequisites

Before launch: implement an isolated phase runner and independent event/cost/
source/environment auditor; verify deadline selection, cache charging, call
restoration/isolation, exact source files and interrupted-pair handling using
synthetic stubs. Reuse current cached adapter certificates only if their source
hashes still match. Do not edit active screen.py, ruin_repair_absolute.py,
awake300_train.py, or calibration code to implement this study.

Implementation status on20October,12:26–12:36 Fortaleza: cached30_phase.py
implements the isolated phase and cached30_phase_audit.py independently
reconstructs eligibility and scalar feasibility/costs. verify_cached30_phase.py
has21 synthetic checks, including a raw improvement timely for the inner clock
but late for the parent deadline, an event exactly on that deadline, no-call
fallback, accounting/source/route corruptions, and restoration after failures.
An exclusive attempt marker persists after interruption or a post-search
assertion; such a phase cannot be silently retried. These checks launch zero
optimizers and do not establish a search benefit. The shared-lock batch owner,
environmental/resume audit and batch fixtures are still required before any
real execution. cached30_starts.py/JSON now binds all27starts to the completed
60-second audit and approved TRAIN input hashes, with independent scalar route
checks. cached30_pairs.py provides a pair loop without a launch entrypoint;
eight orchestration stub checks pass, including interruption/sleep preservation
and resume without duplicate calls. Its future owner must acquire the shared
lock and enforce prerequisites; that owner and an independent batch auditor
are not implemented yet. Nothing has been added to the active
collection or the planned paired70 candidates.
