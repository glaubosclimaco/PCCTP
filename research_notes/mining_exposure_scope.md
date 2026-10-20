# Scope of the60-second MDM comparison

Inspection of all54 completed native raw exports in the corrected collection:
MDM arm performed at least one mining call in2/27 runs (1-33iterations);
hybrid initialization performed none in27runs (1-5iterations). The60-second
comparison is therefore against the recovered MDM implementation under a short
deadline, mostly before its mining phase activates. It does NOT isolate or
establish superiority over the data-mining component.

Preserve all cases; do not select only mining-active runs. Report this exposure
in the final reproducibility/limitations analysis. The original fixed70 study
and long-budget sensitivity are relevant complements, not substitutes for the
charged-time comparison. Before broad publication claims, evaluate longer
horizons or the published stopping criterion and account for runtime/quality.
The current primary gate and component ablation remain unchanged.


## Timing refinement from passive export bounds (19 October)

The cumulative counts above describe native exports, not timestamped mining
calls. Cooperative final iterations can overrun the horizon. Do not state that
both observed activations occurred within the60-second eligible interval.
The passive profiler finds observed calls in2/27MDM runs, but proves a timely
call from a checkpoint exported before the deadline in only1/27. The other
run's timing is unresolved within this instrumentation. All27hybrid starts
have zero calls in completed exports. This refinement changes no route cost,
quality comparison or previously frozen gate.

profile_mining_exposure.py reports lower bounds from exports with seconds<=
native_budget and upper bounds from completed cumulative exports. If no final
export exists, the upper bound stays unknown. Final calls and log return counts
are reported separately; invocation does not establish successful useful mining.
Three synthetic fixtures passed: completed overrun, truncated raw with a late
checkpoint, and completed zero-call run. No active optimizer was modified.
Use the same passive analysis for the300-second experiment, keeping partial
counts distinct from the final81-run denominator.
