# Initialization and engine timing sensitivity, observed before ablation

Original and corrected repetitions use the same search parameters but fresh
native runs. They must not be described as identical search trajectories with
only different output extraction. In rd400 seeds42/43, original native10 had
no timely feasible export and hybrid started at14143 (262visits); corrected
runs exported11298/12282 before the native deadline (256/261visits).
Thus later quality differences also reflect initial-state differences.

Corrected rd400 first two hybrid traces had475LKHcalls:241successful exits,
233subprocess timeouts and one without successful exit. The earlier pure
conjecture that calls always time out is false for this corrected collection.
These two runs are a partial observation, not a universal timing measurement.

Fixed short initialization and0.2s ordering limits can expose performance
sensitivity to execution timing. Report this limitation and all fallback cases.
Do not tune away rd400, combine best repetitions, or selectively choose better
seeds. The upcoming conditional ablation uses EXACT saved starts and remaining
allowances, separating component effects from this initialization variation.
Any future adaptive initialization policy needs its own frozen protocol and
all-case evaluation; no such new policy is implemented in the current batch.
