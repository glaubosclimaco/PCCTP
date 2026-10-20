# TRAIN experiment: mining exposure and time allocation at 300 seconds

Defined before any outcome at this horizon. The 60-second validation gate failed
despite mean improvement; no change to that gate or reopening of the twelve test
geometries. Joint insertion's conditional TRAIN gate also failed. This experiment
addresses a distinct unresolved limitation: little mining activation in the
60-second MDM control, and the opportunity cost of switching away from MDM early.

Nine unchanged hardTRAIN instances; three fresh seeds70,71,72. Three arms each
receive300seconds total including charged feasible fallback and MDM initialization:
(1) corrected elite-bounded original MDM300; (2) MDM10 + unchanged randomized
removal/greedy-repair/LKH remainder; (3) MDM60 + the same search remainder.
No MIP insertion integration, parameter adaptation, offline route pool or learned
selector. Uniform removal6, fixed0.2-second LKH calls, 1% acceptance/reset20 stay
unchanged. This is an exploratory time-allocation experiment, not a new algorithm
claim. Preserve decoder and elite-bound correction disclosures.

81executions sequentially, deterministically shuffled arm order per instance/seed,
no concurrent optimizer, no optional stopping. Each raw record retains absolute
incumbent event timestamps; eligibility is at parent_start+300s. Native cooperative
cleanup may overrun, hard safety timeout330s; late improvement is ineligible. Retain
the common timely feasible fallback if native/postprocessing returns no timely
improvement. Report actual wall durations and failures; do not equate incumbent
horizon with process termination time. No timed run uses an offline solution.

For EACH hybrid versus MDM300, primary effect is100*(mean_seed_MDM-mean_seed_H)/
mean_seed_MDM per geometry, averaged equally over9geometries. Report both complete
comparisons. A candidate's TRAIN gate requires mean>=0.2%, >=6positive geometries,
<=1negative geometry, no feasibility regression/missing/unaudited runs. Report
hybrid60 versus hybrid10 descriptively as an allocation contrast, not another
opportunity to replace the primary gate. No inferential p-value from these many
exploratory attempts. Every scheduled seed and loss remains in the report.

Record mining calls and iteration exposure from native raw/checkpoint/log output
for each arm. A broad comparison against mining-enabled MDM requires evidence of
actual activation; a longer nominal horizon alone does not establish it. If few
runs mine, disclose it and investigate the published stopping rule separately.
Passing TRAIN only justifies a new prospective evaluation design; it does not
authorize treating previously used validation geometries as untouched confirmation.
Failing preserves a negative result and does not trigger a test run.

Before launch preserve all source/dependency/binary/protocol hashes and source split
roles, and finish earlier queues. Native working directories are seed-specific;
there is only one sequential optimizer queue. Test instances remain unread.
