# Passive call diagnostics of completed 60-second studies

Completed 60-second studies only; passive cumulative counters, not per-call timing, causal attribution, or successful-tour counts.

Counters classify polish invocations; normal exit is not a successful or improving tour.

| Study / geometry | Invocations | Mean per run | Normal exits | Timeouts | Other/absent exit |
|---|---:|---:|---:|---:|---:|
| TRAIN: a280 | 1320 | 440.0 | 1303 | 15 | 2 |
| TRAIN: gil262 | 1540 | 513.3 | 1532 | 7 | 1 |
| TRAIN: gr229 | 1529 | 509.7 | 1516 | 11 | 2 |
| TRAIN: kroA200 | 2452 | 817.3 | 2448 | 3 | 1 |
| TRAIN: pr226 | 1978 | 659.3 | 1975 | 2 | 1 |
| TRAIN: rd400 | 735 | 245.0 | 457 | 277 | 1 |
| TRAIN: si175 | 2948 | 982.7 | 2946 | 2 | 0 |
| TRAIN: ts225 | 1763 | 587.7 | 1756 | 7 | 0 |
| TRAIN: tsp225 | 5513 | 1837.7 | 5510 | 3 | 0 |
| Reused validation: d198 | 1927 | 385.4 | 1818 | 107 | 2 |
| Reused validation: kroB100 | 11363 | 2272.6 | 11359 | 4 | 0 |
| Reused validation: kroE100 | 23079 | 4615.8 | 23077 | 2 | 0 |
| Reused validation: pr124 | 20112 | 4022.4 | 20111 | 0 | 1 |
| Reused validation: pr152 | 10453 | 2090.6 | 10449 | 4 | 0 |
| Reused validation: pr264 | 2059 | 411.8 | 2030 | 28 | 1 |
| Reused validation: pr299 | 3395 | 679.0 | 3390 | 5 | 0 |

No optimizer, active 300-second outcome or reserved test input was read.

The other/absent-exit bucket cannot distinguish crashes from setup-only calls. Per-call output validity, improvements and CPU totals cannot be reconstructed from these counters. No operator or budget is changed in response to this report.


## Interpretation restricted to the observed records

In TRAIN rd400, 277/735 invocations report a subprocess timeout (37.69%),
accounting for 277/327 TRAIN timeouts (84.71%). Its end-to-end mean cost
reduction was nevertheless positive in the completed 60-second study;
these counters do not quantify how many gains a longer LKH call would add.
The validation timeout proportion is 150/72388 (0.2072%). In particular,
the four tied validation geometries have at most four timeouts each, while
recording many normal exits. Widespread subprocess timeouts therefore do
not describe the observed execution status in those tied cases. This is
not a causal explanation of their ties or evidence of support optimality.

A future allocation comparison could vary the per-call ordering allowance
against the number of destruction/repair attempts, with all setup charged.
It would need a separate frozen design, controls and seeds. This diagnostic
does not schedule that experiment, change the active 300-second study or
replace the pending completed-70-iteration reference comparison. Reused
validation outcomes cannot support a new untouched-confirmation claim.
