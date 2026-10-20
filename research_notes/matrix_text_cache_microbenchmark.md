# Routing-matrix text cache: isolated prototype measurement

No optimization experiment was executed. The current 300-second collector and the prospective paired70 candidates do not import this prototype.

Nine previously audited seed42 TRAIN starting supports; five interleaved repetitions per serializer after checking exact output equality against the actual current polish writer. LKH subprocess invocation was replaced by an in-process capture stub.

| Geometry | Support | Cache build (ms) | Original serialization (ms) | Cached serialization (ms) | Ratio | Estimated reuse break-even |
|---|---:|---:|---:|---:|---:|---:|
| a280 | 197 | 49.68 | 34.27 | 7.73 | 4.43 | 2 |
| gil262 | 195 | 44.92 | 34.44 | 8.44 | 4.08 | 2 |
| gr229 | 163 | 34.65 | 25.39 | 6.53 | 3.89 | 2 |
| kroA200 | 127 | 26.38 | 14.70 | 3.50 | 4.20 | 3 |
| pr226 | 157 | 33.21 | 25.14 | 5.75 | 4.37 | 2 |
| rd400 | 256 | 141.03 | 61.77 | 12.82 | 4.82 | 3 |
| si175 | 123 | 20.18 | 14.22 | 3.14 | 4.53 | 2 |
| ts225 | 158 | 37.01 | 22.73 | 5.10 | 4.46 | 3 |
| tsp225 | 77 | 8.97 | 5.71 | 1.19 | 4.80 | 2 |

Times summarize only string assembly, excluding file writing, route validation, process startup, and LKH optimization. Cache build is a one-time measured cost, not free preprocessing. The break-even count is a calculation from these medians, not an end-to-end observation. Full-matrix cache requires finite integral costs and O(n²) string storage. Measurements occurred while the ongoing experiment was running; they are engineering diagnostics, not controlled performance claims.

The exact output match prevents changing the TSP costs or numbering in the examined cases. It does not certify a future integrated search, its deadline accounting, memory expense or solution-quality benefit. Any integrated variant needs its own source freeze, startup accounting, independent audits and equal-total-time comparison. No current method, seed, stopping rule or gate is changed. Reserved test inputs remain unopened.

Source records and all individual timings: experiments/lkh/matrix_text_cache_benchmark.json. Independent start-hash/arithmetic reconciliation: matrix_text_cache_benchmark_audit.json.
