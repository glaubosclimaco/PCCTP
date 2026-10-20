# Completed validation of the frozen randomized removal/repair hybrid

The 70 scheduled executions completed: seven validation geometries, five seeds
45–49, two arms, 60 seconds total per run. Independent scalar and absolute-event
auditor PASS: 466 routes, no missing runs or feasibility failures.

| Instance | MDM mean | Hybrid mean | Reduction (%) |
|---|---:|---:|---:|
| g_d198_VT39_T119_W40_50 | 13606.4 | 13626.8 | -0.149929 |
| g_kroB100_VT20_T60_W20_25 | 18406.0 | 18406.0 | 0.000000 |
| g_kroE100_VT20_T20_W60_50 | 10586.0 | 10586.0 | 0.000000 |
| g_pr124_VT24_T25_W75_75 | 34359.0 | 34359.0 | 0.000000 |
| g_pr152_VT29_T31_W92_75 | 49110.0 | 49110.0 | 0.000000 |
| g_pr264_VT52_T159_W53_75 | 46320.2 | 45055.6 | 2.730126 |
| g_pr299_VT101_T99_W99_75 | 31459.2 | 30272.8 | 3.771234 |

Equal-geometry primary mean: 0.907347%. Two wins, four ties,
one loss. The prespecified screen required at least three wins, so the gate
FAILED despite a positive aggregate gain and improvements on pr264 and pr299.
Do not change the gate or claim the method is uniformly ineffective. The original
twelve test geometries remain unopened for this configuration.

The gain uses ratios of five-seed mean costs, then averages the seven geometry
effects equally. It is not a pooled route-cost reduction or a seed-level inference.
Validation geometries had already appeared in an earlier different LKH study;
this is limited screening, not untouched confirmation. No significance or novel
best-known solution is established. Header targets are not treated as verified BKS.

Training previously gave 1.886205% across nine hard geometries and conditional
fixed-support ablation favored the hybrid in nine of nine geometry means. Those
exploratory results remain valid in their stated scope; validation limits the
claim they support. The next repair experiment uses TRAIN states only and does
not retune the current method on these validation outcomes.

Raw records, native checkpoints, protocol/source/binary hashes and a completed
file manifest are retained in experiments/lkh/ruin_validation60. Audit and summary
are stored beside that directory. All late incumbents remain excluded; process
cleanup can exceed the nominal 60-second incumbent horizon.
