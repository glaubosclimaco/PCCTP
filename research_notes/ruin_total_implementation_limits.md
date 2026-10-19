# Implementation limits discovered in total-time diagnostic

Original outputs preserved. The conservative event sensitivity currently has
a large loss on rd400: native MDM10 exported no predeadline feasible route,
so hybrid starts from much poorer feasible fallback. Continuing MDM60 gave
better solutions. All three ruin trajectories stabilize13975, not10523header.
This refutes extrapolating strong-start offline gains to arbitrary starts.
Do not exclude rd400 or replace it by its expensive70-iteration incumbent.

A further possible bottleneck is0.2s LKH calls including dense matrix-file
construction on large supports. No timing evidence yet proves how often calls
time out. Add passive counts in prospective absolute collector: completed,
timeout and other exit failures. No change in search decisions or time budget.
Do not lengthen calls based on this conjecture during the frozen repetition.

After current and corrected batches, required component control is fixed-support
LKH with matched total-time accounting, plus native MDM. Before new trials
state their hypotheses separately; do not attribute routing-only gain to
coverage/prize repair without this ablation.
