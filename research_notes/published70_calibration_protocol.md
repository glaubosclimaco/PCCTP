# Reference calibration: full70-iteration MDM before a paper-scale comparison

This is a reference calibration, not the final ten-seed comparison and not a
reproduction of all107published evaluation instances. It follows the user's
request to align stopping criteria with the papers. The current300-second study
continues unchanged as an anytime supplement and must finish before this queue.

Use the same nine hardTRAIN instances, fresh seed80, MDM original policy,70full
iterations, recovered decoder-initialized/occupied-elite-bounded binary. Keep
the legacy elite capacity10, absolute support2 and stability parameter5; disclose
these recovered implementation settings and fixes. The article reports support
20%; absolute2 corresponds to20% of full capacity10, not necessarily20% of a
partially occupied bank. Do not silently claim all implementation details are
identical to the published executable. Preserve exact hashes and argument list.

Disable the experimental positive-time budget (PCCTP_BUDGET=0). The unchanged
native code then executes numIt=70. Safety process timeout43200seconds is a
censoring limit, never an eligible completed reference. Require completed=true,
iterations=70, final snapshot iteration70, clean exit and all original-input
route costs/feasibility verified. Preserve checkpoints, logs and partial native
exports. A timeout, crash or malformed export is a technical/censored result,
not a70-iteration cost and not a reason to drop the instance from the report.

Prioritize the previously censored gil262 and rd400 references, then the other
seven in lexical order. No candidate search runs during this calibration. Use
private seed-specific native work directories and one optimizer at a time.
Record native final wall time, wrapper wall time, child CPU usage, mining calls,
iteration checkpoints and input/source/binary hashes. Untimed native anytime
fields being empty are expected and must not be interpreted as no feasible route.

Do not compare seconds from the old paper's processor directly to this machine.
These full references estimate local expense and verify activation/termination.
They cannot alone establish significance or superiority. A completed calibration
is a prerequisite for a separately frozen main comparison with ten declared seeds,
all methods on the same machine, and a candidate total-time rule explicitly tied
to or calibrated from the full reference. Ten candidate iterations do not equal
ten MDM iterations; no extra candidate polishing is free. The final larger suite
and reserved data remain subject to prospective design.

If the reference cannot finish within this safety limit, report it and assess
computational feasibility rather than replacing it with the300-second baseline.
Do not tune the algorithm using reserved data or reinterpret short-budget wins
as proof of superiority over the completed published stopping criterion.
