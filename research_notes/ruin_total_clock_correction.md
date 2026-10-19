# Prospective event-clock collector correction

Original ruin_total60 batch preserved unchanged, including overly strict return
acceptance. Conservative reconstructed events reported separately as sensitivity.
Prepare a fresh repetition with identical cases/seeds/parameters and explicit
CLOCK_MONOTONIC absolute timestamps on every accepted post event. Wrapper keeps
any improvement timestamp<=parent_start+60, regardless of subroutine cleanup
return time. This corrects extraction, not search policy. No tuning of arms,
removal count, acceptance threshold, native budget or gate. Do not mix the two
runs or count them as independent geometries. Save parent/post starts to audit.

ruin_repair_absolute.py and ruin_total60_absolute.py are separate source files
so active original worker and its source hashes stay unchanged. Guard requires
original completion; do not launch concurrently. After fresh batch audit absolute
event times and source hashes, original raw feasibility and all costs.

Before prospective repetition, removed irrelevant offline-pool metadata reads
when initial_override is supplied. The search parameters and supplied route
are unchanged; IO overhead remains charged. New collector passed exact deadline
boundary checks. Independent auditor reconstructs fallback/native choice and
all timestamp-qualified improvements from original scalar input.
