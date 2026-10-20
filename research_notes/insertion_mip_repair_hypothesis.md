# Possible next TRAIN-only hypothesis: joint insertion MIP repair

Prepared while the frozen validation runs; a tiny-instance prototype is now
implemented, but is NOT integrated into the method currently being evaluated.
There is no benchmark performance result for this alternative.

Question: can a small joint insertion model choose covering/prize-restoring
visits more effectively than the current randomized ratio repair? This is
classical MIP-assisted repair in spirit; no generic novelty is asserted.

Let Q be the retained simple cycle after destroying optional visits, C the
absent optional candidates, E(Q) its consecutive gaps. Binary z[v,e] inserts
v in gap e=(a,b), costing delta[v,e]=c[a,v]+c[v,b]-c[a,b]. Require each v
in at most one gap and each gap to receive at most one vertex. For each
uncovered demand w impose sum(z[v,e]:v covers w)>=1; require collected
inserted prize >=remaining prize deficit. Minimize sum(delta*z).

Validity: replace each selected original gap by a two-edge path through its
assigned candidate. Distinct candidates and distinct gaps produce one simple
cycle, even for adjacent gaps. Costs add EXACTLY, without triangle inequality.
Retained mandatory visits survive; coverage/quota rows imply feasibility.
This model is exact only for this restricted insertion neighborhood; it is
not a global PCCTP certificate or an unrestricted repair optimum.

Restriction: at most one insertion per gap can exclude feasible full repairs.
Never infer original infeasibility when this submodel is infeasible. Keep
original best incumbent and greedy fallback. A possible candidate reduction
is3cheapest gaps per absent vertex plus a matching that reinserts removed
vertices into distinct gaps, when such a matching exists. Any matching warm
start must be validated; original removal restoration is NOT automatically
representable if several removed visits shared a gap.

Before any study: enumerate all assignments on tiny retained cycles to check
objective, coverage, prize and expansion; include negative insertion costs
and adjacent-gap assignments. Compare same destruction samples, total repair
and LKH time, greedy versus joint repair and an equal-size random-slot
control. Include model construction and solver time. Freeze limits before
collecting hardTRAIN outcomes; do not tune on current validation.

Decision to implement is deferred until validation and the existing controls
are complete. This note supplies a concrete alternative if current generalization
fails, rather than silently increasing the same search budget.

Prototype addendum: the model was implemented only for independent tiny
verification while validation continued unchanged.21models matched exact
assignment checks, including negative insertion coefficients and a feasible
original cycle whose restricted insertion model is infeasible. No hardTRAIN,
validation or test experiment of this prototype has been run. Integration
decision remains deferred until the frozen validation completes.

Boundary verification addendum: verify_insertion_mip_boundaries.py independently
enumerates assignments without using the shared route validator or row builder.
All60additional small models matched, covering 2500 assignments, varying retained
cycle order/size, zero prizes/costs, empty demand neighborhoods, unavailable slots,
missing candidate keys and no absent vertices. An expired feasible identity repair
was correctly ineligible. This is finite implementation evidence only; it supplies
neither a global exact solver nor evidence of improvement on PCCTP benchmarks.
