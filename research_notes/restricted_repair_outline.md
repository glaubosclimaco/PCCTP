# Structure of the restricted-repair methods draft

- Specify incumbent, released optional visits U, and free edge set F jointly.
- Keep the incumbent as a mathematical fallback, independent of solver return.
- Distinguish restricted objective bounds from global PCCTP bounds.
- Describe fixed-path components, mandatory endpoints and fixed-zero edges.
- Preserve original endpoint incidences, coverage, prizes and fixed edge cost.
- Explain contraction and expansion of integer tours, including closed-cycle case.
- State that representation equivalence does not imply novelty or acceleration.

Draft remains subject to final source and formulation audit. Existing independent
contraction checks are recorded in experiments/ml_lns/verify_fragment.py; they
support finite cases and do not replace the argument for general instances.
