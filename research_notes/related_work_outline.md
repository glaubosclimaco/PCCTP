# Related-work outline before prose drafting

- Problem-specific progression: reactiveGRASP2021; exact/reductions/MIP-hybrid2023;
  mining-enhancedGRASP2026. No new claim for GRASP+MIP or preprocessing alone.
- Route ordering: Helsgaun2000 motivates LKH; installed3.0.6 provenance is a
  separate artifact, not established by the2000 paper.
- Destruction/repair: Ropke/Pisinger2006 adaptive LNS supplies broad context,
  but current fixed random-removal implementation has no adaptive weights.
- MIP neighborhoods: localbranching andRINS are classical controls, already
  investigated here, not new contributions. Need avoid overgeneralizing their
  train failures under current SCF formulation/budgets.
- Tour merging: Cook/Seymour2003 is prior art; edge union optimization is not
  inherently new. Exclusion of edges restricts scope of solver bounds.
- Potential paper question: whether explicit coverage/prize repair adds useful
  value beyond strong route ordering under charged initialization and equal
  budgets. Answer depends on complete audited experiments and component control.
- No absence claim from targeted searches; verify full supplied papers before
  submission, especially strongest benchmark results and precise baseline design.
