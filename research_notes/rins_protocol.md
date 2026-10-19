# RINS after GRASP/MDM: frozen exploratory protocol

Nine hardTRAIN, exact same starts as edge_pool/full60. Total60s includes model
setup, LP and subMIP. LP cap10s; if not optimal no agreement fixing (explicit
full-model fallback). Fix only x/y binaries agreeing with incumbent to1e-6,
never continuous flows. Standard HiGHS presolve,threads1,seed42. LP fixed
variables/values recorded, incumbent must remain feasible. External native
clock late-discard and original raw route auditing remain mandatory.

Gate: >=3 improving geometries,mean>=0.2%, plus>=3 better costs than full60.
Offline input generation excluded; no end-to-end claim. Bounds/status local
to fixed neighborhood, not global. No validation/test access. Do not loosen
tolerance or increase LP budget after observing hardTRAIN outcomes.

Classical RINS: Danna,Rothberg,LePape,MathProg102:71-90(2005),
DOI10.1007/s10107-004-0518-7. Primary publisher abstract checked. No novelty
claim for RINS itself. One6vertex model independently enumerated after
agreement fixing; incumbent preserved and restricted optimum matched.
