# GRASP/MDM + local branching, frozen TRAIN diagnostic

Nine hardTRAIN, same original edge_pool starts,60s extra including setup.
Radius12 on Hamming distance of edge AND visit binaries, no continuous flow
variables counted. All graph edges permitted, no predetermined released vertices.
Warmstart and fallback preserved. Compare with completed full_graph60 and
restricted pool. Gate>=3improvements on entry and mean>=0.2%, plus at least
three lower costs than full_graph60. No validation/test access.

This is a classical local branching implementation, not a new method claim.
Fischetti/Lodi2003 DOI10.1007/s10107-003-0395-5 verified author/institution
metadata. RINS Danna/Rothberg/LePape2005 DOI10.1007/s10107-004-0518-7
verified publisher abstract: uses relaxation to construct MIP neighborhood.
RINS is a subsequent distinct candidate, not implemented or claimed here.

Local branching differs from fixed released-vertex repair: solver chooses
which binary variables change within a radius. Three6vertex neighborhoods
(radius0/4/12) matched independent exhaustive enumeration before this study.
Bounds and Optimal status are local to Hamming neighborhood. Input generation
is excluded, hence offline not full algorithm time comparison.
