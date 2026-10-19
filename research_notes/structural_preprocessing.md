# Structural preprocessing investigation

Implemented independently: indispensable prize vertices, singleton coverage,
coverage already satisfied by mandatory vertices, duplicate coverage and
inclusion dominance of coverage rows. No metric assumption, no vertex deletion.
If N1 is contained in N2, requiring a visit in N1 implies coverage of N2.
Indispensable prize uses only visitable vertices, never demand prizes.

Finite exhaustive subset check:20320 subsets in140 random systems, seed42.
All nine hardTRAIN cases lose coverage rows; five gain mandatory vertices.
Four have all demands already covered by mandatory vertices. These counts
do not establish runtime improvement: standard solver presolve may reproduce
these reductions. Keep presolve enabled in both arms of any comparison.

Legacy regras_reducao.cpp contains rules01 and02 for mandatory coverage and
indispensable prize. No novelty claimed. include/Regras.h declares R1-R14
with missing-implementation comment. Legacy rule02 tests type of i within
the j-loop: investigate reachability and effect before calling it a result.
Do not modify user sources or silently change the baseline.

Next experiment after full_graph60: same full model, starts and60s including
preprocessing, all nineTRAIN, default HiGHS presolve. Report setup and solver
work plus quality; count reductions alone are insufficient. Test stays closed.
