from itertools import combinations
import numpy as np
n=9
blocks=[{0}]+[{2*i+1,2*i+2} for i in range(4)]
g={v:k for k,B in enumerate(blocks) for v in B}
x={}
for i,j in combinations(range(n),2):
    x[i,j]=.625 if g[i]==g[j] else (.25 if 0 in (i,j) else .0625)
y=[1]+[.625]*8
demands=[blocks[i]|blocks[j] for i,j in combinations(range(1,5),2)]
cut=lambda S:sum(v for (i,j),v in x.items() if (i in S)!=(j in S))
issues=[]
for mask in range(1,2**n-1):
    S={i for i in range(n) if mask>>i&1}
    val=cut(S)
    U=S if 0 not in S else set(range(n))-S
    if val+1e-8<2*max(y[i] for i in U):issues.append(('SEC',S,val))
    if any(D<=U for D in demands) and val<2-1e-8:issues.append(('cover cut',S,val))
    crossing=[e for e in x if (e[0] in S)!=(e[1] in S)]
    for F in combinations(crossing,3):
        if len(set(sum((list(e) for e in F),[])))<6:continue
        rhs=2*sum(x[e] for e in F)-2
        if val+1e-8<rhs:issues.append(('2match',S,F,val,rhs))
print('issues',issues[:5], 'count',len(issues))
print('partition LHS',sum(v for (i,j),v in x.items() if g[i]!=g[j]),'RHS',4)
print('min demand coverage',min(sum(y[i] for i in D) for D in demands))
print('degrees',[(i,sum(v for e,v in x.items() if i in e),2*y[i]) for i in range(n)])
