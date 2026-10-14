"""Exhaustive: for every port t and every nondegenerate plane <u_t,r>, number of odd nbrs s with u_s.r=0
(necessary for E_t perp E_s). Exact, F2."""
from csp import *
from collections import Counter
mx=Counter()
for i in range(n):
    best=max(sum(1 for j in odd[i] if d(U[j],r)==0) for r in range(1,512) if okplane(U[i],r))
    mx[best]+=1
print('per-port max free odd degree histogram',dict(mx))
