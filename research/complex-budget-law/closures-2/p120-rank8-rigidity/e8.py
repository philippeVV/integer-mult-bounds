"""E8 positive-root (P120) model. Roots in doubled coordinates (Z^8, norm 8); G = Gram/4 (diag 2, offdiag 0,+-1)."""
import itertools, numpy as np
def roots():
    R=[]
    for i,j in itertools.combinations(range(8),2):
        for a,b in itertools.product((1,-1),repeat=2):
            v=[0]*8; v[i],v[j]=2*a,2*b; R.append(tuple(v))
    for s in itertools.product((1,-1),repeat=8):
        if s.count(-1)%2==0: R.append(s)
    P=[]
    for r in R:
        if tuple(-x for x in r) not in P: P.append(r)
    return P
P=roots(); n=len(P)
R=np.array(P,dtype=np.int64)          # doubled coords
G=(R@R.T)//4                           # Gram, entries 2,0,+-1
F=(G%2==1)                             # free graph = non-orthogonal pairs (56-regular)
ORTH=(G==0)
