"""Sub-root systems R of E8 (A8, D8, E7+A1, A4^2, D4^2) as subsets of P120, via simple-root bases and closure.
T_R = {(a,b) in F : a-b or a+b in R}.  Reports block structure (Q(R)-cosets), degrees, independence number lower bounds."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, itertools, sys
from e8 import P, R, G, F, n
idx={p:i for i,p in enumerate(P)}
def key(v):
    v=tuple(int(x) for x in v); m=tuple(-x for x in v)
    return idx.get(v, idx.get(m))
def closure(simple):
    """all E8 roots in the Z-span of the given simple roots (sub-root system generated: use Weyl-closure)"""
    S=[np.array(s) for s in simple]; found={key(s) for s in S}
    frontier=list(found)
    while frontier:
        new=[]
        for a in frontier:
            for s in S:
                ra=R[a]; ip=int(ra@s)//4
                v=ra-ip*s
                k=key(v)
                if k is not None and k not in found: found.add(k); new.append(k)
        frontier=new
    return sorted(found)
# extended E8 Dynkin: build E8 simple roots (doubled coords, Bourbaki)
def vec(*x): return np.array([int(2*t) for t in x])
a1=vec(.5,-.5,-.5,-.5,-.5,-.5,-.5,.5); a2=vec(1,1,0,0,0,0,0,0)
a3=vec(-1,1,0,0,0,0,0,0); a4=vec(0,-1,1,0,0,0,0,0); a5=vec(0,0,-1,1,0,0,0,0)
a6=vec(0,0,0,-1,1,0,0,0); a7=vec(0,0,0,0,-1,1,0,0); a8=vec(0,0,0,0,0,-1,1,0)
simple=[a1,a2,a3,a4,a5,a6,a7,a8]
assert all(key(s) is not None for s in simple)
# highest root and its negative a0 (extended node). Bourbaki marks: 2,3,4,6,5,4,3,2
marks=[2,3,4,6,5,4,3,2]
theta=sum(m*s for m,s in zip(marks,simple)); a0=-theta
assert key(a0) is not None
# Borel-de Siebenthal: remove node i with mark m -> maximal-rank subsystem.
# extended diagram: a0 - a8 - a7 - a6 - a5 - a4 - a3 - a1 ; a2 attached to a4.
ext={0:a0,1:a1,2:a2,3:a3,4:a4,5:a5,6:a6,7:a7,8:a8}
def sub(remove): return closure([ext[i] for i in ext if i not in remove])
SUBS={'D8':sub({1}),'A8':sub({2}),'E7+A1':sub({8}),'A4^2':sub({5}),'D4^2':None}
# D4+D4: D8's D4xD4 (remove node of mark 2 in extended D8) -- take roots of D8 = +-e_i+-e_j with i,j in first 4 or last 4 of the D8 coords
D8=SUBS['D8']
def coordsD8():
    # D8 sublattice basis-free: use direct definition in standard coords: roots (2,2,0..)-type with support in {0..3} or {4..7}
    return [i for i,p in enumerate(P) if sum(1 for x in p if x)==2 and (all(k<4 for k,x in enumerate(p) if x) or all(k>=4 for k,x in enumerate(p) if x))]
SUBS['D4^2']=coordsD8()
def Tgraph(Rset):
    inR=np.zeros(n,bool); inR[list(Rset)]=True
    T=np.zeros((n,n),bool)
    for a in range(n):
        for b in range(n):
            if F[a,b]:
                c=key(R[a]-G[a,b]*R[b])   # third root a -+ b
                T[a,b]=inR[c]
    return T
if __name__=='__main__':
    for name,Rs in SUBS.items():
        T=Tgraph(Rs); assert (T==T.T).all()
        # components
        comp=-np.ones(n,int); c=0
        for s in range(n):
            if comp[s]>=0: continue
            st=[s]; comp[s]=c
            while st:
                a=st.pop()
                for b in np.nonzero(T[a])[0]:
                    if comp[b]<0: comp[b]=c; st.append(b)
            c+=1
        sizes=sorted(np.bincount(comp).tolist(),reverse=True)
        print(name,'|R+|=%d'%len(Rs),'d mean %.2f max %d min %d'%(T.sum(1).mean(),T.sum(1).max(),T.sum(1).min()),'components',sizes[:8], '...' if len(sizes)>8 else '')
