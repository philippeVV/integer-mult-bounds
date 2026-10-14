"""Frame skeleton (r in W: every frame term cost 8, total 64) + corrections.
Correction X: X_aa=0, X_ab=0 on even pairs, X_ab=-G_ab/2 on lost odd pairs, free on plane-orth pairs.
Float ALS masked minrank of X (and of full masked identity) for r in W maximising free edges."""
import sys, json, random, math
import numpy as np
from sa_class import *
FL=sorted(set(emap)-set(P)); wl=[x for x in FL]  # frame labels in same order as frame
W=[[r for r in RC[a] if all(dot(r,wl[k])==0 for k in Ks[a])] for a in range(n)]
print('W sizes',set(len(w) for w in W))
def free_count(r,a): return sum(1 for b in ODD[a] if not lost(r,a,b))
def sa(seed,iters=60000):
    rng=random.Random(seed); r=[rng.choice(W[a]) for a in range(n)]
    tot=sum(free_count(r,a) for a in range(n))//2
    for it in range(iters):
        T=2*(0.05/2)**(it/iters); a=rng.randrange(n); old=r[a]; f0=free_count(r,a)
        r[a]=rng.choice(W[a]); d=free_count(r,a)-f0
        if d>=0 or rng.random()<math.exp(d/T): tot+=d
        else: r[a]=old
    return r,tot
def als(req, tgt, d, iters=300, seed=0, lam=1e-9):
    rng=np.random.default_rng(seed); U=rng.normal(size=(n,d)); V=rng.normal(size=(n,d))
    for it in range(iters):
        for a in range(n):
            m=req[a]; A=V[m]; U[a]=np.linalg.solve(A.T@A+lam*np.eye(d), A.T@tgt[a,m])
        for b in range(n):
            m=req[:,b]; A=U[m]; V[b]=np.linalg.solve(A.T@A+lam*np.eye(d), A.T@tgt[m,b])
    R=(U@V.T-tgt)*req
    return float(np.sqrt((R**2).sum()))
if __name__=="__main__":
    seed=int(sys.argv[1]); r,tot=sa(seed)
    print('free odd edges',tot, 'lost', 2912//1-tot if False else sum(len(ODD[a]) for a in range(n))//2-tot, flush=True)
    # frame cost
    fc=[term_cost(r,D[k],D[k]) for k in range(8)]; print('frame term costs',fc)
    G=np.array([[ip(Rt[a],Rt[b]) for b in range(n)] for a in range(n)],dtype=float)
    F=np.array(free_matrix(r)); req=~F
    tgtX=np.zeros((n,n)); 
    for a in range(n):
        for b in ODD[a]:
            if lost(r,a,b): tgtX[a,b]=-G[a,b]/2
    tgtI=np.eye(n)
    json.dump({'r':r},open('rW_%d.json'%seed,'w'))
    for d in [int(x) for x in sys.argv[2].split(',')]:
        print('d',d,'resid corrections',round(als(req,tgtX,d),5),'resid masked identity',round(als(req,tgtI,d),5),flush=True)
    