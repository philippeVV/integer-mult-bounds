"""Exact single-column replacement neighbourhood.
Replace column d of G_L by an arbitrary v in F4096^8:  dim W(v) = base_d + #{eta : v in V_eta}, V_eta affine of codim s_eta.
(column of H_P handled by transposing roles: swap (H_P,G_L) <-> (G_L,H_P) since M^T = G_L Z H_P^T.)"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir
import numpy as np, random
from f4k import mul, inv, apow, LOG, EXP
from lib import Y, nullspace, rank
import fast
T=np.arange(273)
def Mstack(H,G,skip=None):
    A=fast.bmul(H[None,:,:],fast.POW[:,None,:])
    M=np.zeros((273,H.shape[0],G.shape[0]),dtype=np.int64)
    for d in range(G.shape[1]):
        if d==skip: continue
        M^=fast.bmul(A[:,:,d][:,:,None],G[:,d][None,None,:])
    return M
def systems(H,G,d):
    A=Mstack(H,G,skip=d); hd=H[:,d]; base=0; eqs=[]
    for t in range(273):
        K=nullspace(np.hstack([A[t],hd[:,None]]))   # rows (m,t)
        s=K.shape[0]; base+=s-1
        sc=apow(-15*t*int(Y[d]))
        eqs.append((K[:,:8].copy(), mul(K[:,8],sc)))
    return base,eqs
def count(eqs,v):
    c=0
    for Mm,r in eqs:
        if Mm.shape[0]==0: c+=1; continue
        lhs=np.bitwise_xor.reduce(mul(Mm,v[None,:]),axis=1)
        if np.array_equal(lhs,r): c+=1
    return c
def solve(rows,rhs):
    """solve rows v = rhs; returns (v, rank) with free vars random, or None if inconsistent"""
    A=np.hstack([np.array(rows),np.array(rhs)[:,None]]).astype(np.int64); nr=A.shape[0]; r=0; piv=[]
    for c in range(8):
        if r>=nr: break
        nz=np.nonzero(A[r:,c])[0]
        if len(nz)==0: continue
        i=r+nz[0]; A[[r,i]]=A[[i,r]]; A[r]=mul(inv(int(A[r,c])),A[r])
        col=A[:,c].copy(); col[r]=0; rr=np.nonzero(col)[0]
        if len(rr): A[rr]^=mul(col[rr][:,None],A[r][None,:])
        piv.append(c); r+=1
    if np.any(A[r:,8]!=0): return None,r
    v=np.zeros(8,dtype=np.int64)
    free=[c for c in range(8) if c not in piv]
    for c in free: v[c]=random.randrange(4096)
    for i,c in enumerate(piv):
        v[c]=A[i,8]^int(np.bitwise_xor.reduce(mul(A[i,free],v[free]))) if free else A[i,8]
    return v,r
def ransac(eqs,iters,rng):
    idx=[t for t in range(273) if eqs[t][0].shape[0]>0]
    best=(-1,None)
    for it in range(iters):
        rows=[];rhs=[]; rk=0; used=[]
        for t in rng.sample(idx,len(idx)):
            Mm,r=eqs[t]
            v,rk2=solve(rows+list(Mm),rhs+list(r))
            if v is None: continue
            rows+=list(Mm); rhs+=list(r); used.append(t); rk=rk2
            if rk>=8: break
        v,_=solve(rows,rhs)
        c=count(eqs,v)
        if c>best[0]: best=(c,v)
    return best
