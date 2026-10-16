"""Projective-plane Tanner codes a la Braverman-He, any GF(p^e), any twist.
Points p, lines l in T (normalized), edge (p,l) iff p.l=0.
Point local code: homogeneous deg kP-1 poly in coords of l w.r.t. basis of p-perp.
Line local code: Q(lam^s, mu^s) deg kL-1, s = Frobenius^j, coords of p w.r.t. basis of l-perp.
"""
import sys, itertools, numpy as np
from gf import GF

def build(p,e,kP,kL,j):
    F=GF(p,e); q=F.q
    T=[(0,0,1)]+[(0,1,a) for a in range(q)]+[(1,a,b) for a in range(q) for b in range(q)]
    idx={t:i for i,t in enumerate(T)}
    M,A=F.mul,F.add
    def dot(x,y):
        r=0
        for a,b in zip(x,y): r=int(A[r,M[a,b]])
        return r
    N=len(T)
    nb=[[l for l in range(N) if dot(T[pp],T[l])==0] for pp in range(N)]  # lines through point pp
    edges=[(pp,l) for pp in range(N) for l in nb[pp]]
    eid={ed:i for i,ed in enumerate(edges)}
    # lines' point lists (same incidence since symmetric dot)
    def coords(z,y):
        perp=[t for t in range(N) if dot(T[z],T[t])==0]
        u,v=T[perp[0]],T[perp[1]]
        for a,b in [(0,1),(0,2),(1,2)]:
            det=F.sub[M[u[a],v[b]],M[u[b],v[a]]]
            if det: break
        di=F.inv[det]
        # solve lam*u+mu*v=y on coords a,b
        lam=M[di, F.sub[M[y[a],v[b]],M[y[b],v[a]]]]
        mu=M[di, F.sub[M[u[a],y[b]],M[u[b],y[a]]]]
        lam,mu=int(lam),int(mu)
        # check
        for c in range(3): assert int(A[M[lam,u[c]],M[mu,v[c]]])==y[c]
        return lam,mu
    def mono(lam,mu,k):
        return [int(M[F.pw(lam,k-1-r),F.pw(mu,r)]) for r in range(k)]
    # point local generator coefficient map: value on edge (pp,l) = sum_r a[pp,r]*mono_r(coords of l wrt pp)
    ptc={}
    for pp in range(N):
        for l in nb[pp]:
            lam,mu=coords(pp,T[l]); ptc[(pp,l)]=mono(lam,mu,kP)
    lnc={}
    for l in range(N):
        for pp in nb[l]:
            lam,mu=coords(l,T[pp]); lnc[(pp,l)]=mono(F.frob(lam,j),F.frob(mu,j),kL)
    return dict(F=F,q=q,T=T,N=N,nb=nb,edges=edges,eid=eid,ptc=ptc,lnc=lnc,kP=kP,kL=kL)

def tanner_W(G):
    """basis of W (rows, over edges) from point-coefficient parametrization."""
    F,N,nb,kP,kL=G['F'],G['N'],G['nb'],G['kP'],G['kL']
    q=G['q']; rows=[]
    nunk=N*kP
    for l in range(N):
        pts=nb[l]
        Gl=np.array([[G['lnc'][(pp,l)][r] for pp in pts] for r in range(kL)],dtype=np.int16)
        H=F.nullspace(Gl)  # rows h with Gl h = 0 -> h is parity check? need h with h.c=0 for c in rowspace(Gl): Gl h^T=0
        # value vector w_pts = Cmat @ a ; constraint H w = 0
        for h in H:
            row=np.zeros(nunk,dtype=np.int16)
            for t,pp in enumerate(pts):
                if h[t]==0: continue
                coef=G['ptc'][(pp,l)]
                for r in range(kP):
                    row[pp*kP+r]=F.add[row[pp*kP+r],F.mul[h[t],coef[r]]]
            rows.append(row)
    Mx=np.array(rows,dtype=np.int16)
    K=F.nullspace(Mx)   # point-coefficient vectors
    # map to edge values
    E=len(G['edges'])
    Cm=np.zeros((nunk,E),dtype=np.int16)
    for i,(pp,l) in enumerate(G['edges']):
        for r in range(kP): Cm[pp*kP+r,i]=G['ptc'][(pp,l)][r]
    # W = K @ Cm over GF
    Wb=np.zeros((K.shape[0],E),dtype=np.int16)
    for r in range(nunk):
        c=K[:,r]
        nz=np.nonzero(c)[0]
        if len(nz)==0: continue
        Wb[nz]=F.add[Wb[nz], F.mul[c[nz][:,None],Cm[r][None,:]]]
    return Wb

if __name__=="__main__":
    p,e,kP,kL,j=map(int,sys.argv[1:6])
    G=build(p,e,kP,kL,j)
    Wb=tanner_W(G)
    E=len(G['edges']); S=Wb.shape[0]
    rk=G['F'].rank(Wb) if S else 0
    gap= 5*S/(E-S) if S<E else float('inf')
    print(f"q={G['q']} kP={kP} kL={kL} twist=Frob^{j}: N={G['N']} E={E} dimW={S} (edge-rank {rk}) S/E={S/E:.4f} nondet ratio 5S/(E-S)={gap:.4f} delta_nd={1-1/gap if gap>0 else 0:.4f}")
