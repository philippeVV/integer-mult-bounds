"""Circulant model of the Singer Tanner code. edge (i,d): point i, line i+d, d in D'.
C_P, C_L: 8x16 generator matrices over F4096 (columns indexed by D' order).
dim W = sum_{omega^273=gamma} dim(C_P ∩ Lambda_omega C_L)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir
import numpy as np
from f4k import *
from singer import singer
n,D,PW=singer()
def U_coords():
    """alpha^d (d in D) in F16-basis (u1,u2) of U=ker Tr: returns lam,mu as F4096 elements of subfield."""
    us=[PW[d] for d in D]
    u1,u2=us[0],us[1]
    sub=[0]+[int(EXP[273*k]) for k in range(15)]
    co={}
    for a in sub:
        for b in sub:
            co[int(mul(a,u1))^int(mul(b,u2))]=(a,b)
    return {d:co[PW[d]] for d in D}
def gen_twisted(Dp,frob_line=2,deg=8):
    co=U_coords()
    def pw_(x,k): return 1 if k==0 else (0 if x==0 else apow(LOG[x]*k))
    def fr(x,j):
        for _ in range(j): x=int(mul(x,x))
        return x
    CP=np.array([[int(mul(pw_(co[d][0],deg-1-r),pw_(co[d][1],r))) for d in Dp] for r in range(deg)])
    CL=np.array([[int(mul(pw_(fr(co[d][0],frob_line),deg-1-r),pw_(fr(co[d][1],frob_line),r))) for d in Dp] for r in range(deg)])
    return CP,CL
def dimW(CP,CL,Dp,gamma_log=0,per=False):
    # omega = alpha^{(gamma_log + 4095 t)/273}... omegas with omega^273 = alpha^gamma_log: omega = alpha^{e}, 273 e = gamma_log mod 4095
    es=[e for e in range(4095) if (273*e-gamma_log)%4095==0]
    assert len(es)==273
    tot=0; spec=[]
    for e in es:
        lam=np.array([apow(e*d) for d in Dp])
        M=np.vstack([CP,mul(CL,lam[None,:])])
        k=16-rank(M); tot+=k; spec.append(k)
    return (tot,spec,es) if per else tot
if __name__=="__main__":
    import sys
    d0=int(sys.argv[1]); Dp=[d for d in D if d!=d0]
    for j in [2,1,0]:
        CP,CL=gen_twisted(Dp,j)
        res=[(g,dimW(CP,CL,Dp,273*g)) for g in range(15)]
        print("frob",j,res,flush=True)
