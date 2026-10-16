"""expQE core: delay-network view of the Singer Tanner code on PG(2,16)-M(d0=59).
dim W = sum_{eta in mu_273} nullity( H_P diag(eta^{y_d}) G_L^T ),  y_d = 118 d - 3 (centred)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir
import numpy as np
from f4k import EXP, LOG, mul, inv, apow, rank, ORD
from singer import singer
n,D,PW=singer()
d0=59; Dp=[d for d in D if d!=d0]
Y=np.array([((118*d-3+136)%273)-136 for d in Dp])
MU=[apow(15*t) for t in range(273)]          # mu_273
def det(M):
    M=np.array(M,dtype=np.int64); k=M.shape[0]; dt=1
    for c in range(k):
        nz=np.nonzero(M[c:,c])[0]
        if len(nz)==0: return 0
        i=c+nz[0]
        if i!=c: M[[c,i]]=M[[i,c]]
        p=int(M[c,c]); dt=int(mul(dt,p)); M[c]=mul(inv(p),M[c])
        col=M[:,c].copy(); col[c]=0; rows=np.nonzero(col)[0]
        if len(rows): M[rows]^=mul(col[rows][:,None],M[c][None,:])
    return dt
def nullspace(M):
    """rows spanning {x: M x^T = 0}"""
    M=np.array(M,dtype=np.int64); r=0; nr,nc=M.shape; piv=[]
    for c in range(nc):
        if r>=nr: break
        nz=np.nonzero(M[r:,c])[0]
        if len(nz)==0: continue
        i=r+nz[0]
        if i!=r: M[[r,i]]=M[[i,r]]
        M[r]=mul(inv(int(M[r,c])),M[r])
        col=M[:,c].copy(); col[r]=0; rows=np.nonzero(col)[0]
        if len(rows): M[rows]^=mul(col[rows][:,None],M[r][None,:])
        piv.append(c); r+=1
    free=[c for c in range(nc) if c not in piv]
    B=np.zeros((len(free),nc),dtype=np.int64)
    for t,f in enumerate(free):
        B[t,f]=1
        for i,pc in enumerate(piv): B[t,pc]=M[i,f]
    return B
def lampow(lam,ys):
    if lam==0: raise ValueError
    return np.array([apow(int(LOG[lam])*int(y)) for y in ys])
def Mmat(HP,GL,lam):
    return matmul(mul(HP,lampow(lam,Y)[None,:]),GL.T)
def matmul(A,B):
    A=np.asarray(A); B=np.asarray(B)
    out=np.zeros((A.shape[0],B.shape[1]),dtype=np.int64)
    for k in range(A.shape[1]): out^=mul(A[:,k][:,None],B[k][None,:])
    return out
def spectrum(GP,GL,pts=MU):
    HP=nullspace(GP)
    return [8-rank(Mmat(HP,GL,l)) for l in pts]
def dimW(GP,GL): return sum(spectrum(GP,GL))
def detpoly(GP,GL):
    """p(x)=x^{348} det(H_P X G_L^T), coefficients c[0..784] (exact via DFT over F4096^*)"""
    HP=nullspace(GP); vals={}
    lams=[apow(t) for t in range(ORD)]
    v=np.array([det(Mmat(HP,GL,l)) for l in lams]); v=mul(v,np.array([apow(348*t) for t in range(ORD)]))
    c=[]
    for k in range(ORD):
        s=0
        w=mul(v,np.array([apow(-k*t) for t in range(ORD)]))
        s=int(np.bitwise_xor.reduce(w)); c.append(s)
    return c
