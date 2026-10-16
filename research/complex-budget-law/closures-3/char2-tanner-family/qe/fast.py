"""batched ranks over F4096 + fast dim W for the d0=59 Singer graph."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir
import numpy as np
from f4k import EXP, LOG, mul, apow
from lib import Y, nullspace, Dp
def bmul(a,b):
    r=EXP[(LOG[a]+LOG[b])%4095]
    return np.where((a==0)|(b==0),0,r)
def brank(M):
    M=M.copy(); B,R,C=M.shape; r=np.zeros(B,dtype=np.int64); ar=np.arange(B)
    rows=np.arange(R)[None,:]
    for c in range(C):
        col=M[:,:,c]; ok=(col!=0)&(rows>=r[:,None])
        has=ok.any(1)
        if not has.any(): continue
        i=np.argmax(ok,1)
        b=ar[has]; ib=i[has]; rb=r[has]
        tmp=M[b,ib].copy(); M[b,ib]=M[b,rb]; M[b,rb]=tmp
        piv=M[b,rb,c]; invp=EXP[(-LOG[piv])%4095]
        M[b,rb]=bmul(invp[:,None],M[b,rb])
        f=M[b,:,c].copy(); f[np.arange(len(b)),rb]=0
        M[b]^=bmul(f[:,:,None],M[b,rb][:,None,:])
        r[has]+=1
    return r
T=np.arange(273)
POW=np.array([[apow(15*t*int(y)) for y in Y] for t in T])  # eta_t^{y_d}, eta_t=alpha^{15t}
def spec(GP,GL):
    HP=nullspace(GP)
    # M_t = HP diag(POW[t]) GL^T : (273,8,8)
    A=bmul(HP[None,:,:],POW[:,None,:])            # (273,8,16)
    M=np.zeros((273,HP.shape[0],GL.shape[0]),dtype=np.int64)
    for d in range(GL.shape[1]): M^=bmul(A[:,:,d][:,:,None],GL[:,d][None,None,:])
    return GL.shape[0]-brank(M)
def dimW(GP,GL): return int(spec(GP,GL).sum())
