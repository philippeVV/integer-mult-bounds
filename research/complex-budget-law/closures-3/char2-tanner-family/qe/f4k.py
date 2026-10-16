"""F_4096 log/antilog arithmetic (poly x^12+x^6+x^4+x+1, alpha=x), rank of small matrices."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir
import numpy as np
from singer import gf_mul
Q=4096; ORD=4095
EXP=np.zeros(2*ORD,dtype=np.int64); LOG=np.full(Q,-1,dtype=np.int64)
v=1
for k in range(ORD): EXP[k]=v; EXP[k+ORD]=v; LOG[v]=k; v=gf_mul(v,2)
def mul(a,b):
    a=np.asarray(a); b=np.asarray(b)
    r=EXP[(LOG[a]+LOG[b])%ORD]
    return np.where((a==0)|(b==0),0,r)
def inv(a): return int(EXP[(-LOG[a])%ORD])
def apow(k): return int(EXP[k%ORD])
def rank(M):
    M=np.array(M,dtype=np.int64); r=0; nr,nc=M.shape
    for c in range(nc):
        if r>=nr: break
        nz=np.nonzero(M[r:,c])[0]
        if len(nz)==0: continue
        i=r+nz[0]
        if i!=r: M[[r,i]]=M[[i,r]]
        M[r]=mul(inv(int(M[r,c])),M[r])
        col=M[:,c].copy(); col[r]=0; rows=np.nonzero(col)[0]
        if len(rows): M[rows]^=mul(col[rows][:,None],M[r][None,:])
        r+=1
    return r
