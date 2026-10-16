"""Bit-exact check of schedule784.json with Singer-invariant local codes C_P (all points), C_L (all lines) over F4096.
Simulates the causal network: session values on the 784 R arcs are free symbols; each vertex applies its local bijection
(in-arcs -> out-arcs). Phi maps session symbols -> values recomputed at the R arcs' tails.
Sound & exact  <=>  Phi = I  <=>  dim W = 784.  In general dim W = 784 - rank(Phi - I)  (W ~ fixed points of the monodromy).
Also reports gate invertibility at every vertex."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir
import json, sys, numpy as np
from f4k import mul, inv, rank
from lib import nullspace, Dp, n
S=json.load(open(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'data', 'schedule784.json')))   # expNI schedule784.json
assert S['d0']==59 and len(S['R'])==784
V=2*n; Didx={d:k for k,d in enumerate(Dp)}
edges=[(i,n+(i+d)%n) for i in range(n) for d in Dp]           # ei = 16*i + class index
cls=[k for i in range(n) for k in range(16)]
arcs=[tuple(a) for a in S['arcs']]; R=S['R']; Rk={e:k for k,e in enumerate(R)}
pos={v:t for t,v in enumerate(S['order'])}
inc=[[] for _ in range(V)]
for ei,(a,b) in enumerate(edges): inc[a].append(ei); inc[b].append(ei)
def solve_inv(A):
    """inverse of square A over F4096 (None if singular)"""
    k=A.shape[0]; M=np.hstack([A,np.eye(k,dtype=np.int64)])
    for c in range(k):
        nz=np.nonzero(M[c:,c])[0]
        if len(nz)==0: return None
        i=c+nz[0]; M[[c,i]]=M[[i,c]]; M[c]=mul(inv(int(M[c,c])),M[c])
        col=M[:,c].copy(); col[c]=0; rr=np.nonzero(col)[0]
        if len(rr): M[rr]^=mul(col[rr][:,None],M[c][None,:])
    return M[:,k:]
def matmul(A,B):
    out=np.zeros((A.shape[0],B.shape[1]),dtype=np.int64)
    for k in range(A.shape[1]): out^=mul(A[:,k][:,None],B[k][None,:])
    return out
def run(GP,GL,verbose=True):
    val=np.zeros((len(edges),784),dtype=np.int64); have=[False]*len(edges)
    for e,k in Rk.items(): val[e,k]=1; have[e]=True        # head of session arc sees the source symbol
    Phi=np.zeros((784,784),dtype=np.int64); cache={}; bad=0
    for v in S['order']:
        G=GP if v<n else GL
        ins=[e for e in inc[v] if arcs[e][1]==v]; outs=[e for e in inc[v] if arcs[e][0]==v]
        assert len(ins)==8==len(outs)
        I=[cls[e] for e in ins]; O=[cls[e] for e in outs]
        key=(v<n,tuple(I),tuple(O))
        if key not in cache:
            Ai=solve_inv(G[:,I].T.copy())
            Ao=solve_inv(G[:,O].T.copy())
            cache[key]=None if Ai is None or Ao is None else matmul(G[:,O].T,Ai)
        T=cache[key]
        if T is None: bad+=1; continue
        assert all(have[e] for e in ins)
        out=matmul(T,val[ins])
        for r,e in enumerate(outs):
            if e in Rk: Phi[:,Rk[e]]=out[r]       # tail of a session arc: recomputed value (sink)
            else: val[e]=out[r]; have[e]=True
    if bad: print("non-invertible gates:",bad); return None
    D=Phi.T.copy(); D[np.arange(784),np.arange(784)]^=1
    rk=rank(D); dW=784-rk
    if verbose: print(f"gate types {len(cache)} all bijective; rank(Phi-I)={rk}  => dim W = {dW};  sound/exact: {rk==0}")
    return dW
if __name__=="__main__":
    from circ import gen_twisted
    from f4k import apow
    CP,CL=gen_twisted(Dp,2)
    GL=mul(CL,np.array([apow(9*d) for d in Dp])[None,:])
    import fast; fw=fast.dimW(CP,GL); print("Fourier dim W:",fw)
    dw=run(CP,GL)
    print('PASS_EXPQE_SCHED' if fw==772 and dw==772 else 'FAIL_EXPQE_SCHED')
