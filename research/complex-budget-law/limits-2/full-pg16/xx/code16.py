"""Hermitian F16 Tanner code of expNC (pg.py 2 4 8 9 2: point codes [17,8], line codes [17,9] with Frobenius x^4 twist),
transported to the Singer labelling used here (point i ~ line N+(i+d)%N, d in D).  Singer point i <-> projective point
alpha^{-i} in F_4096 = F_16^3 (basis 1, alpha, alpha^2), lines matched by point sets.  Returns per-vertex local generator
rows (one k_v-vector per incident edge, in nb order) as uint8 arrays over pg.py's GF(16) tables."""
import numpy as np, pickle, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from singer import singer, gf_mul
from pg import build
def make():
    fn=os.path.join(os.path.dirname(os.path.abspath(__file__)),'code16.pkl')
    if os.path.exists(fn): return pickle.load(open(fn,'rb'))
    N,D,pw=singer(); G=build(2,4,8,9,2); F=G['F']; T=G['T']; idx={t:i for i,t in enumerate(T)}
    # F16 inside F4096: find beta in subfield with beta^4+beta+1=0
    sub=[0]+[pw[(273*k)%4095] for k in range(15)]
    def pw_(y,n):
        r=1
        for _ in range(n): r=gf_mul(r,y)
        return r
    beta=[b for b in sub if b and (pw_(b,4)^b^1)==0][0]
    phi=[0]*16
    for e in range(16):
        v=0
        for i in range(4):
            if e>>i&1: v^=pw_(beta,i)
        phi[e]=v
    assert len(set(phi))==16
    # check phi is a field hom
    for a in range(16):
        for b in range(16): assert phi[int(F.mul[a,b])]==gf_mul(phi[a],phi[b]) and phi[int(F.add[a,b])]==phi[a]^phi[b]
    al=pw[1]; al2=pw[2]; coords={}
    for c0 in range(16):
        for c1 in range(16):
            for c2 in range(16): coords[phi[c0]^gf_mul(phi[c1],al)^gf_mul(phi[c2],al2)]=(c0,c1,c2)
    assert len(coords)==4096
    def norm(v):
        f=[x for x in v if x][0]; inv=int(F.inv[f]); return tuple(int(F.mul[inv,x]) for x in v)
    pt=[idx[norm(coords[pw[(-i)%4095]])] for i in range(N)]   # Singer point i -> pg point index
    assert len(set(pt))==N
    ptsOn={}   # pg line L -> frozenset of pg points
    for pp in range(N):
        for L in G['nb'][pp]: ptsOn.setdefault(L,set()).add(pp)
    inv_={frozenset(s):L for L,s in ptsOn.items()}
    ln=[inv_[frozenset(pt[(l-d)%N] for d in D)] for l in range(N)]
    assert len(set(ln))==N
    nb=[[] for _ in range(2*N)]
    for i in range(N):
        for d in D: l=(i+d)%N; nb[i].append(N+l); nb[N+l].append(i)
    rows=[]
    for v in range(2*N):
        if v<N: rows.append(np.array([G['ptc'][(pt[v],ln[u-N])] for u in nb[v]],dtype=np.uint8))
        else: rows.append(np.array([G['lnc'][(pt[u],ln[v-N])] for u in nb[v]],dtype=np.uint8))
    out=dict(N=N,D=D,nb=nb,rows=rows,add=F.add.astype(np.uint8),mul=F.mul.astype(np.uint8),inv=F.inv.astype(np.uint8))
    pickle.dump(out,open(fn,'wb')); return out
if __name__=="__main__":
    C=make(); print('built', C['N'], [r.shape for r in C['rows'][:1]], [r.shape for r in C['rows'][-1:]])
