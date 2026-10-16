"""Bigraded Hilbert function h(7,7) of F2[P0..2, L0..2]/(relations) for the twisted incidence:
twist exponent pair (a,b): relations sum_i P_i^a L_i^b for all (a,b)=(2^s,2^t) with s-t = c mod 4 and both exponents <= 7.
Rank over F2 (field independent)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir
import itertools, numpy as np
def monos(k): return [m for m in itertools.product(range(k+1),repeat=3) if sum(m)==k]
def h(c,DP=7,DL=7):
    rels=[(2**s,2**t) for s in range(3) for t in range(3) if (s-t)%4==c]
    tgt={(p,l):i for i,(p,l) in enumerate(itertools.product(monos(DP),monos(DL)))}
    rows=[]
    for (a,b) in rels:
        if a>DP or b>DL: continue
        for p in monos(DP-a):
            for l in monos(DL-b):
                r=np.zeros(len(tgt),dtype=np.uint8)
                for i in range(3):
                    pp=list(p); pp[i]+=a; ll=list(l); ll[i]+=b
                    r[tgt[(tuple(pp),tuple(ll))]]^=1
                rows.append(r)
    if not rows: return len(tgt),rels
    M=np.array(rows); rk=0
    for col in range(M.shape[1]):
        nz=np.nonzero(M[rk:,col])[0]
        if len(nz)==0: continue
        i=rk+nz[0]; M[[rk,i]]=M[[i,rk]]
        rr=np.nonzero(M[:,col])[0]; rr=rr[rr!=rk]; M[rr]^=M[rk]; rk+=1
        if rk==M.shape[0]: break
    return len(tgt)-rk,rels
vals=[]
for c,name in [(0,'untwisted'),(1,'x^2 / x^8'),(2,'x^4 (Hermitian)'),(3,'x^8 / x^2')]:
    r=h(c); vals.append(r[0]); print(name,r)
print('PASS_EXPQE_HILB' if vals==[512,708,772,708] else 'FAIL_EXPQE_HILB')   # = measured dim W of each twist (expQE table)
