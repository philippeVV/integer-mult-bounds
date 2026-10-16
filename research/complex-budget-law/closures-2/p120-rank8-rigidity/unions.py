"""Unions of Weyl-conjugate sub-root-system patterns T = T_R1 u T_R2 (u T_R3), d <= 30.
Certified lower bound: alpha(T) (CP-SAT; uses proven upper bound on alpha only if optimal) -> identity principal minor."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, sys, json, itertools
from e8 import R, G, F, n
from subsys import SUBS, Tgraph, key
from alpha import alpha
def weyl(rng, L=40):
    """random Weyl element as a permutation of P120 (signs dropped)"""
    perm=np.arange(n)
    for _ in range(L):
        s=R[rng.integers(n)]
        # reflection: v -> v - (v.s/4)*s  (doubled coords, |s|^2 = 8 -> 2 v.s/|s|^2 = v.s/4)
        perm=np.array([key(R[p]-(int(R[p]@s)//4)*s) for p in perm])
    return perm
def conj(Rs,perm): return sorted({int(perm[i]) for i in Rs})
if __name__=='__main__':
    rng=np.random.default_rng(int(sys.argv[1]) if len(sys.argv)>1 else 0)
    combos=[('A4^2','A4^2'),('D4^2','D4^2'),('A4^2','D4^2'),('A8','A4^2'),('A8','D4^2'),('A8','A8'),('A4^2','A4^2','A4^2'),('D4^2','D4^2','A4^2'),('E7+A1','A4^2'),('D8','A4^2')]
    out=[]
    for combo in combos:
        for rep in range(int(sys.argv[2]) if len(sys.argv)>2 else 4):
            T=np.zeros((n,n),bool); Rsets=[]
            for k,name in enumerate(combo):
                Rs=SUBS[name] if k==0 else conj(SUBS[name],weyl(rng))
                Rsets.append(Rs); T|=Tgraph(Rs)
            d=T.sum(1).mean()
            if d>30: print(combo,'d=%.1f >30 skip'%d,flush=True); continue
            lo,hi,opt=alpha(T,list(range(n)))
            print(combo,rep,'d=%.2f maxdeg %d alpha %d (bound %d, optimal %s)'%(d,T.sum(1).max(),lo,hi,opt),flush=True)
            out.append(dict(combo=combo,rep=rep,d=float(d),alpha=lo,alpha_ub=hi,R=Rsets))
    json.dump(out,open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'unions_%s.json'%(sys.argv[1] if len(sys.argv)>1 else '0')),'w'))
    ok = all(r['alpha'] == r['alpha_ub'] for r in out) and min(r['alpha'] for r in out) >= 15
    print('min certified alpha over unions with d <= 30:', min(r['alpha'] for r in out))
    print('PASS_UNIONS' if ok else 'FAIL_UNIONS')
