"""soundness (info-set split) + routing distance bound, expQV/sound_m style, natural code on H(q) (q prime)."""
import sys, time, random, collections, numpy as np
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from tanner import *
from rankp import rank_modp
q=int(sys.argv[1]); seed=int(sys.argv[2])
H=build(q); K=GF(q,1); T=Tanner(H,K); k=T.kP
Mx=T.check_matrix()
Kern=nullspace(K,Mx) if q==3 else None
if Kern is None:
    # nullspace via flint-free path: small enough for q=5? use numba rref on int16 (11718^2) -- ok in time
    Kern=nullspace(K,Mx)
S=Kern.shape[0]; E=T.E
edges=[(p,l) for p in range(T.NP) for l in T.pcols[p]]
Wb=np.zeros((S,E),dtype=np.int16)
for ei,(p,l) in enumerate(edges):
    g=T.GP[p][:,T.pidx[(p,l)]]; acc=np.zeros(S,dtype=np.int16)
    for r in range(k):
        if g[r]: acc=K.add[acc,K.mul[g[r],Kern[:,p*k+r]]]
    Wb[:,ei]=acc
rng=random.Random(seed); perm=list(range(E)); rng.shuffle(perm)
_,piv=rref(K,Wb[:,perm]); assert len(piv)==S
R=sorted(perm[c] for c in piv); split={edges[ei]:t for t,ei in enumerate(R)}
nunk=T.NP*k+S; rows=np.zeros((Mx.shape[0],nunk),dtype=np.int16); r=0
for l in range(T.NL):
    for h in T.HL[l]:
        for t,p in enumerate(T.lcols[l]):
            if h[t]==0: continue
            if (p,l) in split: c=T.NP*k+split[(p,l)]; rows[r,c]=K.add[rows[r,c],h[t]]
            else:
                col=T.GP[p][:,T.pidx[(p,l)]]; rows[r,p*k:(p+1)*k]=K.add[rows[r,p*k:(p+1)*k],K.mul[h[t],col]]
        r+=1
sd=nunk-rank_modp(rows,q)
Rs=set(R); adj=collections.defaultdict(list)
for ei,(p,l) in enumerate(edges):
    if ei in Rs: continue
    adj[p].append(T.NP+l); adj[T.NP+l].append(p)
def dist(a,b):
    d={a:0}; dq=collections.deque([a])
    while dq:
        x=dq.popleft()
        if x==b: return d[x]
        for y in adj[x]:
            if y not in d: d[y]=d[x]+1; dq.append(y)
    return 10**9
ds=[dist(edges[ei][0],T.NP+edges[ei][1]) for ei in R]; Ep=E-S; tot=sum(ds)
print(f"H({q}) nat seed={seed}: E={E} S={S} split-dim={sd} sound={sd==S} dist={dict(collections.Counter(ds))} routing<= {Ep}/{tot} delta_nd>={1-Ep/tot:.4f}",flush=True)
A=np.zeros((T.NP+T.NL,)*2,dtype=np.int8)
for p,l in edges: A[p,T.NP+l]=A[T.NP+l,p]=1
np.save(f"A{q}.npy",A)
