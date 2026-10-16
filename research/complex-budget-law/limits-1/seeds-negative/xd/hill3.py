"""H(3) [4,2] gates over F9: hill-climb on per-vertex MDS gates starting from the natural code (S=128)."""
import sys, time, random
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from tanner import *
seed=int(sys.argv[1]); iters=int(sys.argv[2]); K=GF(3,2)
H=build(3); T=Tanner(H,K); rng=random.Random(seed)
S=T.dimW(); best=S; print("start",S,flush=True); t0=time.time()
for it in range(iters):
    side=rng.random()<0.5; v=rng.randrange(T.NP)
    if rng.random()<0.5:
        G=randmds(K,2,4,rng)
    else:
        G=(T.GP[v] if side else nullspace(K,T.HL[v])).copy(); c=rng.randrange(4); G[:,c]=K.mul[rng.randrange(1,K.q),G[:,c]]
    if side: old=T.GP[v]; T.GP[v]=G
    else: old=T.HL[v]; T.HL[v]=nullspace(K,G)
    S=T.dimW()
    if S>=best: best=S
    else:
        if side: T.GP[v]=old
        else: T.HL[v]=old
    if it%100==0: print(it,best,f"{time.time()-t0:.0f}s",flush=True)
print("final",best,best/T.E)
