import sys, time, random
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from tanner import *
H=build(4); NP,NL=incidence(H)
adj=graph(H)
for mseed in [0,1,2]:
    M=matching(H,random.Random(mseed))
    if mseed==0:
        g=girth(graph(H,M),NP+NL); print("H(4)-M girth",g,flush=True)
    for (p,e) in [(2,2),(2,4)]:
        K=GF(p,e)
        for sp,sl in [(0,0),(0,1),(1,0),(1,1)]:
            if K.q==4 and sp==1 and sl==1: continue
            t=time.time(); T=Tanner(H,K,drop=M,sp=sp,sl=sl,mode='nat')
            S=T.dimW()
            print(f"q=4-M m{mseed} K=F{K.q} sp={sp} sl={sl} mds={T.mds_ok()} E={T.E} S={S} S/E={S/T.E:.4f} a={a_of(S,T.E):.3f} ({time.time()-t:.1f}s)",flush=True)
        if mseed==0:
            for mode in ['randgauge','randgate']:
                T=Tanner(H,K,drop=M,mode=mode,seed=0); S=T.dimW()
                print(f"q=4-M m{mseed} K=F{K.q} {mode} E={T.E} S={S} S/E={S/T.E:.4f}",flush=True)
