import sys, time
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from tanner import *
H=build(3); NP,NL=incidence(H)
for (p,e) in [(3,1),(3,2),(3,3)]:
    K=GF(p,e)
    for mode in ['nat','randgauge','randgate']:
        for dP,dL in ([(0,0),(1,0),(0,1),(1,1)] if mode=='nat' else [(0,0)]):
            for seed in ([0] if mode=='nat' else [0,1,2]):
                if mode!='nat' and K.q==3 and mode=='randgate': pass
                t=time.time(); T=Tanner(H,K,mode=mode,seed=seed,dualP=bool(dP),dualL=bool(dL))
                S=T.dimW()
                print(f"q=3 K=F{K.q} {mode} dualP={dP} dualL={dL} seed={seed} mds={T.mds_ok()} E={T.E} S={S} S/E={S/T.E:.4f} a={a_of(S,T.E):.3f} ({time.time()-t:.1f}s)",flush=True)
