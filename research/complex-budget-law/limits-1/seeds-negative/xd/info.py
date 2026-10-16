import time
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from tanner import *
for q,cfgs in [(2,[(1,2),(2,1)]),(3,[(1,3),(3,1),(2,2)]),(4,[(2,3),(3,2)])]:
    H=build(q)
    for p,e in ([(2,1),(2,2)] if q==2 else [(3,1),(3,2)] if q==3 else [(2,2),(2,4)]):
        K=GF(p,e)
        for kP,kL in cfgs:
            for sp,sl in ([(0,0),(0,1)] if (K.q>q) or q==4 else [(0,0)]):
                T=Tanner(H,K,sp=sp,sl=sl,kP=kP,kL=kL); S=T.dimW()
                print(f"full H({q}) K=F{K.q} [deg={T.deg}] kP={kP} kL={kL} sp={sp} sl={sl} E={T.E} S={S} S/E={S/T.E:.4f} a={a_of(S,T.E):.3f}",flush=True)
