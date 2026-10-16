import sys, time
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from tanner import *
from rankp import rank_modp
q=int(sys.argv[1]); t=time.time()
H=build(q); NP,NL=incidence(H); print("built",NP,NL,f"{time.time()-t:.0f}s",flush=True)
if len(sys.argv)>2: print("girth",girth(graph(H),NP+NL),flush=True)
K=GF(q,1); T=Tanner(H,K); print("tanner",T.deg,T.kP,"mds",T.mds_ok(),f"{time.time()-t:.0f}s",flush=True)
Mx=T.check_matrix(); print("matrix",Mx.shape,f"{time.time()-t:.0f}s",flush=True)
r=rank_modp(Mx,q); S=Mx.shape[1]-r
print(f"H({q}) natural [{T.deg},{T.kP}] E={T.E} S={S} S/E={S/T.E:.4f} a={a_of(S,T.E):.3f} ({time.time()-t:.0f}s)",flush=True)
assert q != 5 or S == 2263, S   # limits-1 replay: H(5) natural [6,3] over F5
print('PASS_EXPXD_H%d' % q)
