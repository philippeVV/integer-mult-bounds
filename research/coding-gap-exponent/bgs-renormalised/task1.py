import sys; sys.path.insert(0, __file__.rsplit('/',1)[0])
from bh import CORE, SESS
from cp import solve
import numpy as np
core_s = [(s, t) for q, s, t in SESS]
E = list(CORE); nv = 182; ext_s = []
for q, s, t in SESS:
    E += [(s, nv), (t, nv+1)]; ext_s.append((nv, nv+1)); nv += 2
for name, (n, EE, ss) in {"core": (182, list(CORE), core_s), "full leaf-ext": (nv, E, ext_s)}.items():
    res, w, y, d, _ = solve(n, EE, ss, verbose=True)
    z = res.fun
    print(f"{name}: m={len(EE)}  MCF* = z* = {z:.13f}  (147/157={147/157:.13f})  min dist-sum check {d.sum():.12f}")
    print(f"   HiGHS vertex dual: min w = {w.min():.3e}  #w<1e-12: {(w<1e-12).sum()}", end='')
    if len(EE) > 735: print(f"  max leaf weight {w[735:].max():.3e}")
    else: print()
    res2, w2, y2, d2, om = solve(n, EE, ss, mode='maxmin', zcap=z*(1+1e-10), verbose=True, method='highs-ipm', tol=1e-5)
    print(f"   max-min weight on optimal face: omega* = {om:.6e}   (uniform 1/785 = {1/785:.6e}); weights spread [{w2.min():.4e},{w2.max():.4e}]")
