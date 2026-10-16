"""lambda_max(M|C) at extreme (reducible: each vertex = direct sum of two 2x2 blocks along a pairing) degree-4 configurations.
M is bilinear in (O_P, O_L) and each O_v ranges in a simplex whose vertices are the 3 pairings (x 1/4), so the sup of
lambda_max over flat degree-4 codes is attained at such configurations (convexity in each side separately)."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from xk_sparse import load_any
name = sys.argv[1]; V, nP, E, g = load_any(name); rng = np.random.default_rng(0)
N = np.zeros((len(V), E))
for i, es in enumerate(V): N[i, es] = 1
U, s, Vt = np.linalg.svd(N); B = Vt[np.sum(s > 1e-9):].T
PAIRS = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]
def Omat(choice, side):
    O = np.zeros((E, E)); rng_ = range(nP) if side == 0 else range(nP, len(V))
    for i in rng_:
        es = V[i]
        for (x, y) in PAIRS[choice[i]]: O[es[x], es[y]] = O[es[y], es[x]] = 0.25
    return O
def lam(choice):
    OP, OL = Omat(choice, 0), Omat(choice, 1)
    M = -3/16*np.eye(E) + 0.5*(OP+OL) - 1.5*(OP@OL + OL@OP)
    return np.linalg.eigvalsh(B.T @ M @ B)[-1]
best = -9
for trial in range(int(sys.argv[2]) if len(sys.argv) > 2 else 20):
    ch = rng.integers(3, size=len(V)); cur = lam(ch); improved = True
    while improved:
        improved = False
        for i in rng.permutation(len(V)):
            for c in range(3):
                if c == ch[i]: continue
                old = ch[i]; ch[i] = c; v = lam(ch)
                if v > cur + 1e-12: cur = v; improved = True
                else: ch[i] = old
    best = max(best, cur); print(name, 'trial', trial, 'local-max lambda_max(M|C) over pairings =', round(cur, 6), flush=True)
print(name, 'BEST', best)
