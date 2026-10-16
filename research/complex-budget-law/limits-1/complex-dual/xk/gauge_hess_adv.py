"""Adversary over flat degree-4 codes for lambda_min(Gam) (want: always < 0).
For a flat rank-2 projector on C^4, X_v = 4(|P|^{o2} - I/4) = w1 M1 + w2 M2 + w3 M3 (the 3 perfect matchings), w in simplex
(checked numerically below); we search the whole simplex product (a relaxation if some w are not realisable)."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from gauge_hess import *
from scipy.optimize import minimize
rng = np.random.default_rng(1)
# realisability check of the matching structure
for _ in range(5):
    q = np.linalg.qr(rng.standard_normal((4, 2)) + 1j * rng.standard_normal((4, 2)))[0]
    d = np.zeros(4)
    for _ in range(3000):
        Y_ = orth(np.exp(d)[:, None] * q); a = np.sum(abs(Y_)**2, 1); d += 0.5*(np.log(.5)-np.log(a))
    P = orth(np.exp(d)[:, None] * q); O = 4 * abs(P @ P.conj().T)**2
    print('flat 4x4: |P01|^2=|P23|^2?', np.isclose(O[0,1], O[2,3]), ' |P02|=|P13|?', np.isclose(O[0,2], O[1,3]), ' |P03|=|P12|?', np.isclose(O[0,3], O[1,2]), ' sum', round(O[0,1]+O[0,2]+O[0,3], 9))
name = sys.argv[1]; V, nP, E, g = load_any(name)
PAIRS = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]
Ms = []
for i, es in enumerate(V):
    mats = []
    for pr in PAIRS:
        m = np.zeros((E, E))
        for (x, y) in pr: m[es[x], es[y]] = m[es[y], es[x]] = 1
        mats.append(m)
    Ms.append(mats)
def build(W):
    X = np.zeros((E, E)); Y = np.zeros((E, E))
    for i in range(len(V)):
        T = sum(W[i, k] * Ms[i][k] for k in range(3))
        if i < nP: X += T
        else: Y += T
    return X, Y
def f(theta):
    W = np.exp(theta.reshape(-1, 3)); W /= W.sum(1, keepdims=True); X, Y = build(W); return lam_min(X, Y)[0], W
best = -9
for trial in range(int(sys.argv[2]) if len(sys.argv) > 2 else 6):
    th = rng.standard_normal(3 * len(V)) * (3 if trial % 2 else 0.5)
    r = minimize(lambda t: -f(t)[0], th, method='Powell' if len(th) < 0 else 'L-BFGS-B', options=dict(maxiter=60))
    val, W = f(r.x); best = max(best, val)
    print(f'{name} trial {trial}: max lambda_min(Gam) found = {val:.6f}  (max weight per vertex: median {np.median(W.max(1)):.3f})', flush=True)
# extremes: random pairings + greedy
for trial in range(3):
    ch = rng.integers(3, size=len(V)); W = np.eye(3)[ch]; cur = lam_min(*build(W))[0]; imp = True
    while imp:
        imp = False
        for i in rng.permutation(len(V)):
            for c in range(3):
                if c == ch[i]: continue
                old = ch[i]; ch[i] = c; W = np.eye(3)[ch]; v = lam_min(*build(W))[0]
                if v > cur + 1e-12: cur = v; imp = True
                else: ch[i] = old
    best = max(best, cur); print(f'{name} reducible extreme (greedy) lambda_min = {cur:.6f}', flush=True)
print(name, 'BEST max lambda_min(Gam) =', best)
