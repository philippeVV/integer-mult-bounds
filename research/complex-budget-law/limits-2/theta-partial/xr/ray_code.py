"""Decompose a saved code into (flat code C0, relative gauge delta) and scan Phi(t) = inf_x Theta3(C0; x, x + t delta), t in [0, 1.5]; t = 1 is the code itself.
usage: ray_code.py GRAPH codes.npy"""
import sys, numpy as np
from model import *
from xk_sparse import flat_gauge, gauge_value_grad
from scipy.optimize import minimize
def inner_min(S, Qs, J, starts, maxiter=20000):
    best = None
    for st in starts:
        r = minimize(lambda z: gauge_value_grad(S, Qs, J, z), st, jac=True, method='L-BFGS-B', options=dict(maxiter=maxiter, gtol=1e-13, ftol=1e-16))
        if best is None or r.fun < best.fun: best = r
    return float(best.fun), best.x
name, f = sys.argv[1], sys.argv[2]
V, nP, E, g = load_any(name); S = Struct(V, nP, E)
Q = [orth(np.asarray(q, complex)) for q in np.load(f, allow_pickle=True)]
Q0 = []; fv = np.zeros((2, E))
for i, (es, q) in enumerate(zip(V, Q)):
    d = np.zeros(len(es))
    for _ in range(6000):
        Yq = orth(np.exp(d)[:, None] * q); a = np.sum(abs(Yq)**2, 1); d += 0.5*(np.log(.5) - np.log(a))
    Q0.append(orth(np.exp(d)[:, None] * q)); fv[0 if i < nP else 1, es] = -d
delta = fv[1] - fv[0]
print(f'{name} {f}: max flat err {max(abs(np.sum(abs(q)**2,1)-.5).max() for q in Q0):.1e}, |delta|_rms {np.sqrt(np.mean(delta**2)):.3f}', flush=True)
x = np.zeros(E); vals = []
for t in ([0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0, 1.1, 1.25, 1.5] if E < 500 else [0.25, 0.5, 0.75, 1.0, 1.25]):
    Qs = [Q0[i] * (np.exp(t * delta[es])[:, None] if i >= nP else 1) for i, es in enumerate(V)]
    st = [x, np.zeros(E), -t * delta / 2, flat_gauge(S, Qs, 0, 300)[0], flat_gauge(S, Qs, 1, 300)[0]] if E < 500 else ([x, -t * delta / 2] if t > 0.3 else [-t * delta / 2])
    val, x = inner_min(S, Qs, 3, st, maxiter=5000); vals.append(val); print(f'  t={t}: Phi={val:.8f}', flush=True)
print('monotone' if all(vals[i+1] <= vals[i] + 1e-7 for i in range(len(vals)-1)) else 'NON-MONOTONE', 'max', max(vals))
