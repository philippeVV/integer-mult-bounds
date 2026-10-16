"""Tools for Conjecture M / linear-rule tests.
f(u, v) = Theta_3 with point codes e^u A0_p and line codes e^v B0_l (C0 flat). Split gradients g_u, g_v (log-scaling)."""
import numpy as np
from model import *
from xk_sparse import theta_grad
from scipy.optimize import minimize

def fuv(S, Q0, u, v):
    Ys = []
    for i, (es, q) in enumerate(zip(S.V, Q0)):
        s = (u if i < S.nP else v)[es]; Ys.append(orth(np.exp(s)[:, None] * q))
    th, m, Gb = theta_grad(S, Ys, 3)
    gu = np.zeros(S.E); gv = np.zeros(S.E)
    for i, (es, Y, G) in enumerate(zip(S.V, Ys, Gb)):
        Pi = Y @ Y.conj().T; h = np.real(np.diag(G @ Pi + Pi @ G - 2 * Pi @ G @ Pi))
        (gu if i < S.nP else gv)[es] += h
    return th, gu, gv

def inner(S, Q0, d, starts, gtol=1e-13, maxiter=30000):
    """precise inf_x f(x, x+d); returns value, x, gu, gv"""
    best = None
    for st in starts:
        def fg(x):
            th, gu, gv = fuv(S, Q0, x, x + d); return th, gu + gv
        r = minimize(fg, st, jac=True, method='L-BFGS-B', options=dict(maxiter=maxiter, gtol=gtol, ftol=1e-17, maxcor=30))
        if best is None or r.fun < best.fun: best = r
    x = best.x; th, gu, gv = fuv(S, Q0, x, x + d)
    return th, x, gu, gv

def setup(name, base, seed):
    V, nP, E, g = load_any(name); S = Struct(V, nP, E); rng = np.random.default_rng(seed)
    if base == 'sic': W = np.full((len(V), 3), 1 / 3)
    elif base.startswith('bdry'):   # bdry<frac>: point at fraction frac of the way from SIC to the Gamma-boundary along a random Dirichlet(0.7) w
        frac = float(base[4:]); W0 = np.full((len(V), 3), 1 / 3); W1 = rng.dirichlet([0.7] * 3, size=len(V))
        o = np.ones(E) / np.sqrt(E); P = np.eye(E) - np.outer(o, o)
        def lam(Wt):
            X, Y = XYmats(V, nP, E, Wt); G = xdelta(Kmat(X, Y))[0]; return np.linalg.eigvalsh(P @ G @ P + np.outer(o, o))[0]
        lo, hi = 0.0, 1.0
        for _ in range(30):
            mid = (lo + hi) / 2
            if lam((1 - mid) * W0 + mid * W1) > 0: lo = mid
            else: hi = mid
        W = (1 - frac * lo) * W0 + frac * lo * W1
    else: W = rng.dirichlet([float(base)] * 3, size=len(V))
    Q0 = [flat_code(W[i], rng) for i in range(len(V))]
    X, Y = XYmats(V, nP, E, W); K = Kmat(X, Y); G, C, D = xdelta(K)
    o = np.ones(E) / np.sqrt(E); P = np.eye(E) - np.outer(o, o)
    lmin = np.linalg.eigvalsh(P @ G @ P + np.outer(o, o))[0]
    return S, Q0, W, G, C, D, lmin, rng
