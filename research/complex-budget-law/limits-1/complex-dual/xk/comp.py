"""Complementary gauge: minimiser of the convex two-sided capacity Phi(x) = sum_v logdet(Q_v^H e^{2x} Q_v) - 2 sum_e x_e.
Stationarity: a_e + b_e = 1 on every edge (E equations, E unknowns: NOT overdetermined)."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from xk_sparse import *
from scipy.optimize import minimize

def phi(S, Qs, x):
    val = -2 * x.sum(); g = -2 * np.ones(S.E)
    for es, Q in zip(S.V, Qs):
        Y = np.exp(x[es])[:, None] * Q
        Gm = Y.conj().T @ Y; val += np.real(np.linalg.slogdet(Gm)[1])
        lev = np.real(np.sum((Y @ np.linalg.inv(Gm)) * Y.conj(), axis=1))
        g[es] += 2 * lev
    return val, g

def comp_gauge(S, Qs, x0=None, tol=1e-12):
    x0 = np.zeros(S.E) if x0 is None else x0
    r = minimize(lambda x: phi(S, Qs, x), x0, jac=True, method='L-BFGS-B', options=dict(maxiter=5000, gtol=tol, ftol=1e-15))
    x = r.x - r.x.mean(); a, b = leverages(S, Qs, x)
    return x, float(abs(a + b - 1).max())
