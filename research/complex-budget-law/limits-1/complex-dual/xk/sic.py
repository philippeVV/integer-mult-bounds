"""SIC (tetrahedral) local codes: rows = qubit SIC vectors / sqrt2  -> flat, equiangular |P_ij|^2 = 1/12, MDS.
Test: gauge Hessian Gam (positive definite?), global inf over gauges, and the code/gauge Schur complement."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from maxmin_hess import *
from gauge_hess import XY_from_codes, lam_min
from scipy.optimize import minimize
w = np.exp(2j * np.pi / 3)
SIC = np.array([[1, 0], [1/np.sqrt(3), np.sqrt(2/3)], [1/np.sqrt(3), np.sqrt(2/3) * w], [1/np.sqrt(3), np.sqrt(2/3) * w**2]]) / np.sqrt(2)
def sic_codes(V, rng, randomize=True):
    out = []
    for es in V:
        Y = SIC.copy()
        if randomize:
            Y = Y[rng.permutation(4)] * np.exp(2j * np.pi * rng.random(4))[:, None]
        out.append(orth(Y))
    return out
if __name__ == '__main__':
    name = sys.argv[1]; V, nP, E, g = load_any(name); S = Struct(V, nP, E); rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 0)
    Qs = sic_codes(V, rng)
    P = Qs[0] @ Qs[0].conj().T; print('SIC |P|^2 =', np.round(abs(P)**2, 4).tolist())
    X, Y = XY_from_codes(V, nP, E, Qs); l, _ = lam_min(X, Y); print(name, 'SIC lambda_min(Gam) =', l, flush=True)
    print('Theta_3 at gauge 0 =', theta(S, Qs, 3))
    starts = [np.zeros(E), flat_gauge(S, Qs, 0)[0]] + [rng.standard_normal(E) * s for s in (0.3, 1, 2, 3, 3)]
    for st in starts:
        r = minimize(lambda x: gauge_value_grad(S, Qs, 3, x), st, jac=True, method='L-BFGS-B', options=dict(maxiter=5000, gtol=1e-12))
        print('  inf from start: %.9f  spread %.2f' % (r.fun, np.ptp(r.x)), flush=True)
    if E <= 200 and len(sys.argv) < 4:
        test(S, Qs, name + ' SIC')
