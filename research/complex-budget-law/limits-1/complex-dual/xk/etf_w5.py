import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from xk_sparse import *
from gauge_hess import XY_from_codes, lam_min
from scipy.optimize import minimize
V, nP, E, g = load_any('W5'); S = Struct(V, nP, E); rng = np.random.default_rng(0)
Qs = [orth(np.asarray(q, complex)) for q in np.load('init_W5_etf_n0_s0.npy', allow_pickle=True)]
X, Y = XY_from_codes(V, nP, E, Qs); print('W5 ICO lambda_min(Gam) =', lam_min(X, Y)[0], ' Theta(0) =', theta(S, Qs, 3), flush=True)
for st in [np.zeros(E), flat_gauge(S, Qs, 0)[0]] + [rng.standard_normal(E) * s for s in (0.3, 1, 3)]:
    r = minimize(lambda x: gauge_value_grad(S, Qs, 3, x), st, jac=True, method='L-BFGS-B', options=dict(maxiter=3000, gtol=1e-12))
    print('  inf from start: %.9f spread %.2f' % (r.fun, np.ptp(r.x)), flush=True)
