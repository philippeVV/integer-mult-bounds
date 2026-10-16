"""Exact coupled homogeneous test at n=6 for the profile (1/3,1/3,1/3,2/3,2/3,2/3), started from the split-LP worst P."""
import numpy as np, sys
from scipy.optimize import minimize
exec(open('coupled.py').read().split('rng = np.random.default_rng(0)')[0])
rng = np.random.default_rng(3); n = 6
best = -9
for pi in [np.array([3, 4, 5, 0, 1, 2]), np.array([4, 5, 3, 2, 0, 1])]:
    for s in range(60):
        z0 = rng.standard_normal(n * n)
        r = minimize(lambda z: -G(proj(z, n), pi)[0] + 300 * G(proj(z, n), pi)[1], z0, method='BFGS', options=dict(maxiter=600))
        P = proj(r.x, n); gv, pen = G(P, pi)
        if pen < 1e-6 and gv > best: best = gv; Pb = P
print('n=6 coupled homogeneous: max excess per coordinate', best, 'lev', np.round(np.real(np.diag(Pb)), 4))
