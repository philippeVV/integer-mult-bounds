import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from comp import *
from local import local_moments
from scipy.optimize import minimize
name, f = sys.argv[1], sys.argv[2]
V, nP, E, g = load_any(name); S = Struct(V, nP, E)
Q = [orth(np.asarray(q, complex)) for q in np.load(f, allow_pickle=True)]
xs = np.load(f.replace('.npy', '_gauge.npy')) if os.path.exists(f.replace('.npy', '_gauge.npy')) else np.zeros(E)
xc, err = comp_gauge(S, Q); a, b = leverages(S, Q, xc)
m = local_moments(S, Q, xc)[0]
print('comp gauge: err', err, 'theta3', theta(S, Q, 3, xc), 'F_1/2 = 4m3-4m2+m1 =', 4*m[3]-4*m[2]+m[1], 'lev range', a.min(), a.max())
rng = np.random.default_rng(0)
starts = [xs, xc, np.zeros(E), flat_gauge(S, Q, 0)[0], flat_gauge(S, Q, 1)[0]] + [xc + 0.3 * rng.standard_normal(E) for _ in range(3)]
best = None
for st in starts:
    r = minimize(lambda x: gauge_value_grad(S, Q, 3, x), st, jac=True, method='L-BFGS-B', options=dict(maxiter=20000, gtol=1e-13, ftol=1e-16))
    a, b = leverages(S, Q, r.x)
    print('  start -> %.10f  |g| %.1e  xspread %.2f  lev [%.3f,%.3f]' % (r.fun, np.abs(r.jac).max(), r.x.max() - r.x.min(), a.min(), a.max()), flush=True)
    if best is None or r.fun < best.fun: best = r
print(name, 'inf over gauges (unconstrained, multistart):', best.fun, ' 1/8 =', 0.125)
