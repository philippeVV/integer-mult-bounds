"""Numerical evidence only: smallest d with X Y^T = I off Gamma_E (diag 1, zero on non-edges), by ALS."""
import json, sys, numpy as np
from gamma import edges, n
def als(Z, d, iters=400, seed=0):
    rng = np.random.default_rng(seed)
    X = rng.standard_normal((n, d)); Y = rng.standard_normal((n, d))
    for it in range(iters):
        for a in range(n):
            idx = [a] + Z[a]; A = Y[idx]; b = np.zeros(len(idx)); b[0] = 1
            X[a] = np.linalg.lstsq(A, b, rcond=None)[0]
        for c in range(n):
            idx = [c] + Zt[c]; A = X[idx]; b = np.zeros(len(idx)); b[0] = 1
            Y[c] = np.linalg.lstsq(A, b, rcond=None)[0]
    M = X @ Y.T; err = 0.0
    for a in range(n):
        err = max(err, abs(M[a, a] - 1), max((abs(M[a, b]) for b in Z[a]), default=0))
    return err
f = sys.argv[1]; data = json.load(open(f)); r = (data[0] if isinstance(data, list) else data)['r']
E = set(edges(r)); E |= {(j, i) for i, j in E}
Z = [[b for b in range(n) if b != a and (a, b) not in E] for a in range(n)]
Zt = Z  # symmetric pattern
for d in map(int, sys.argv[2:]):
    print(f, 'd', d, 'max residual', '%.2e' % als(Z, d), flush=True)
