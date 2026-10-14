"""Levenberg-Marquardt for masked low rank: N = X Y^T (or sym V V^T), N_aa=1, N_ab=0 on forbidden pairs.
usage: lm.py mode K seed [sym]   mode: line | file.json (plane r)"""
import sys, json, numpy as np
from scipy.optimize import least_squares
from scipy.sparse import lil_matrix
from common import *
mode, K, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]); sym = len(sys.argv) > 4 and sys.argv[4] == 'sym'
if mode == 'line':
    F = (D[Pa][:, Pa] == 0); np.fill_diagonal(F, False)
else:
    d = json.load(open(mode)); d = d[0] if isinstance(d, list) else d; F = orthmask(d['r'])
rows = [(a, a) for a in range(n)] + [(a, b) for a in range(n) for b in range(n) if a != b and not F[a, b] and (not sym or a < b)]
ra = np.array([p[0] for p in rows]); rb = np.array([p[1] for p in rows]); tgt = (ra == rb).astype(float)
def unpack(z):
    if sym: V = z.reshape(n, K); return V, V
    return z[:n * K].reshape(n, K), z[n * K:].reshape(n, K)
def fun(z):
    X, Y = unpack(z); return (X[ra] * Y[rb]).sum(1) - tgt
def jac(z):
    X, Y = unpack(z); m = len(rows)
    J = lil_matrix((m, z.size))
    from scipy.sparse import coo_matrix
    I = np.repeat(np.arange(m), K); kk = np.tile(np.arange(K), m)
    if sym:
        c1 = ra[:, None] * K + np.arange(K); c2 = rb[:, None] * K + np.arange(K)
        v1 = Y[rb]; v2 = X[ra]
        return coo_matrix((np.concatenate([v1.ravel(), v2.ravel()]), (np.concatenate([I, I]), np.concatenate([c1.ravel(), c2.ravel()]))), shape=(m, z.size)).tocsr()
    c1 = ra[:, None] * K + np.arange(K); c2 = n * K + rb[:, None] * K + np.arange(K)
    return coo_matrix((np.concatenate([Y[rb].ravel(), X[ra].ravel()]), (np.concatenate([I, I]), np.concatenate([c1.ravel(), c2.ravel()]))), shape=(m, z.size)).tocsr()
best = None
for rs in range(int(sys.argv[5]) if len(sys.argv) > 5 else 3):
    rng = np.random.default_rng(seed * 100 + rs)
    z0 = rng.standard_normal((1 if sym else 2) * n * K) / np.sqrt(K) * 1.5
    s = least_squares(fun, z0, jac=jac, method='trf', max_nfev=400, xtol=1e-12, ftol=1e-12)
    c = 2 * s.cost
    best = c if best is None else min(best, c)
    print(mode[-30:], 'K', K, 'sym' if sym else 'gen', 'restart', rs, 'res %.3g' % c, 'max %.2g' % np.abs(s.fun).max(), flush=True)
    if c < 1e-12: break
