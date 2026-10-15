"""Numerical masked min-rank (Levenberg-Marquardt, floating point, uncertified) for F8-line ports and paired planes.
usage: c_lm.py f8:<trial>|paired K restarts [sym]"""
import sys, json, numpy as np
from scipy.optimize import least_squares
from scipy.sparse import coo_matrix
from common import *
import f8
mode, K, R = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]); sym = len(sys.argv) > 4
P = load_ports(); n = len(P)
if mode.startswith('f8:'):
    _, free = f8.build(1, int(mode[3:])); F = np.array(free)
else:
    cert = json.load(open(MAIN + '/research/paired-plane-port-matching-seed6c563a/CERTIFICATE.json'))
    pl = {}
    for i, j in cert['matching']: pl[i] = pl[j] = (P[i], P[j])
    F = np.array([[a != b and all(dot(x, y) == 0 for x in pl[a] for y in pl[b]) for b in range(n)] for a in range(n)])
rows = [(a, a) for a in range(n)] + [(a, b) for a in range(n) for b in range(n) if a != b and not F[a, b] and (not sym or a < b)]
ra = np.array([p[0] for p in rows]); rb = np.array([p[1] for p in rows]); tgt = (ra == rb).astype(float)
def unpack(z):
    if sym: V = z.reshape(n, K); return V, V
    return z[:n * K].reshape(n, K), z[n * K:].reshape(n, K)
def fun(z):
    X, Y = unpack(z); return (X[ra] * Y[rb]).sum(1) - tgt
def jac(z):
    X, Y = unpack(z); m = len(rows); I = np.repeat(np.arange(m), K)
    c1 = ra[:, None] * K + np.arange(K)
    c2 = (rb[:, None] * K + np.arange(K)) if sym else (n * K + rb[:, None] * K + np.arange(K))
    return coo_matrix((np.concatenate([Y[rb].ravel(), X[ra].ravel()]), (np.concatenate([I, I]), np.concatenate([c1.ravel(), c2.ravel()]))), shape=(m, z.size)).tocsr()
best = 1e9
for rs in range(R):
    rng = np.random.default_rng(1000 * K + rs)
    z0 = rng.standard_normal((1 if sym else 2) * n * K) / np.sqrt(K) * 1.5
    s = least_squares(fun, z0, jac=jac, method='trf', max_nfev=300, xtol=1e-12, ftol=1e-12)
    best = min(best, 2 * s.cost)
    if best < 1e-12: break
print(mode, 'K', K, 'sym' if sym else 'gen', 'best residual %.3g' % best, 'free pairs', int(F.sum()) // 2, flush=True)
