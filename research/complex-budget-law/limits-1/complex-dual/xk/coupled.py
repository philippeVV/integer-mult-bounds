"""Homogeneous coupled test: all vertices use the same projector P (n=4), with the line's coordinate pi(i) on colour i,
a o pi = 1 - a.  Exact comp-gauge value per (point,line) pair of E*F_{1/2}/4 = sum_v tr K(K-1/2)^2 - sum_e ab(ab-1/2)^2 + 3 sum_e u w.
Target E/32 -> n/32 per pair. Maximise over P with a penalty enforcing a o pi = 1 - a."""
import numpy as np
from scipy.optimize import minimize
t = 0.5
def proj(z, n):
    k = n // 2; Y = (z[:n*k] + 1j * z[n*k:]).reshape(n, k); Q, _ = np.linalg.qr(Y); return Q @ Q.conj().T
def pieces(P):
    n = P.shape[0]; a = np.real(np.diag(P)); b = 1 - a; K = P @ np.diag(b) @ P
    u = np.real(np.diag(K)) - a**2 * b; Kt = K - t * np.eye(n)
    return a, b, u, np.real(np.trace(K @ Kt @ Kt))
def G(P, pi):
    a, b, u, tk = pieces(P); n = len(a)
    return (2 * tk + np.sum(-a*b*(a*b-t)**2 + 3*u*u[pi]) - n/32) / n, np.sum((a[pi] - (1 - a))**2)
rng = np.random.default_rng(0)
for n, pi in [(4, np.array([2, 3, 0, 1])), (4, np.array([1, 0, 3, 2])), (6, np.array([3, 4, 5, 0, 1, 2]))]:
    best = -1
    for s in range(40):
        r = minimize(lambda z: -G(proj(z, n), pi)[0] + 200 * G(proj(z, n), pi)[1], rng.standard_normal(n*n), method='BFGS')
        P = proj(r.x, n); g, pen = G(P, pi)
        if pen < 1e-6 and g > best: best = g; Pb = P
    a, b, u, tk = pieces(Pb)
    print(n, pi, 'max coupled excess per coord', best, 'lev', np.round(a, 4), 'u', np.round(u, 4))
