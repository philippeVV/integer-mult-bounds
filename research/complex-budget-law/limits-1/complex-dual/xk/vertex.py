"""Per-vertex inequality test. For a rank-n/2 projector A on C^n: a = diag A, b = 1 - a, K = A D_b A,
u_e = K_ee - a_e^2 b_e.   Psi_t(A) = tr K (K-t)^2 + sum_e [ (3/2) u_e^2 - (1/2) a_e b_e (a_e b_e - t)^2 ].
Claim V: Psi_t(A) <= n (1-t)^2 / 16 (t = 1/2: n/64), equality at flat-leverage projectors."""
import numpy as np, sys
from scipy.optimize import minimize
def proj(z, n):
    k = n // 2; Y = (z[:n*k] + 1j * z[n*k:]).reshape(n, k); Q, _ = np.linalg.qr(Y); return Q @ Q.conj().T
def psi(P, t=0.5, lam=None):
    n = P.shape[0]; a = np.real(np.diag(P)); b = 1 - a; K = P @ np.diag(b) @ P
    u = np.real(np.diag(K)) - a**2 * b; Kt = K - t * np.eye(n)
    return np.real(np.trace(K @ Kt @ Kt)) + np.sum(1.5 * u**2 - 0.5 * a * b * (a * b - t)**2)
if __name__ == '__main__':
    t = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5
    rng = np.random.default_rng(0)
    for n in [2, 4, 6, 8]:
        best = -1
        for s in range(30 if n < 8 else 10):
            z0 = rng.standard_normal(n * n)
            r = minimize(lambda z: -psi(proj(z, n), t), z0, method='BFGS')
            if -r.fun > best: best = -r.fun; Pb = proj(r.x, n)
        a = np.real(np.diag(Pb))
        print(f'n={n} t={t}: max Psi = {best:.6f}  target n(1-t)^2/16 = {n*(1-t)**2/16:.6f}  ratio {best/(n*(1-t)**2/16):.4f}  lev at max {np.round(np.sort(a),3)}', flush=True)
