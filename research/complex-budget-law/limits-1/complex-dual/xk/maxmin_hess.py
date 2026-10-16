"""Second-order max-min test of Conjecture Theta at a two-sided flat code C0 (Theta_3(C0, gauge 0) = 1/8 exactly).
Variables: code tangent W_v in C^{(n-k) x k} (Y_v -> Y_v + Yperp_v W_v), gauge y in R^E.
f(W, y) = Theta_3(codes(W), gauge y).  H = Hessian at 0 (finite differences of analytic gradients).
Conjecture-consistent iff  max_W min_y  of the quadratic model <= 0:
  if H_yy has a negative eigenvalue (non-trivial gauge direction) the gauge can push Theta below 1/8 locally;
  else Schur complement S = H_WW - H_Wy H_yy^+ H_yW must be <= 0 (else codes with Theta* > 1/8 exist nearby)."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'xg'))
import xg_core as D
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from xk_sparse import *

def setup(V, Qs):
    perp = []
    for Y in Qs:
        n, k = Y.shape; U, s, Vh = np.linalg.svd(Y, full_matrices=True); perp.append(U[:, k:])
    return perp

def grad_all(S, Qs, perp, w, y):
    """w: real vector of all W_v (re, im), y: gauge. returns f, grad_w, grad_y."""
    off = 0; Xs = []; Ws = []
    for Y, Pp in zip(Qs, perp):
        n, k = Y.shape; m = (n - k) * k
        W = (w[off:off + m] + 1j * w[off + m: off + 2 * m]).reshape(n - k, k); off += 2 * m; Ws.append(W)
    Ys = []; Rs = []; Dd = []
    for es, Y, Pp, W in zip(S.V, Qs, perp, Ws):
        d = np.exp(y[es]); X = d[:, None] * (Y + Pp @ W); Qx, R = np.linalg.qr(X); Ys.append(Qx); Rs.append(R); Dd.append(d)
    th, Zs, m = code_grad(S, Ys, 3)
    gw = []
    for Z, R, d, Pp in zip(Zs, Rs, Dd, perp):
        G = Pp.conj().T @ (d[:, None] * Z) @ np.linalg.inv(R).conj().T
        gw.append(np.r_[G.real.ravel(), G.imag.ravel()])
    gw = np.concatenate(gw)
    # gauge gradient at the perturbed code
    Qw = [Y + Pp @ W for Y, Pp, W in zip(Qs, perp, Ws)]
    _, gy = gauge_value_grad(S, Qw, 3, y)
    return th, gw, gy

def hessian(S, Qs, h=1e-5):
    perp = setup(S.V, Qs); nw = sum(2 * (Y.shape[0] - Y.shape[1]) * Y.shape[1] for Y in Qs); N = nw + S.E
    H = np.zeros((N, N)); z = np.zeros(N)
    for i in range(N):
        zp = z.copy(); zp[i] = h; zm = z.copy(); zm[i] = -h
        _, a1, b1 = grad_all(S, Qs, perp, zp[:nw], zp[nw:]); _, a2, b2 = grad_all(S, Qs, perp, zm[:nw], zm[nw:])
        H[:, i] = (np.r_[a1, b1] - np.r_[a2, b2]) / (2 * h)
    H = (H + H.T) / 2
    f0, gw0, gy0 = grad_all(S, Qs, perp, z[:nw], z[nw:])
    return f0, np.r_[gw0, gy0], H, nw

def test(S, Qs, label):
    f0, g0, H, nw = hessian(S, Qs)
    Hww, Hwy, Hyy = H[:nw, :nw], H[:nw, nw:], H[nw:, nw:]
    ey = np.linalg.eigvalsh(Hyy)
    scale = abs(H).max()
    tol = 1e-6 * scale
    negy = ey[0] < -tol
    if negy:
        verdict = 'gauge has a descent direction (H_yy min eig %.3e): local Theta* < 1/8 near C0 except the flat manifold' % ey[0]
        Smax = None
    else:
        Hp = np.linalg.pinv(Hyy, rcond=1e-8)
        Sm = Hww - Hwy @ Hp @ Hwy.T; Smax = np.linalg.eigvalsh((Sm + Sm.T) / 2)[-1]
        verdict = 'H_yy >= 0 (min %.2e); Schur complement max eig %.3e (%s)' % (ey[0], Smax, 'OK <= 0' if Smax <= tol else 'POSITIVE: counterexample direction')
    print(f'{label}: f0 = {f0:.10f}, |grad| = {abs(g0).max():.1e}, max|H| = {scale:.2e}; {verdict}', flush=True)
    return ey, Smax

if __name__ == '__main__':
    name = sys.argv[1]; V, nP, E, g = load_any(name); S = Struct(V, nP, E); rng = np.random.default_rng(0)
    def flatten(q):
        d = np.zeros(q.shape[0])
        for _ in range(3000):
            Y = orth(np.exp(d)[:, None] * q); a = np.sum(abs(Y)**2, 1); d += 0.5*(np.log(.5) - np.log(a))
        return orth(np.exp(d)[:, None] * q)
    test(S, D.fourier_codes(V, rng), name + ' Fourier')
    for s in range(2):
        test(S, [flatten(q) for q in D.rand_codes(V, rng)], name + f' random flat #{s}')
