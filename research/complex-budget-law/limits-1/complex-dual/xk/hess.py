"""Second-order (Hessian) form of the comp-gauge functional at two-sided flat codes.
Exact identity (girth >= 8, comp gauge, delta = a - 1/2 on point side, b - 1/2 = -delta):
  E (F - 1/8)/4 = Q(delta) + O(|delta|^3),  Q = delta^T M delta,
  M = -(3/16) I + (1/2)(O_P + O_L) - (3/2)(O_P O_L + O_L O_P),   O_v = |A_v|^{o2} - diag(a^2) (off-diagonal part).
delta ranges over C = {star sums zero at every point and line}. Report lambda_max(M | C)."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'xg'))
import xg_core as D
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from xk_sparse import load_any, orth
def flatten(q):
    d = np.zeros(q.shape[0])
    for _ in range(2000):
        Y = orth(np.exp(d)[:, None] * q); a = np.sum(abs(Y)**2, 1); d += 0.5*(np.log(.5) - np.log(a))
    return orth(np.exp(d)[:, None] * q)
def lam(V, nP, E, Qs):
    OP = np.zeros((E, E)); OL = np.zeros((E, E)); Cn = []
    for i, (es, q) in enumerate(zip(V, Qs)):
        P = q @ q.conj().T; O = abs(P)**2; np.fill_diagonal(O, 0)
        (OP if i < nP else OL)[np.ix_(es, es)] = O
        r = np.zeros(E); r[es] = 1; Cn.append(r)
    M = -3/16*np.eye(E) + 0.5*(OP+OL) - 1.5*(OP@OL + OL@OP)
    N = np.array(Cn); U, s, Vt = np.linalg.svd(N); rk = np.sum(s > 1e-9); B = Vt[rk:].T   # basis of C
    w = np.linalg.eigvalsh(B.T @ M @ B); return w[-1], B.shape[1]
if __name__ == '__main__':
    rng = np.random.default_rng(0)
    for name in sys.argv[1:]:
        V, nP, E, g = load_any(name)
        for kind in ['fourier', 'rand-flat', 'rand-flat', 'rand-flat']:
            Qs = D.fourier_codes(V, rng) if kind == 'fourier' else [flatten(q) for q in D.rand_codes(V, rng)]
            l, dim = lam(V, nP, E, Qs); print(name, kind, 'dim C =', dim, ' lambda_max(M|C) =', round(l, 6), flush=True)
