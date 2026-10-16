"""Cutting-plane search for a per-vertex certificate of Conjecture Theta at the complementary gauge.
Vertex functional (t fixed):  Psi_c(A) = tr K(K-t)^2 + sum_e [ (3/2) e^{sig(a_e)} u_e^2 - (1/2) ab(ab-t)^2 + eta_c(a_e) ]
eta_c(x) = sum_k c_k (2x-1)^{2k+1}  (odd about 1/2, cancels over each edge since b = 1-a).
Goal: max_A Psi_c(A) <= n(1-t)^2/16 for all rank-n/2 projectors A (n in NS)."""
import numpy as np, sys, json
from scipy.optimize import minimize, linprog
t = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5
NS = [int(v) for v in (sys.argv[2] if len(sys.argv) > 2 else '4').split(',')]
DEG = int(sys.argv[3]) if len(sys.argv) > 3 else 7
CB = float(sys.argv[5]) if len(sys.argv) > 5 else 1.0
SIG = float(sys.argv[4]) if len(sys.argv) > 4 else 0.0   # sigma(x) = SIG*(2x-1)
rng = np.random.default_rng(0)
def proj(z, n):
    k = n // 2; Y = (z[:n*k] + 1j * z[n*k:]).reshape(n, k); Q, _ = np.linalg.qr(Y); return Q @ Q.conj().T
def parts(P):
    n = P.shape[0]; a = np.clip(np.real(np.diag(P)), 0, 1); b = 1 - a; K = P @ np.diag(b) @ P
    u = np.real(np.diag(K)) - a**2 * b; Kt = K - t * np.eye(n)
    base = np.real(np.trace(K @ Kt @ Kt)) + np.sum(1.5 * np.exp(SIG * (2*a-1)) * u**2 - 0.5 * a * b * (a * b - t)**2) - n * (1 - t)**2 / 16
    feats = np.array([np.sum((2 * a - 1) ** (2 * k + 1)) for k in range(DEG)])
    return base, feats
def val(P, c):
    b, f = parts(P); return (b + f @ c) / P.shape[0]
pool = []; c = np.zeros(DEG)
for rnd in range(int(__import__("os").environ.get("ROUNDS", "25"))):
    worst = -1; new = []
    for n in NS:
        for s in range(12):
            z0 = rng.standard_normal(n * n) if (s < 8 or not pool) else None
            if z0 is None:
                cand = [p for p in pool if p.shape[0] == n]
                if not cand: z0 = rng.standard_normal(n * n)
                else:
                    P0 = cand[rng.integers(len(cand))]; w, U = np.linalg.eigh(P0); Y = U[:, n//2:]
                    Y = Y + 0.05 * (rng.standard_normal(Y.shape) + 1j * rng.standard_normal(Y.shape)); z0 = np.concatenate([Y.real.ravel(), Y.imag.ravel()])
            r = minimize(lambda z: -val(proj(z, n), c), z0, method='BFGS', options=dict(maxiter=400))
            P = proj(r.x, n); v = val(P, c); new.append(P); worst = max(worst, v)
    pool += new
    # LP: minimise z s.t. (base_i + f_i c)/n_i <= z ; |c_k| <= 1
    B = np.array([parts(P)[0] / P.shape[0] for P in pool]); Fm = np.array([parts(P)[1] / P.shape[0] for P in pool])
    A_ub = np.hstack([Fm, -np.ones((len(pool), 1))]); res = linprog(np.r_[np.zeros(DEG), 1], A_ub=A_ub, b_ub=-B, bounds=[(-CB, CB)] * DEG + [(None, None)])
    c = res.x[:DEG]
    print(json.dumps(dict(round=rnd, worst_before_LP=round(worst, 6), lp_value=round(res.x[-1], 6), pool=len(pool), c=[round(v, 4) for v in c])), flush=True)
np.save(f'split_t{t}_n{"-".join(map(str,NS))}_sig{SIG}.npy', c)
vals = sorted(((val(P, c), np.round(np.sort(np.real(np.diag(P))), 3).tolist()) for P in pool), key=lambda z: -z[0])
for v, a in vals[:6]: print('worst', round(v, 7), a)
