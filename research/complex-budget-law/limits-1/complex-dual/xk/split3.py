"""Local LP duality certificate for Conjecture Theta at the complementary gauge (t = 1/2).
At a+b = 1 (girth >= 8):  E F/4 = sum_v tr K_v (K_v - 1/2)^2 + sum_e [ 3 u_e w_e - e(a_e) ],  e(x) = x(1-x)(x(1-x)-1/2)^2,
with K_v, u_e (resp. w_e) functions of the gauged vertex code alone. Seek h(x,u) with
  (EDGE)   3uw - e(x) <= h(x,u) + h(1-x,w)             for x in [0,1], u, w in [0, x(1-x)]
  (VERTEX) tr K(K-1/2)^2 + sum_e h(a_e, u_e) <= n/64      for every rank-n/2 projector on C^n.
Then F <= 1/8 for every code on every girth-8 graph with vertex degrees in NS (at the comp gauge).
h(x,u) = sum c_ij s^i v^j, s = 2x-1, v = 8u-1."""
import numpy as np, sys, json, os
from scipy.optimize import minimize, linprog
NS = [int(v) for v in (sys.argv[1] if len(sys.argv) > 1 else '4').split(',')]
DI, DJ = (int(sys.argv[2]), int(sys.argv[3])) if len(sys.argv) > 3 else (6, 3)
CB = float(sys.argv[4]) if len(sys.argv) > 4 else 50
ROUNDS = int(os.environ.get('ROUNDS', '20'))
IJ = [(i, j) for i in range(DI + 1) for j in range(DJ + 1)]
def feat(x, u):
    s = 2 * np.asarray(x) - 1; v = 8 * np.asarray(u) - 1
    return np.stack([s**i * v**j for i, j in IJ], axis=-1)
e = lambda x: x * (1 - x) * (x * (1 - x) - 0.5) ** 2
# edge grid
xs = np.unique(np.r_[np.linspace(0, 1, 81), np.geomspace(1e-4, 0.05, 25), 1 - np.geomspace(1e-4, 0.05, 25)]); rows = []; rhs = []
for x in xs:
    m = x * (1 - x); us = np.linspace(0, m, 15)
    for u in us:
        for w in us:
            rows.append(-(feat(x, u) + feat(1 - x, w))); rhs.append(e(x) - 3 * u * w)
EA = np.array(rows); Eb = np.array(rhs)
def proj(z, n):
    k = n // 2; Y = (z[:n*k] + 1j * z[n*k:]).reshape(n, k); Q, _ = np.linalg.qr(Y); return Q @ Q.conj().T
def vparts(P):
    n = P.shape[0]; a = np.clip(np.real(np.diag(P)), 0, 1); b = 1 - a; K = P @ np.diag(b) @ P
    u = np.real(np.diag(K)) - a**2 * b; Kt = K - 0.5 * np.eye(n)
    return np.real(np.trace(K @ Kt @ Kt)) - n / 64, feat(a, np.clip(u, 0, None)).sum(0)
def vval(P, c):
    b, f = vparts(P); return (b + f @ c) / P.shape[0]
rng = np.random.default_rng(int(os.environ.get('SEED', '0'))); pool = []; c = np.zeros(len(IJ))
# start c at the AM-GM split: h = 1.5 u^2 - e(x)/2  (approximate by least squares on grid)
for rnd in range(ROUNDS):
    worst = -1
    for n in NS:
        for s in range(10):
            cand = [p for p in pool if p.shape[0] == n]
            if s >= 6 and cand:
                P0 = cand[rng.integers(len(cand))]; w_, U = np.linalg.eigh(P0); Y = U[:, n//2:] + 0.05 * (rng.standard_normal((n, n//2)) + 1j * rng.standard_normal((n, n//2)))
                z0 = np.concatenate([Y.real.ravel(), Y.imag.ravel()])
            else: z0 = rng.standard_normal(n * n)
            r = minimize(lambda z: -vval(proj(z, n), c), z0, method='BFGS', options=dict(maxiter=300))
            P = proj(r.x, n); pool.append(P); worst = max(worst, vval(P, c))
    VB = np.array([vparts(P)[0] / P.shape[0] for P in pool]); VF = np.array([vparts(P)[1] / P.shape[0] for P in pool])
    nv = len(IJ)
    A_ub = np.vstack([np.hstack([VF, -np.ones((len(pool), 1))]), np.hstack([EA, np.zeros((len(EA), 1))])])
    b_ub = np.r_[-VB, Eb]
    res = linprog(np.r_[np.zeros(nv), 1], A_ub=A_ub, b_ub=b_ub, bounds=[(-CB, CB)] * nv + [(None, None)], method='highs')
    if res.status != 0: print('LP status', res.status, res.message); break
    c = res.x[:nv]
    print(json.dumps(dict(round=rnd, worst_before=round(worst, 7), lp=round(res.x[-1], 8), pool=len(pool))), flush=True)
np.save(f'split3_n{"-".join(map(str, NS))}_{DI}_{DJ}.npy', c)
vals = sorted(((vval(P, c), np.round(np.sort(np.real(np.diag(P))), 3).tolist()) for P in pool), key=lambda z: -z[0])
for v, a in vals[:5]: print('worst', round(v, 8), a)
# edge check on a fine grid
xs = np.linspace(0, 1, 201); mx = -1
for x in xs:
    m = x * (1 - x); us = np.linspace(0, m, 41); U, Wg = np.meshgrid(us, us)
    viol = 3 * U * Wg - e(x) - feat(x, U) @ c - feat(1 - x, Wg) @ c; mx = max(mx, viol.max())
print('fine-grid edge max violation', mx)
