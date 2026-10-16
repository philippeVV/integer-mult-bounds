"""Probe the positive Hessian direction: near-reducible flat MDS codes on W3, leverages pushed to 1/2 + eps*delta (points),
1/2 - eps*delta (lines) by per-vertex scaling (changes the code, not a gauge). Compute F and Theta_3 at the comp gauge
(= gauge 0 here) and inf over gauges."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from comp import *
from local import local_moments
from scipy.optimize import minimize
name = sys.argv[1]; mix = float(sys.argv[2]); V, nP, E, g = load_any(name); S = Struct(V, nP, E)
rng = np.random.default_rng(int(sys.argv[3]) if len(sys.argv) > 3 else 0)
PAIRS = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]
N = np.zeros((len(V), E))
for i_, es_ in enumerate(V): N[i_, es_] = 1
_, sv, Vt = np.linalg.svd(N); B = Vt[np.sum(sv > 1e-9):].T
def Omat_from(Qs, side):
    O = np.zeros((E, E)); r = range(nP) if side == 0 else range(nP, len(V))
    for i in r:
        es = V[i]; P = Qs[i] @ Qs[i].conj().T; o = abs(P)**2; np.fill_diagonal(o, 0); O[np.ix_(es, es)] = o
    return O
def scale_to(q, target, it=600):
    d = np.zeros(q.shape[0])
    for _ in range(it):
        pass
    for _ in range(it):
        Y = orth(np.exp(d)[:, None] * q); a = np.sum(abs(Y)**2, 1); d += 0.5*(np.log(target) - np.log(a))
    return orth(np.exp(d)[:, None] * q)
# reducible pairing config maximising lambda (greedy as in hess_red)
def lam_of(Qs):
    OP, OL = Omat_from(Qs, 0), Omat_from(Qs, 1); M = -3/16*np.eye(E) + 0.5*(OP+OL) - 1.5*(OP@OL+OL@OP)
    w, U_ = np.linalg.eigh(B.T @ M @ B); return w[-1], B @ U_[:, -1]
def red_codes(ch):
    Qs = []
    for i, es in enumerate(V):
        Y = np.zeros((4, 2), complex)
        for k, (x, y) in enumerate(PAIRS[ch[i]]):
            ph = np.exp(2j*np.pi*rng.random()); Y[x, k] = 1/np.sqrt(2); Y[y, k] = ph/np.sqrt(2)
        G = rng.standard_normal((4, 2)) + 1j*rng.standard_normal((4, 2))
        Qs.append(scale_to(orth(Y + mix*G), 0.5*np.ones(4)))
    return Qs
def lam_ideal(ch):
    OP = np.zeros((E, E)); OL = np.zeros((E, E))
    for i_, es in enumerate(V):
        O = OP if i_ < nP else OL
        for (x, y) in PAIRS[ch[i_]]: O[es[x], es[y]] = O[es[y], es[x]] = 0.25
    M = -3/16*np.eye(E) + 0.5*(OP+OL) - 1.5*(OP@OL+OL@OP); return np.linalg.eigvalsh(B.T @ M @ B)[-1]
ch = rng.integers(3, size=len(V)); cur = lam_ideal(ch); improved = True
while improved:
    improved = False
    for i in rng.permutation(len(V)):
        for c in range(3):
            if c == ch[i]: continue
            old = ch[i]; ch[i] = c; v = lam_ideal(ch)
            if v > cur + 1e-12: cur = v; improved = True
            else: ch[i] = old
print('ideal reducible lambda', cur, flush=True)
Qs = red_codes(ch)
lmax, dvec = lam_of(Qs); dvec = dvec / abs(dvec).max()
print(f'{name} mix={mix}: flat near-reducible codes, lambda_max(M|C) = {lmax:.5f}; min 2x2 Pluecker |minor|^2 = {min(0.25-np.max(np.triu(abs(q@q.conj().T)**2,1)) for q in Qs):.4f}', flush=True)
for eps in [float(v) for v in os.environ.get("EPS", "0.02,0.05,0.1").split(",")]:
    Qe = [scale_to(q, 0.5 + (eps if i < nP else -eps)*dvec[es]) for i, (es, q) in enumerate(zip(V, Qs))]
    a, b = leverages(S, Qe); m = local_moments(S, Qe)[0]; F = 4*m[3] - 4*m[2] + m[1]
    th0 = theta(S, Qe, 3); xc, err = comp_gauge(S, Qe); thc = theta(S, Qe, 3, xc)
    st = [xc, np.zeros(E), flat_gauge(S, Qe, 0)[0], flat_gauge(S, Qe, 1)[0]] + [xc + 0.5*rng.standard_normal(E) for _ in range(4)]
    best = min((minimize(lambda x: gauge_value_grad(S, Qe, 3, x), s0, jac=True, method='L-BFGS-B', options=dict(maxiter=5000, gtol=1e-12)) for s0 in st), key=lambda r: r.fun)
    print(f'  eps={eps}: max|a+b-1| at gauge 0 = {abs(a+b-1).max():.1e}; F(0) = {F:.6f}; Theta3(0) = {th0:.6f}; Theta3(comp) = {thc:.6f}; inf_D Theta3 = {best.fun:.6f} (spread {np.ptp(best.x):.2f})', flush=True)
