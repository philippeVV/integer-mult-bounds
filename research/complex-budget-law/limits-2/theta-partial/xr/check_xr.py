"""expXR checker: Theorem L (second-order Conjecture Theta on the whole two-sided-flat manifold, any degrees, girth >= 8).
Prints PASS_EXPXR if all checks pass."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
import sys, os, numpy as np
from model import *
import xg_core as Dc
from gauge_hess import Gam as Gam_xk
ok = True
def chk(name, cond, info=''):
    global ok; ok &= bool(cond); print(('PASS ' if cond else 'FAIL ') + name, info, flush=True)
def inert(M, tol=1e-9):
    ev = np.linalg.eigvalsh((M + M.T) / 2); s = abs(ev).max()
    return (int((ev > tol * s).sum()), int((ev < -tol * s).sum()), int((abs(ev) <= tol * s).sum()))
def H0(X, Y):
    I = np.eye(len(X)); return np.block([[4 * Y, 3 * (I + Y @ X)], [3 * (I + X @ Y), 4 * X]])
def flatten(q, it=4000):
    d = np.zeros(q.shape[0])
    for _ in range(it):
        Yq = orth(np.exp(d)[:, None] * q); a = np.sum(abs(Yq) ** 2, 1); d += 0.5 * (np.log(.5) - np.log(a))
    return orth(np.exp(d)[:, None] * q)
def XY_codes(V, nP, E, Qs):
    X = np.zeros((E, E)); Y = np.zeros((E, E))
    for i, (es, q) in enumerate(zip(V, Qs)):
        P = q @ q.conj().T; O = 4 * abs(P) ** 2; np.fill_diagonal(O, 0); (X if i < nP else Y)[np.ix_(es, es)] = O
    return X, Y
rng = np.random.default_rng(7)

# L1: second-order model vs exact Theta_3 (point codes scaled by e^u, line codes by e^v), degree 4 and degree 6
for name in ['W3', 'W5']:
    V, nP, E, g = load_any(name); S = Struct(V, nP, E)
    Q0 = [flatten(q) for q in Dc.rand_codes(V, rng)]; X, Y = XY_codes(V, nP, E, Q0); K = Kmat(X, Y)
    th0 = theta(S, Q0, 3, np.zeros(E)); u, v = rng.standard_normal(E), rng.standard_normal(E); r = []
    for tau in [4e-3, 2e-3]:
        Qs = [Q0[i] * np.exp(tau * (u if i < nP else v)[es])[:, None] for i, es in enumerate(V)]
        z = tau * np.r_[u, v]; r.append((theta(S, Qs, 3, np.zeros(E)) - 0.125) / (4 / E * z @ K @ z))
    chk(f'L1 {name}: Theta3(flat)=1/8 and (Theta3-1/8)/model -> 1', abs(th0 - .125) < 1e-12 and abs(r[1] - 1) < 2 * abs(r[0] - 1) + 1e-9 and abs(r[1] - 1) < 0.05, f'th0-1/8={th0-.125:.1e} ratios={np.round(r, 4)}')
    Kxx = xdelta(K)[0]; chk(f'L1b {name}: diagonal block of K = expXK Gamma', np.allclose(Kxx, Gam_xk(X, Y), atol=1e-12))

# L2: Lemma A. H0 singular iff T=[[0,X],[Y,0]] has eigenvalue (-2 +- i sqrt5)/3 (needs unimodular non-real eig of XY)
lam0 = (-2 + 1j * np.sqrt(5)) / 3
chk('L2a: roots of 3l^2+4l+3 are unimodular, cos(2theta)=-1/9', abs(3 * lam0**2 + 4 * lam0 + 3) < 1e-14 and abs(abs(lam0) - 1) < 1e-14 and abs((lam0**2).real + 1/9) < 1e-14)
phi = np.arccos(-1 / 9) / 2; R = np.array([[np.cos(phi), -np.sin(phi)], [np.sin(phi), np.cos(phi)]])
Xs = np.diag([1., -1.]); Ys = R @ np.diag([1., -1.]) @ R.T
chk('L2b: reflections at the critical angle make H0 singular (sharpness)', abs(np.linalg.det(H0(Xs, Ys))) < 1e-12, f'det={np.linalg.det(H0(Xs, Ys)):.1e}')
worst = 1e9; inert_ok = True
for t in range(400):
    n = rng.integers(2, 8)
    def sc(neg_ok):
        Q = np.linalg.qr(rng.standard_normal((n, n)))[0]; e = rng.uniform(-1 if neg_ok else -0.999, 1, n)
        if rng.random() < 0.5: e[rng.integers(n)] = 1.0      # eigenvalue +1 allowed
        return (Q * e) @ Q.T
    X_, Y_ = sc(False), sc(True); ev = np.linalg.eigvalsh(H0(X_, Y_)); worst = min(worst, abs(ev).min())
    inert_ok &= ((ev > 0).sum() == n and (ev < 0).sum() == n)
chk('L2c: X > -I, |X|,|Y| <= 1 random: H0 nonsingular with inertia (n,n)', inert_ok and worst > 0, f'min |eig| {worst:.2e}')

# L3: H0 (1,1) = 10 (1,1); n_+(H0) = E; n_+(K) <= E-1 on flat codes of W3, W4M, W5 (complex, real, Fourier)
for name in ['W3', 'W4M', 'W5']:
    V, nP, E, g = load_any(name)
    for kind in ['rand', 'real', 'fourier']:
        Qs = Dc.fourier_codes(V, rng) if kind == 'fourier' else [flatten(q) for q in Dc.rand_codes(V, rng, cplx=(kind == 'rand'))]
        X, Y = XY_codes(V, nP, E, Qs); H = H0(X, Y); one = np.ones(2 * E)
        iH = inert(H); iK = inert(Kmat(X, Y))
        spec = np.linalg.eigvalsh(X); 
        chk(f'L3 {name} {kind}: |X|<=1, H0 1=10 1, n+(H0)<=E, n+(K)<=E-1', spec[0] > -1 - 1e-9 and spec[-1] < 1 + 1e-9 and np.allclose(H @ one, 10 * one) and iH[0] <= E and iK[0] <= E - 1,
            f'In(H0)={iH} In(K)={iK} E-1={E-1}')

# L4: corollary: Gamma > 0 on 1^perp  =>  Schur complement <= 0 (SIC on W3, random flat on W5) ; strict off the trivial kernel
for name, kind in [('W3', 'sic'), ('W3', 'dir'), ('W5', 'rand')]:
    V, nP, E, g = load_any(name)
    if name == 'W3':
        W = np.full((len(V), 3), 1 / 3) if kind == 'sic' else rng.dirichlet([6] * 3, size=len(V)); X, Y = XYmats(V, nP, E, W)
    else:
        X, Y = XY_codes(V, nP, E, [flatten(q) for q in Dc.rand_codes(V, rng)])
    K = Kmat(X, Y); ev, es_, Sm = schur(K); iG = inert(xdelta(K)[0]); nneg = int((es_ < -1e-9 * abs(es_).max()).sum())
    chk(f'L4 {name} {kind}: Gamma>0 on 1^perp, Schur<=0, #neg = E-V+1', iG == (E - 1, 0, 1) and es_[-1] < 1e-12 and nneg == E - len(V) + 1,
        f'In(Gamma)={iG} Schur max={es_[-1]:.1e} #neg={nneg} (E-V+1={E-len(V)+1})')
print('PASS_EXPXR' if ok else 'FAIL_EXPXR')
