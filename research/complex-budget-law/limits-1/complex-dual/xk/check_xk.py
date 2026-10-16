"""expXK checker:  uv run -q --no-project --python 3.12 --with numpy --with scipy python -I -B check_xk.py   (prints PASS_EXPXK; ~2 min)
K1 sparse Theta_3 + gradients == dense expXG code (W3).
K2 tree expansion of m3 (Lemma 1) == exact trace, random codes + random gauges (W3, W5).
K3 Theorem 2: if both sides are flat at some gauge, m = (1/4, 3/16, 5/32) and Theta_3 = 1/8 exactly, for ANY codes (girth >= 8).
K4 Theorem 3: complementary gauge (a_e + b_e = 1) exists (convex capacity) and dim W is unchanged.
K5 Reduction identity at the comp gauge: E*F/4 = sum_v tr K_v(K_v-1/2)^2 + sum_e [3 u_e w_e - e(a_e)],  F = 4m3-4m2+m1 >= Theta_3 >= s.
K6 W5 adversary end point: inf_D Theta_3 = 0.1249992 < 1/8 (multistart, unconstrained), comp gauge value < 1/8.
K7 analytic gauge Hessian Gam at flat codes == finite differences; Fourier W3 has a gauge-descent direction.
K8/K9 ETF codes (SIC on W3, icosahedral on W5): Gam > 0 on 1-perp and Theta_3^* = 1/8 exactly (Conjecture Theta is tight).
K10 the complementary gauge alone fails: its Hessian M|C has a positive eigenvalue at a reducible configuration.
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'xg'))
import xg_core as D
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from comp import *
from local import local_moments
ok = True
def check(c, msg):
    global ok; print(('ok   ' if c else 'FAIL ') + msg, flush=True); ok &= bool(c)
V, nP, E, g = load_any('W3'); S = Struct(V, nP, E); rng = np.random.default_rng(7)
Q = D.rand_codes(V, rng); x = rng.standard_normal(E) * 0.5
t1, g1 = D.gauge_value_grad(V, nP, E, Q, 3, x); t2, g2 = gauge_value_grad(S, Q, 3, x)
check(abs(t1 - t2) < 1e-12 and abs(g1 - g2).max() < 1e-12, f'K1 sparse == dense on W3: Theta {t2:.10f}, grad err {abs(g1-g2).max():.1e}')
for name in ['W3', 'W5']:
    V, nP, E, g = load_any(name); S = Struct(V, nP, E)
    Q = D.rand_codes(V, rng); x = rng.standard_normal(E) * 0.7
    Ys = [orth(np.exp(x[es])[:, None] * q) for es, q in zip(V, Q)]
    me = theta_grad(S, Ys, 3)[1]; ml = local_moments(S, Q, x)[0]
    check(abs(me - ml).max() < 1e-12, f'K2 {name}: exact m1..m3 {np.round(me[1:], 6)} == tree expansion (err {abs(me-ml).max():.1e})')
    # K3: per-vertex Barthe flattening (each vertex code rescaled on its own coordinates) -> flat codes, non-Fourier
    Qf = []
    for es, q in zip(V, Q):
        d = np.zeros(len(es))
        for _ in range(3000):
            Y = orth(np.exp(d)[:, None] * q); a = np.sum(abs(Y) ** 2, 1); d += 0.5 * (np.log(0.5) - np.log(a))
        Qf.append(orth(np.exp(d)[:, None] * q))
    a, b = leverages(S, Qf); mf = local_moments(S, Qf)[0]
    check(abs(a - .5).max() < 1e-9 and abs(b - .5).max() < 1e-9 and abs(mf[1:] - [1/4, 3/16, 5/32]).max() < 1e-9 and abs(theta(S, Qf, 3) - 1/8) < 1e-9,
          f'K3 {name}: two-sided flat random (non-Fourier) codes: m = {np.round(mf[1:], 8)}, Theta_3 = {theta(S, Qf, 3):.10f} = 1/8')
    xc, err = comp_gauge(S, Q); a, b = leverages(S, Q, xc)
    n0 = D.nullity(V, nP, E, Q)[0] if E < 400 else None; n1 = D.nullity(V, nP, E, [orth(np.exp(xc[es])[:, None] * q) for es, q in zip(V, Q)])[0] if E < 400 else None
    check(err < 1e-6 and n0 == n1, f'K4 {name}: comp gauge found, max|a+b-1| = {err:.1e}, lev range [{a.min():.3f},{a.max():.3f}], dim W {n0} -> {n1}')
    m, a, b, u, w = local_moments(S, Q, xc); F = 4 * m[3] - 4 * m[2] + m[1]
    tk = 0.0
    for i, (es, q) in enumerate(zip(V, Q)):
        Y = orth(np.exp(xc[es])[:, None] * q); P = Y @ Y.conj().T; d = 1 - np.real(np.diag(P)); K = P @ np.diag(d) @ P
        tk += np.real(np.trace(K @ (K - .5 * np.eye(len(es))) @ (K - .5 * np.eye(len(es)))))
    ee = a * b * (a * b - .5) ** 2; rhs = tk + np.sum(3 * u * w - ee)
    check(abs(E * F / 4 - rhs) < 1e-7 and F >= theta(S, Q, 3, xc) - 1e-12, f'K5 {name}: E F/4 = {E*F/4:.8f} = local sum {rhs:.8f}; F = {F:.6f} >= Theta_3(comp) = {theta(S, Q, 3, xc):.6f}')
V, nP, E, g = load_any('W5'); S = Struct(V, nP, E)
f = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'adv_W5_box3_s1_adv_W5_box3_s0_rand.npy')
Q = [orth(np.asarray(q, complex)) for q in np.load(f, allow_pickle=True)]
xc, err = comp_gauge(S, Q); tc = theta(S, Q, 3, xc)
from scipy.optimize import minimize
r = minimize(lambda x: gauge_value_grad(S, Q, 3, x), xc, jac=True, method='L-BFGS-B', options=dict(maxiter=20000, gtol=1e-13, ftol=1e-16))
check(r.fun < 0.125 and tc < 0.125, f'K6 W5 adversary end point: Theta_3(comp) = {tc:.9f}, inf_D = {r.fun:.9f} < 1/8')

# K7-K10: gauge Hessian at two-sided flat codes; ETF codes; comp-gauge failure near reducible codes
from gauge_hess import XY_from_codes, lam_min
from sic import sic_codes
V, nP, E, g = load_any('W3'); S = Struct(V, nP, E)
Qf = D.fourier_codes(V, np.random.default_rng(0)); X, Y = XY_from_codes(V, nP, E, Qf); l, v = lam_min(X, Y)
tau = 3e-3; tv = gauge_value_grad(S, Qf, 3, tau * v)[0] - 0.125
check(l < 0 and abs(tv - 4 / E * l * tau**2) < 1e-3 * abs(tv), f'K7 analytic gauge Hessian (W3 Fourier): lambda_min(Gam) = {l:.6f} < 0; predicted {4/E*l*tau**2:.4e} vs true {tv:.4e} -> gauge descends below 1/8')
Qs = sic_codes(V, np.random.default_rng(1)); X, Y = XY_from_codes(V, nP, E, Qs); l, _ = lam_min(X, Y)
r = minimize(lambda x: gauge_value_grad(S, Qs, 3, x), flat_gauge(S, Qs, 0)[0] + 0.3 * np.random.default_rng(2).standard_normal(E), jac=True, method='L-BFGS-B', options=dict(maxiter=5000, gtol=1e-13))
check(l > 0 and abs(r.fun - 0.125) < 1e-8, f'K8 W3 SIC (tetrahedral ETF) codes: lambda_min(Gam) = {l:.6f} > 0 (strict local min), Theta_3 inf from perturbed start = {r.fun:.10f} = 1/8 (attained)')
Vw, nPw, Ew, gw = load_any('W5')
Qi = [orth(np.asarray(q, complex)) for q in np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'init_W5_etf_n0_s0.npy'), allow_pickle=True)]
Xw, Yw = XY_from_codes(Vw, nPw, Ew, Qi); lw, _ = lam_min(Xw, Yw); Sw = Struct(Vw, nPw, Ew)
check(lw > 0 and abs(theta(Sw, Qi, 3) - 0.125) < 1e-12, f'K9 W5 icosahedral ETF(3,6) codes: Theta_3(0) = {theta(Sw, Qi, 3):.12f} = 1/8, lambda_min(Gam) = {lw:.6f} > 0')
# K10 comp-gauge counterexample direction at a reducible (pairing) configuration
PAIRS = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]
N_ = np.zeros((len(V), E))
for i_, es in enumerate(V): N_[i_, es] = 1
_, sv, Vt = np.linalg.svd(N_); Bc = Vt[np.sum(sv > 1e-9):].T
rng_ = np.random.default_rng(0); ch = rng_.integers(3, size=len(V))
def lam_red(ch):
    OP = np.zeros((E, E)); OL = np.zeros((E, E))
    for i_, es in enumerate(V):
        O = OP if i_ < nP else OL
        for (x_, y_) in PAIRS[ch[i_]]: O[es[x_], es[y_]] = O[es[y_], es[x_]] = 0.25
    M = -3/16*np.eye(E) + 0.5*(OP+OL) - 1.5*(OP@OL+OL@OP); return np.linalg.eigvalsh(Bc.T @ M @ Bc)[-1]
cur = lam_red(ch)
for i_ in rng_.permutation(len(V)):
    for c_ in range(3):
        old = ch[i_]; ch[i_] = c_; v_ = lam_red(ch)
        if v_ > cur: cur = v_
        else: ch[i_] = old
check(cur > 0, f'K10 comp-gauge Hessian M|C at a reducible degree-4 configuration on W3: lambda_max = {cur:.5f} > 0 (comp gauge alone cannot prove Theta)')

print('PASS_EXPXK' if ok else 'FAILED')
