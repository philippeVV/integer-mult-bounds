"""expXR part 2 checker: linear-rule fallback refuted; Danskin derivative of Phi; precise ray monotonicity. Prints PASS_EXPXR2."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
import sys, numpy as np
from xr_tools import *
from xk_sparse import theta_grad
import xg_core as Dc
ok = True
def chk(name, cond, info=''):
    global ok; ok &= bool(cond); print(('PASS ' if cond else 'FAIL ') + name, info, flush=True)
# L5: the linear Schur gauge rule exceeds 1/8 on a Gamma > 0 base at |delta|_rms = 0.6, while the true inf over gauges stays below 1/8
S, Q0, W, G, C, D, lmin, rng = setup('W3', '4', 5); E = S.E
L = -np.linalg.pinv(G, rcond=1e-10) @ C
d = np.load(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), 'linsph_W3_4_5_r0.6.npy'))
Ys = [orth(np.exp((L @ d + (d if i >= S.nP else 0))[es])[:, None] * q) for i, (es, q) in enumerate(zip(S.V, Q0))]
PP, PL = Dc.projectors(S.V, S.nP, E, Ys); th_dense = Dc.markov_ref(Dc.moments(PP, PL, 3), 3)
th_sp = fuv(S, Q0, L @ d, L @ d + d)[0]; inf = inner(S, Q0, d, [L @ d, -d / 2, np.zeros(E)])[0]
chk('L5: Gamma>0 base, linear rule gives Theta3 > 1/8 (dense, independent), true inf < 1/8', lmin > 0 and th_dense > 0.125 + 1e-4 and abs(th_dense - th_sp) < 1e-10 and inf < 0.125 - 1e-3,
    f'lambda_min(Gamma)={lmin:.2e} |d|_rms={np.linalg.norm(d)/np.sqrt(E):.2f} lin={th_dense:.8f} inf={inf:.8f}')
# L6: Danskin: Phi'(1) = <g_v, delta> at the precise gauge optimum (finite differences of precisely solved Phi)
S, Q0, W, G, C, D, lmin, rng = setup('W3', 'sic', 0); E = S.E; L = -np.linalg.pinv(G, rcond=1e-10) @ C
d = rng.standard_normal(E); d *= 0.5 * np.sqrt(E) / np.linalg.norm(d)
th, x, gu, gv = inner(S, Q0, d, [L @ d, -d / 2]); h = 1e-4
p = inner(S, Q0, (1 + h) * d, [x])[0]; m_ = inner(S, Q0, (1 - h) * d, [x])[0]
R = gu @ x + gv @ (x + d)
chk('L6: Danskin  Phi\'(1) = <g_v,delta> = radial derivative <g_u,u>+<g_v,v>, and < 0', abs((p - m_) / (2 * h) - gv @ d) < 1e-6 * max(1, abs(gv @ d)) and abs(R - gv @ d) < 1e-8 and R < 0,
    f'FD={(p-m_)/(2*h):.8e} <g_v,d>={gv@d:.8e} R={R:.8e} |grad_x|={np.abs(gu+gv).max():.1e}')
# L7: precise ray (forward/backward continuation) through a Gamma-boundary base is monotone and <= 1/8
S, Q0, W, G, C, D, lmin, rng = setup('W3', 'bdry0.97', 21); E = S.E; L = -np.linalg.pinv(G, rcond=1e-10) @ C
d = rng.standard_normal(E); d *= 0.4 * np.sqrt(E) / np.linalg.norm(d)
ts = np.linspace(0.1, 1.5, 8); vals = np.full(len(ts), np.inf)
for sweep in [range(len(ts)), range(len(ts) - 1, -1, -1)]:
    xp = None
    for i in sweep:
        st = [L @ (ts[i] * d), -ts[i] * d / 2] + ([xp] if xp is not None else [])
        th, xp, _, _ = inner(S, Q0, ts[i] * d, st); vals[i] = min(vals[i], th)
chk('L7: precise ray through a near-Gamma-boundary base is monotone and below 1/8', all(np.diff(vals) <= 1e-10) and vals.max() < 0.125,
    f'lambda_min={lmin:.1e} Phi={np.round(vals, 7)}')
# L8: Theorem P (exact formula at the point-flat gauge, any degrees):  E dm1 = 0, E dm2 = q/4, E dm3 = 3q/8 + C,
#     q/8 + C = sum_p [ 1/4 sum_{i<j} X_ij e_i e_j + 6 sum_{i<j<k} Re(Pi_ij Pi_jk Pi_ki) e_i e_j e_k ],  Theta3 <= 1/8 <=> q/8 + C <= q^2/(2E)
from check_xr_util import flatten, XY_codes
from itertools import combinations
for name in ['W3', 'W5']:
    V, nP, E, g = load_any(name); S = Struct(V, nP, E); rng = np.random.default_rng(11)
    for trial in range(3):
        Qs = [flatten(q) if i < nP else q for i, q in enumerate(Dc.rand_codes(V, rng))]
        if trial == 2:   # point-disjoint non-flat lines only
            used = set(); keep = []
            for i in rng.permutation(range(nP, len(V))):
                pts = {j for j in range(nP) if set(V[j]) & set(V[i])}
                if not (pts & used): used |= pts; keep.append(i)
            Qs = [Qs[i] if (i < nP or i in keep) else flatten(Qs[i]) for i in range(len(V))]
        b = np.zeros(E)
        for i in range(nP, len(V)): b[V[i]] = np.sum(abs(Qs[i]) ** 2, 1)
        eps = b - 0.5; X, _ = XY_codes(V, nP, E, Qs); q = eps @ X @ eps; loc = 0.0
        for i in range(nP):
            es = V[i]; P = Qs[i] @ Qs[i].conj().T; e = eps[es]; n = len(es)
            loc += 0.25 * sum(X[es[a], es[c]] * e[a] * e[c] for a, c in combinations(range(n), 2))
            loc += 6 * sum(np.real(P[a, c] * P[c, k] * P[k, a]) * e[a] * e[c] * e[k] for a, c, k in combinations(range(n), 3))
        th, m, _ = theta_grad(S, Qs, 3)
        Cc = loc - q / 8; pred = np.array([1, .25, 3/16 + q / (4 * E), 5/32 + (3 * q / 8 + Cc) / E])
        crit = loc - q * q / (2 * E)
        good = abs(m - pred).max() < 1e-12 and (np.sign(crit) == np.sign(th - .125) or (abs(crit) < 1e-12 and abs(th - .125) < 1e-12))
        if trial == 2: good = good and abs(th - .125) < 1e-12
        chk(f'L8 {name} trial {trial}: point-flat exact moments + criterion sign' + (' + point-disjoint class gives exactly 1/8' if trial == 2 else ''), good,
            f'moment err {abs(m-pred).max():.1e}, crit {crit:.3e}, Theta3-1/8 {th-.125:.3e}')
print('PASS_EXPXR2' if ok else 'FAIL_EXPXR2')
