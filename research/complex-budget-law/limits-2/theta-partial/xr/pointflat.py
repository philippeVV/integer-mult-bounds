"""Exact formula at the point-flat gauge (all point codes flat, line codes arbitrary, girth >= 8, any degrees):
 E dm1 = 0,  E dm2 = q/4,  E dm3 = 3q/8 + Cc,   q = eps^T X eps,
 Cc = sum_p T_p(eps) - 1/8 sum_e eps_e^3 - 3/8 sum_e (X eps)_e eps_e^2,   T_p = sum_ijk Pi_ij Pi_jk Pi_ki eps_i eps_j eps_k,
 eps = line leverages - 1/2, X = 4|Pi_P|^{o2} - I.  Criterion: Theta3 <= 1/8  <=>  q/8 + Cc <= q^2/(2E)."""
import sys, numpy as np
from model import *
from check_xr_util import flatten, XY_codes
import xg_core as Dc
from xk_sparse import theta_grad
def exact(S, Qs):
    V, nP, E = S.V, S.nP, S.E
    b = np.zeros(E); a = np.zeros(E)
    for i, (es, q) in enumerate(zip(V, Qs)):
        (a if i < nP else b)[es] = np.sum(abs(q) ** 2, 1)
    eps = b - 0.5; X, _ = XY_codes(V, nP, E, Qs); qd = eps @ X @ eps
    T = 0.0
    for i in range(nP):
        es = V[i]; P = Qs[i] @ Qs[i].conj().T; PE = P * eps[es][None, :]; T += np.real(np.trace(PE @ PE @ PE))
    Cc = T - np.sum(eps ** 3) / 8 - 3 / 8 * np.sum((X @ eps) * eps ** 2)
    return a, eps, qd, Cc
name = sys.argv[1]; V, nP, E, g = load_any(name); S = Struct(V, nP, E); rng = np.random.default_rng(3)
for trial in range(4):
    Qs = [flatten(q) if i < nP else q for i, q in enumerate(Dc.rand_codes(V, rng))]
    if trial >= 2:   # only a point-disjoint set of lines non-flat (greedy), the rest flat
        used = set(); nonflat = []
        for i in rng.permutation(range(nP, len(V))):
            pts = {j for j in range(nP) if set(V[j]) & set(V[i])}
            if not (pts & used): used |= pts; nonflat.append(i)
        Qs = [Qs[i] if (i < nP or i in nonflat) else flatten(Qs[i]) for i in range(len(V))]
    a, eps, qd, Cc = exact(S, Qs)
    th, m, _ = theta_grad(S, Qs, 3)
    pred_m = np.array([1, 0.25, 3 / 16 + qd / (4 * E), 5 / 32 + (3 * qd / 8 + Cc) / E])
    crit = qd / 8 + Cc - qd ** 2 / (2 * E)
    print(f'trial {trial}: point lev dev {abs(a-.5).max():.1e}, |eps|max {abs(eps).max():.2f}, moments err {abs(m - pred_m).max():.1e}, '
          f'Theta3 = {th:.10f}, criterion q/8+C-q^2/2E = {crit:.3e} (sign should match Theta3-1/8 = {th-0.125:.3e})' + (f' [#nonflat lines {len(nonflat)}]' if trial >= 2 else ''))
