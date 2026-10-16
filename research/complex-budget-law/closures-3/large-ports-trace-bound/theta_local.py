"""expQL Theorem B (PSD decoders, local theta bound), recomputed on three families.
rank N >= V^2 / sum_a (1 + theta_a), theta_a = Lovasz theta of the NON-free graph induced on a's free neighbourhood N(a).
Families: E8 P120 (free = non-orthogonal root lines; vertex-transitive, one vertex suffices),
PG(2,4) dot and PG(3,4) dot field lines (fam.py; one representative per (degree, local-degree) class).
Expected (expQL table): theta = 14 / 2 / 3, PSD rank >= 8.000 / 5.57 / 16.0, PSD super-frame factor <= 1.875 / 0.958 / 1.000 from the local bound
(expQL's 0.889 for PG(2,4) uses the global theta(F) = 6 instead, rank >= 6).
Needs numpy and cvxpy (Clarabel)."""
import os, sys, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, cvxpy as cp
from lattices import E8
from fam import build


def lovasz(edges, m):
    """theta(G) = max <J,X>, tr X = 1, X_ij = 0 on edges of G, X psd."""
    if m <= 1: return float(m)
    X = cp.Variable((m, m), symmetric=True)
    cons = [X >> 0, cp.trace(X) == 1]
    if edges:
        I = np.array([e[0] for e in edges]); J = np.array([e[1] for e in edges]); cons.append(X[I, J] == 0)
    p = cp.Problem(cp.Maximize(cp.sum(X)), cons); p.solve(solver='CLARABEL')
    return p.value


def e8_free():
    L = E8(); G = L @ L.T; free = G != 0; np.fill_diagonal(free, False); return free, 8


def pg_free(u, n, kind):
    F, B, P = build(u, n, kind)
    V = len(P); free = np.zeros((V, V), dtype=bool)
    for i, j in itertools.combinations(range(V), 2):
        free[i, j] = free[j, i] = B(P[i], P[j]) == 0
    return free, n


def bound(free, omega, transitive=False):
    V = free.shape[0]; deg = free.sum(1); cls = {}
    for a in range(V):
        nb = np.nonzero(free[a])[0]
        cls.setdefault((int(deg[a]), tuple(sorted(free[np.ix_(nb, nb)].sum(1).tolist()))), []).append(a)
    if transitive: cls = {'all': list(range(V))}
    tot = 0.0; worst = 0.0
    for mem in cls.values():
        a = mem[0]; nb = np.nonzero(free[a])[0]; sub = free[np.ix_(nb, nb)]; m = len(nb)
        E = [(i, j) for i in range(m) for j in range(i + 1, m) if not sub[i, j]]   # non-free pairs inside N(a)
        th = lovasz(E, m); tot += len(mem) * (1 + th); worst = max(worst, th)
    rk = V * V / tot
    return V, worst, rk, (V / omega) / rk


ok = True
for name, (free, om), tr, exp in [('E8 P120', e8_free(), True, (14, 8.0, 1.875)),
                                  ('PG(2,4) dot', pg_free(2, 3, 'dot'), False, (2, 5.57, 0.958)),
                                  ('PG(3,4) dot', pg_free(2, 4, 'dot'), False, (3, 16.0, 1.0))]:
    V, th, rk, fac = bound(free, om, tr)
    print(f"{name:12s} V={V:4d} omega={om}: max theta_a = {th:.4f}  PSD rank >= {rk:.3f}  V/omega = {V/om:.3f}  PSD factor <= {fac:.4f}")
    ok &= abs(th - exp[0]) < 1e-3 and abs(rk - exp[1]) < 0.01 and abs(fac - exp[2]) < 2e-3
print('PASS_EXPQL_THETA' if ok else 'FAIL_EXPQL_THETA')
