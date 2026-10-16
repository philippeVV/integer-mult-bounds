"""expXH part 1b: fixed window h <= 9 with multiplicity, signature (8, rho).
Ports x = (v, z): v in E8, z in Z^rho, form v.v' - 4 z.z' (i.e. E8 + <-4>^rho = E8 + (2Z)^rho(-1)), norm 2:
v.v = 2 + 4|z|^2.  Labels (x mod 2) live in E8/2E8 (+1), so h <= 9 automatically; ports in a root class with
multiplicity (any pairwise-orthogonal family).  Valid = pairwise ip in {0,+-1}.  Decoder W/2 = sum of 8 E8-coordinate
terms (standard coordinates, half-integers allowed) + rho z-terms (weight -4).
usage: e8pad.py rho maxz2 cpsat_seconds"""
import sys, itertools, random
import numpy as np
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from lib import evaluate, max_clique

rho, maxz2, secs = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3])
# E8 doubled coordinates y = 2v
def e8_norm(N):  # vectors v with v.v = N, as doubled integer coords
    out = []
    T = 4 * N
    for par in (0, 1):
        vals = [t for t in range(-7, 8) if t % 2 == par and t * t <= T]
        cur = []
        def rec(i, q):
            if i == 8:
                if q == T and sum(cur) % 4 == 0:
                    out.append(cur.copy())
                return
            rem = 8 - i
            if par == 1 and q + rem > T: return
            for t in vals:
                if q + t * t <= T:
                    cur.append(t); rec(i + 1, q + t * t); cur.pop()
        rec(0, 0)
    return np.array(out, dtype=np.int64)

zs = {}
for z in itertools.product(range(-3, 4), repeat=rho):
    q = sum(t * t for t in z)
    if q <= maxz2:
        zs.setdefault(q, []).append(z)
P = []
for q, zl in sorted(zs.items()):
    Vv = e8_norm(2 + 4 * q)
    for z in zl:
        for y in Vv:
            vec = list(y) + list(z)
            # line representative: first nonzero of (z, y) positive
            key = list(z) + list(y)
            nz = next(t for t in key if t != 0)
            if nz > 0:
                P.append(vec)
P = np.array(P, dtype=np.int64)
print(f"rho={rho} maxz2={maxz2} candidates={len(P)}", flush=True)
# integer model for the Gram: doubled E8 coords with weight 1/4, z with weight -4  -> scale: W = (Y Y^T)/4 - 4 Z Z^T
Y, Z = P[:, :8], P[:, 8:]
W = (Y @ Y.T) // 4 - 4 * (Z @ Z.T)
assert np.all(np.diag(W) == 2)
adj = np.isin(W, (-1, 0, 1)); np.fill_diagonal(adj, False)
rng = random.Random(0)
deg = adj.sum(1)
best = []
for t in range(40):
    order = sorted(range(len(P)), key=lambda a: -deg[a] + (rng.random() * len(P) * 0.05 if t else 0))
    cl = []; ok = np.ones(len(P), dtype=bool)
    for a in order:
        if ok[a]: cl.append(a); ok &= adj[a]
    if len(cl) > len(best): best = cl
print('greedy', len(best), flush=True)
from lib import clique_ls
ls = clique_ls(adj, iters=int(secs*1000), restarts=3, init=best)
print('local search', len(ls), flush=True)
if len(ls) > len(best): best = ls
opt = None
if secs > 0 and len(P) <= 2500:
    c2, opt, ub = max_clique(adj, seconds=secs, hint=best)
    print('cpsat', len(c2), 'optimal' if opt else f'ub={ub}', flush=True)
    if len(c2) > len(best): best = c2
S = sorted(best)
# evaluate with a scaled integer representation: use X = [Y/2 (as rationals) , Z] -> evaluate needs integer X & eps;
# use X = [Y, Z] with eps weights (1/4, -4): implement directly
from lib import rank_f2_rows, bits_of
X = P[S]; Wv = W[np.ix_(S, S)]; V = len(S)
B = ((Wv + 1) % 2).astype(np.uint8); np.fill_diagonal(B, 1)
h = rank_f2_rows(bits_of(B)); rk = np.linalg.matrix_rank(Wv.astype(float))
C = 0; det = []
for c in range(8 + rho):
    Sc = np.nonzero(X[:, c])[0]
    if len(Sc):
        r = rank_f2_rows(bits_of(B[np.ix_(Sc, Sc)])); det.append(r); C += r
nz = int((Z[S] != 0).any(1).sum())
print(f"V={V} (z!=0: {nz}) h={h} rank={rk} C={C} det={det} gamma={C / V:.3f} merit={C / V * (h - 1):.2f}")
