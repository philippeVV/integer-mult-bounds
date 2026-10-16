"""R3 library. Lift with support S (rigid side X, R2 assumed): x_t = (r_t, 0) on T = V minus S, x_s = xi_s (c_s, 1) on S.
Row-wise: given X, the completion exists iff x_a notin span X_{orth(a)} for every line a; then y_a ranges over an affine space
of dim 8 - rank X_{orth(a)}.  Good line a: orth(a) minus S spans r_a^perp (rank 7).  At a good line a every pair s, s' in
S cap orth(a) must satisfy (c_s - c_s').r_a = 0 (whether a is in T or in S)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir
import sys, random
import numpy as np
pass   # expSM modules (e8, kfam, maxdef, cliques_lib, orbits, catalog, subsys) are copied into this directory
from e8 import *
from kfam import RT, dot, kernel, normline, CLASS, check, deficient, Fr
from maxdef import mincost
RI = np.array(ROOT, np.int64)             # doubled coords
RF = RI / 2.0
ORTHL = [[b for b in range(NL) if b != a and IP[a][b] == 0] for a in range(NL)]
P_ = 2147483647
def rank_modp(M):
    M = np.array(M, np.int64) % P_; r = 0; rows, cols = M.shape
    for c in range(cols):
        nz = np.nonzero(M[r:, c])[0]
        if len(nz) == 0: continue
        p = r + nz[0]; M[[r, p]] = M[[p, r]]
        inv = pow(int(M[r, c]), P_ - 2, P_); M[r] = (M[r] * inv) % P_
        f = M[:, c].copy(); f[r] = 0
        nzr = np.nonzero(f)[0]
        if len(nzr): M[nzr] = (M[nzr] - (f[nzr, None] * M[r][None, :]) % P_) % P_
        r += 1
        if r == rows: break
    return r
def rank_f(idx): return np.linalg.matrix_rank(RF[list(idx)], 1e-9) if len(idx) else 0
def bad_lines(S):
    S = set(S); return [a for a in range(NL) if rank_f([b for b in ORTHL[a] if b not in S]) < 7]
def csys(S, lines=None):
    """exact (mod-p upper bound on dim, matched by the concurrent lower bound) dimension of the c-system.
    lines: which lines a impose (c_s - c_s').r_a = 0 (default: all good lines).  Returns dim - (8 + |S|) = extra."""
    S = list(S); pos = {s: i for i, s in enumerate(S)}; n = 8 * len(S)
    if lines is None:
        bl = set(bad_lines(S)); lines = [a for a in range(NL) if a not in bl]
    rows = []
    for a in lines:
        ss = [s for s in S if IP[a][s] == 0 and s != a]
        for s in ss[1:]:
            row = np.zeros(n, np.int64); row[8 * pos[s]:8 * pos[s] + 8] = RI[a]; row[8 * pos[ss[0]]:8 * pos[ss[0]] + 8] -= RI[a]
            rows.append(row)
    rk = rank_modp(rows) if rows else 0
    return n - rk - 8 - len(S)
def build(S, c, xi, rnd, pick=None):
    """exact X from (S, c, xi); Y row-wise: random point of each affine solution space. Returns X, Y, freedom dims or None."""
    X = []
    for a in range(NL):
        X.append([xi[a] * v for v in c[a]] + [xi[a]] if a in S else list(RT[a]) + [Fr(0)])
    Y = []; fr = {}
    for a in range(NL):
        K = kernel([X[b] for b in ORTHL[a]])
        vals = [dot(X[a], k) for k in K]
        j = next((i for i, v in enumerate(vals) if v != 0), None)
        if j is None: return None
        y = [x / vals[j] for x in K[j]]
        hom = [[x - vals[i] / vals[j] * z for x, z in zip(K[i], K[j])] for i in range(len(K)) if i != j]
        fr[a] = len(hom)
        for h in hom:
            t = pick(a) if pick else Fr(rnd.randint(-9, 9), rnd.randint(1, 9))
            y = [x + t * z for x, z in zip(y, h)]
        Y.append(y)
    return X, Y, fr
def tov(g):
    v = np.array([float(x) for x in g]); return v / np.linalg.norm(v)
def gauge_floor(X, Y):
    dX, kx = deficient(X); dY, ky = deficient(Y)
    if kx > 1 or ky > 1: return None, dX, dY, (kx, ky)
    return mincost([(tov(g), v) for g, v in dX.items()], [(tov(g), v) for g, v in dY.items()])[0], dX, dY, (kx, ky)
def is_clique(S): return all(IP[a][b] != 0 for a in S for b in S if a != b)
def rank_S(S): return rank_f(S)
