"""kappa-family (canonical lift gauge), exact over Q.
X = [diag(al) R, ka],  Y = [diag(alp) R, kap]   (R = true root coords, r.r = 2)
supp ka = A, supp kap = B, no a in A, b in B with a perp b;  2 al_a alp_a + ka_a kap_a = 1 for every line a.
Every lift (T = V minus S rigid rank 8, S pairwise odd, roots of S independent) is GL9 x torus equivalent to a member
with A = S (proved in README sec 2; checked numerically by chkparam.py)."""
from fractions import Fraction as Fr
import random
from e8 import *
RT = [tuple(Fr(x, 2) for x in r) for r in ROOT]
def dot(u, v): return sum(x * y for x, y in zip(u, v))
CLASS = {h: [a for a in range(NL) if dotq(LAB[a], h)] for h in range(1, 511)}   # h = 511 (all 120 lines) excluded
def member(A, B, rnd, al=None, alp=None):
    """random exact member; al/alp: optional dicts forcing values (e.g. 0) on lines of A & B."""
    A = set(A); B = set(B)
    assert not any(IP[a][b] == 0 for a in A for b in B)
    al = dict(al or {}); alp = dict(alp or {})
    X = []; Y = []; P = {}
    def rq(): return Fr(rnd.choice([-1, 1]) * rnd.randint(1, 9), rnd.randint(1, 9))
    for a in range(NL):
        k = rq() if a in A else Fr(0); kp = rq() if a in B else Fr(0)
        if a in A and a in B:
            x = al.get(a, rq())
            if a in alp: y = alp[a]
            else: y = rq()
            if 2 * x * y == 1: y += 1
            # rescale kap to satisfy the diagonal: k kp = 1 - 2 x y
            kp = (1 - 2 * x * y) / k
            assert kp != 0
        else:
            x = Fr(1) if a not in al else al[a]          # torus normalisation
            y = (1 - k * kp) / (2 * x)
        X.append([x * c for c in RT[a]] + [k]); Y.append([y * c for c in RT[a]] + [kp])
        P[a] = (x, k, y, kp)
    return X, Y, P
def check(X, Y):
    for a in range(NL):
        for b in range(NL):
            if a == b or IP[a][b] == 0:
                v = dot(X[a], Y[b])
                if v != (1 if a == b else 0): return False
    return True
def kernel(rows, n=9):
    """exact right kernel basis of the matrix with given rows (Fractions)."""
    M = [list(r) for r in rows]; piv = []; r = 0
    for c in range(n):
        p = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if p is None: continue
        M[r], M[p] = M[p], M[r]; pv = M[r][c]; M[r] = [x / pv for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c]; M[i] = [x - f * y for x, y in zip(M[i], M[r])]
        piv.append(c); r += 1
        if r == len(M): break
    free = [c for c in range(n) if c not in piv]; K = []
    for fc in free:
        v = [Fr(0)] * n; v[fc] = Fr(1)
        for i, pc in enumerate(piv): v[pc] = -M[i][fc]
        K.append(v)
    return K
def normline(v):
    i = next(j for j, x in enumerate(v) if x != 0); return tuple(x / v[i] for x in v)
def deficient(X):
    """exact: dict line(normalised tuple) -> list of h with X_{O_h} g = 0; error if some kernel has dim > 1."""
    out = {}; kd = 0
    for h, O in CLASS.items():
        K = kernel([X[a] for a in O])
        kd = max(kd, len(K))
        if len(K) == 1: out.setdefault(normline(K[0]), []).append(h)
        elif len(K) > 1: out.setdefault(('dim', len(K)), []).append(h)
    return out, kd
