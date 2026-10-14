"""EXACT probe of the simplest 2-adic gadget: T1 = (1/2) u u^T, T2 = -(1/2)... with u = w + 2x, v = w
(sum (uu^T - vv^T)/2 = 2 (w x^T + x w^T)/... ) -- we test both signs: (uu^T +/- vv^T)/2.
Gain = rank_F2(B o (Mbar_req (x) J)) - (rho(supp u) + rho(supp v)), Mbar = reduction of the sum on required entries.
Random w, x supports (sizes 1..8) over random plane choices; report max gain (planes) and (lines, sanity)."""
import sys, random
from fractions import Fraction as Fr
from base import *
rng = random.Random(int(sys.argv[1])); RC = [rchoices(q) for q in P]
def run(planes_mode, trials):
    best = -99
    for t in range(trials):
        r = [rng.choice(RC[i]) for i in range(n)]
        basis = [[P[a], r[a]] for a in range(n)] if planes_mode else [[P[a]] for a in range(n)]
        if planes_mode:
            F = free_matrix(r); req = lambda a, b: not F[a][b]
        else:
            req = lambda a, b: a == b or dot(P[a], P[b]) == 1
        W = rng.sample(range(n), rng.randint(1, 8)); X = rng.sample(range(n), rng.randint(1, 8))
        u = {a: 0 for a in set(W) | set(X)}
        for a in W: u[a] += 1
        for a in X: u[a] += 2
        u = {a: x for a, x in u.items() if x}
        v = {a: 1 for a in W}
        for sgn in (1, -1):
            M = {}
            for a in u:
                for b in u: M[(a, b)] = M.get((a, b), 0) + Fr(u[a] * u[b], 2)
            for a in v:
                for b in v: M[(a, b)] = M.get((a, b), 0) + sgn * Fr(v[a] * v[b], 2)
            rows = []
            ok = True
            for a in range(n):
                for x in basis[a]:
                    vv = 0; bit = 0
                    for b in range(n):
                        val = M.get((a, b), 0); on = False
                        if req(a, b):
                            if Fr(val).denominator != 1: ok = False
                            on = Fr(val).numerator % 2 == 1
                        for y in basis[b]:
                            if on and dot(x, y): vv |= 1 << bit
                            bit += 1
                    rows.append(vv)
            if not ok: continue
            cost = pair_rank([w for i in u for w in basis[i]], [w for i in u for w in basis[i]]) + pair_rank([w for i in v for w in basis[i]], [w for i in v for w in basis[i]])
            best = max(best, rank2(rows) - cost)
    return best
print('lines max gain', run(False, 300), flush=True)
print('planes max gain', run(True, 300), flush=True)
